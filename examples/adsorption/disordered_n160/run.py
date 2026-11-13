"""Project, track, score and replay the N160 refinement."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import shutil

import numpy as np

from examples.adsorption.disordered_packing.compare import profile_check
from examples.adsorption.disordered_packing.extend import combined_observations
from examples.adsorption.disordered_packing.mesh import compact, load_field, shared_faces
from examples.adsorption.disordered_packing.tracer import observe_cases, tracer_plan
from examples.adsorption.disordered_n128.native import PARENT as N96
from examples.adsorption.disordered_refinement.projection import PARITY_LIMIT, project
from examples.adsorption.dispersion_holdout.run import execute
from examples.adsorption.transport_accuracy.run import validate_record
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT
from pore2process_demo.provenance import source_manifest, verified_source_path

from .native import GRID, HERE, OUTPUT, PACKING, PARENT, PARENT_GRID, PARENT_SHA, n96_archive, parent

SIZES = ("d800", "d450")
PARENT_FOLDER = PARENT/f"disordered_N{PARENT_GRID}"


def handoff(base):
    n96 = n96_archive()
    parity_path = "projection_N64_parity/summary.json"
    parity_file = N96/parity_path
    if sha(parity_file)!=n96["artifacts_sha256"][parity_path]:
        raise ValueError("Archived N64 projector verification changed")
    parity = read_json(parity_file)
    accelerator = HERE.parent/"disordered_refinement/balance.jl"
    if (parity["max_face_difference_relative"]>PARITY_LIMIT
            or parity["source_sha256"]!=sha(accelerator)):
        raise ValueError("Conservative projector requires the archived parity check")
    folder,coords,velocity,info = load_field(base,"disordered",GRID)
    directory = folder/"faces"
    directory.mkdir()
    faces = shared_faces(coords,velocity,info["shape"],info["raw"]["dx_m"])
    coords.astype("<i4",copy=False).tofile(directory/"coords.i32")
    for name,value in faces.items():
        value.tofile(directory/name)
    result,_ = project(directory,info,coords,faces)
    write_json(directory/"summary.json",result)


def cut_index(folder):
    offset = read_json(PARENT_FOLDER/"compact/summary.json")["cut_offset_m"]
    h = read_json(folder/"raw.json")["dx_m"]
    cut = round(offset/h)
    np.testing.assert_allclose(cut*h,offset,rtol=1e-13,atol=1e-15)
    return cut


def trace(base):
    folder = base/f"disordered_N{GRID}"
    if (folder/"tracer_summary.json").exists():
        raise ValueError("Completed tracer data already exist")
    cut = cut_index(folder)
    compact(folder,cut)
    plan = tracer_plan(base,"disordered",GRID,cut)
    path = folder/"tracer_plan.json"
    write_json(path,plan)
    print("N160 tracer inputs frozen",sha(path),flush=True)

    def run_case(item):
        name,case = item
        elapsed = execute(case["command"],folder/(name+".log"),timeout=3600)
        write_json(folder/(name+"_execution.json"),elapsed)
        print(name,elapsed,flush=True)

    jobs = sorted(plan["cases"].items(),key=lambda item:(item[0].endswith("_half_step"),item[0]))
    with ThreadPoolExecutor(max_workers=2) as pool:
        list(pool.map(run_case,jobs))
    values = observe_cases(folder,plan)
    write_json(folder/"tracer_summary.json",dict(cases=values,tracer_plan_sha256=sha(path)))
    print("N160 eight-section tracer complete",folder,flush=True)


def scores(n128,n160,folder):
    grid,time,extended = {},{},{}
    for size in SIZES:
        primary,half = size+"_primary",size+"_half_step"
        time[size] = profile_check(n160[primary]["observed"],n160[half]["observed"],.01)
        for key in (primary,half):
            np.testing.assert_allclose(n128[key]["Pe"],n160[key]["Pe"],rtol=1e-14,atol=0.)
            old,new = (data[key]["observed"] for data in (n128,n160))
            grid[key] = profile_check(old,new,.02)
            grid[key]["outlet_survival_difference"] = new["survival"][-1]-old["survival"][-1]
            grid[key]["outlet_fast_survival_difference"] = new["fast_survival"][-1]-old["fast_survival"][-1]
        additional = folder/(size+"_combined_observations.json")
        if additional.exists():
            both = combined_observations(folder,size)
            if both!=read_json(additional):
                raise ValueError("The additional tracer ensemble differs on replay")
            extended[size] = profile_check(both[primary],both[half],.01)
    return dict(grid=grid,time_step=time,time_step_combined=extended,
                primary_grid_passed=all(grid[size+"_primary"]["passed"] for size in SIZES),
                grid_passed=all(case["passed"] for case in grid.values()),
                time_step_passed=all(extended.get(size,case)["passed"] for size,case in time.items()),
                recovery_passed=all(case["inert_recovery_passed"] for case in n160.values()))


def evaluate(base):
    previous = parent()["analysis"]
    for name in ("plan.json","disordered_centres_m.txt"):
        if sha(base/name)!=sha(PARENT/name):
            raise ValueError("The geometry or carbon inventory differs from N128")
    folder,coords,_,native = load_field(base,"disordered",GRID)
    handoff_info = read_json(folder/"faces/summary.json")
    if native!=read_json(folder/"verified_field.json") or native!=handoff_info["info"]:
        raise ValueError("The accepted N160 field differs on replay")
    np.testing.assert_array_equal(np.fromfile(folder/"faces/coords.i32",dtype="<i4").reshape(-1,3),coords)
    voxel_delta = native["voxel_solid_volume_m3"]/previous["native"]["voxel_solid_volume_m3"]-1
    if abs(voxel_delta)>.005:
        raise ValueError("Voxel carbon volume changed by more than 0.5%")
    path = folder/"tracer_plan.json"
    frozen = read_json(path)
    cut = cut_index(folder)
    current = tracer_plan(base,"disordered",GRID,cut)
    for name,digest in frozen["source_sha256"].items():
        if sha(verified_source_path(name,digest))!=digest:
            raise ValueError("A frozen tracer source changed")
    current["source_sha256"] = frozen["source_sha256"]
    if current!=frozen or compact(folder,cut,write=False)!=frozen["compact"]:
        raise ValueError("The tracer plan or compact flow differs on replay")
    n160 = observe_cases(folder,frozen)
    if read_json(folder/"tracer_summary.json")!=dict(cases=n160,tracer_plan_sha256=sha(path)):
        raise ValueError("N160 tracer observations differ on replay")
    plan128 = read_json(PARENT_FOLDER/"tracer_plan.json")
    for key in ("marker_count","seed","stop_s"):
        if plan128[key]!=frozen[key]:
            raise ValueError(f"The N160 tracer setting differs from N128: {key}")
    for key,case in frozen["cases"].items():
        if case["dt_s"]!=plan128["cases"][key]["dt_s"]:
            raise ValueError("The grid comparison requires matched tracer time steps")
    n128 = previous["N128"]
    return dict(N128=n128,N160=n160,**scores(n128,n160,folder),native=native,
                projection=handoff_info["checks"],compact=frozen["compact"],
                tracer_plan_sha256=sha(path),voxel_solid_change_relative=voxel_delta,
                geometry_sha256=sha(base/"disordered_centres_m.txt"),
                dry_carbon_mass_kg=previous["dry_carbon_mass_kg"])


def seal(base):
    result = evaluate(base)
    files = [HERE/name for name in ("native.py","run.py","protocol.md")]
    files += [PACKING/name for name in ("openlb/capillary3d.cpp","openlb/config.mk","packing.py",
                                       "mesh.py","tracer.py","track.jl","extend.py","compare.py")]
    files += [HERE.parent/name for name in ("disordered_n128/native.py","disordered_refinement/native.py",
             "disordered_refinement/balance.jl","disordered_refinement/projection.py",
             "repeated_geometry/track.jl","repeated_geometry/predict.py",
             "repeated_geometry/diagnostics.py","carbon_contact/track_contacts.jl","path_exchange/track.jl")]
    files += [ROOT/"julia/adsorption/PorousFaceTransport.jl"]
    manifest = source_manifest(files)
    for name in manifest:
        target = base/"source_snapshot"/name
        target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(ROOT/name,target)
    summary = base/"summary.json"
    write_json(summary,dict(schema="pore2process.disordered-n160.v1",run_id=base.name,
        analysis=result,parent_summary_sha256=PARENT_SHA,source_sha256=manifest,
        artifacts_sha256={str(f.relative_to(base)):sha(f) for f in base.rglob("*") if f.is_file() and f!=summary},
        new_adsorption_runs=0,figure7_fits=0))
    write_json(OUTPUT/"latest.json",dict(summary=str(summary.relative_to(ROOT))))
    print(summary,flush=True)


def check(path):
    saved = validate_record(path)
    if evaluate(path.parent)!=saved["analysis"]:
        raise ValueError("The N160 comparison differs on replay")
    a = saved["analysis"]
    print(dict(sources=len(saved["source_sha256"]),artifacts=len(saved["artifacts_sha256"]),
               grid=a["grid_passed"],time=a["time_step_passed"],recovery=a["recovery_passed"]))


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("path",type=Path)
    parser.add_argument("action",choices=("handoff","trace","seal","check"))
    args = parser.parse_args()
    {"handoff":handoff,"trace":trace,"seal":seal,"check":check}[args.action](args.path.resolve())
