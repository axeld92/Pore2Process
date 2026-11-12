"""Run and replay the N64-to-N96 eight-section tracer comparison."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import shutil

import numpy as np

from examples.adsorption.disordered_packing.compare import profile_check
from examples.adsorption.disordered_packing.extend import combined_observations
from examples.adsorption.disordered_packing.mesh import compact, load_field
from examples.adsorption.disordered_packing.tracer import observe_cases, tracer_plan
from examples.adsorption.dispersion_holdout.run import execute
from examples.adsorption.transport_accuracy.run import validate_record
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT
from pore2process_demo.provenance import source_manifest, verified_source_path

from .native import GRID, HERE, OUTPUT, PACKING, PARENT, PARENT_SHA, parent
from .projection import handoff

SIZES = ("d800","d450")


def archived_cut(base):
    archived = read_json(PARENT/"disordered_N64/compact/summary.json")
    info = read_json(base/f"disordered_N{GRID}/faces/summary.json")["info"]
    h = info["raw"]["dx_m"]
    cut = round(archived["cut_offset_m"]/h)
    np.testing.assert_allclose(cut*h,archived["cut_offset_m"],rtol=1e-13,atol=1e-15)
    return cut


def trace(base):
    folder = base/f"disordered_N{GRID}"
    if (folder/"tracer_summary.json").exists():
        raise ValueError("Completed N96 tracer archive already exists")
    cut = archived_cut(base)
    compact(folder,cut)
    plan = tracer_plan(base,"disordered",GRID,cut)
    path = folder/"tracer_plan.json"
    write_json(path,plan)
    print("N96 matched tracer frozen",sha(path),flush=True)

    def run_case(item):
        name,case = item
        elapsed = execute(case["command"],folder/(name+".log"),timeout=3600)
        write_json(folder/(name+"_execution.json"),elapsed)
        print(name,elapsed,flush=True)

    with ThreadPoolExecutor(max_workers=2) as pool:
        jobs = sorted(plan["cases"].items(),key=lambda item:(item[0].endswith("_half_step"),item[0]))
        list(pool.map(run_case,jobs))
    values = observe_cases(folder,plan)
    write_json(folder/"tracer_summary.json",dict(cases=values,tracer_plan_sha256=sha(path)))
    print("N96 eight-section tracer complete",folder,flush=True)


def evaluate(base):
    older = parent()
    for name in ("plan.json","disordered_centres_m.txt","disordered_N64/compact/summary.json",
                 "disordered_N64/tracer_plan.json","disordered_N64/faces/summary.json"):
        if sha(PARENT/name)!=older["artifacts_sha256"][name]:
            raise ValueError(f"Pinned N64 comparison changed: {name}")
    for name in ("plan.json","disordered_centres_m.txt"):
        if sha(base/name)!=sha(PARENT/name):
            raise ValueError("The refinement geometry or inventory differs")
    folder,coords,_,native = load_field(base,"disordered",GRID)
    flow = read_json(folder/"faces/summary.json")
    if native!=flow["info"] or native!=read_json(folder/"verified_field.json"):
        raise ValueError("The accepted N96 field differs on replay")
    np.testing.assert_array_equal(np.fromfile(folder/"faces/coords.i32",dtype="<i4").reshape(-1,3),coords)
    cut = archived_cut(base)
    plan_path = folder/"tracer_plan.json"
    saved = read_json(plan_path)
    current = tracer_plan(base,"disordered",GRID,cut)
    for name,digest in saved["source_sha256"].items():
        if sha(verified_source_path(name,digest))!=digest:
            raise ValueError("A frozen tracer source changed")
    current["source_sha256"] = saved["source_sha256"]
    if current!=saved or compact(folder,cut,write=False)!=saved["compact"]:
        raise ValueError("The N96 tracer plan or compact flow differs")
    observed = observe_cases(folder,saved)
    if read_json(folder/"tracer_summary.json")!=dict(cases=observed,tracer_plan_sha256=sha(plan_path)):
        raise ValueError("N96 tracer observations differ on replay")
    n64 = older["analysis"]["tracer"]["disordered_N64"]
    plan64 = read_json(PARENT/"disordered_N64/tracer_plan.json")
    if (plan64["marker_count"],plan64["seed"],plan64["stop_s"])!=(saved["marker_count"],saved["seed"],saved["stop_s"]):
        raise ValueError("N96 inlet sampling or observation horizon differs from N64")
    voxel_change = native["voxel_solid_volume_m3"]/older["analysis"]["hydraulics"]["disordered_N64"]["info"]["voxel_solid_volume_m3"]-1
    if abs(voxel_change)>.005:
        raise ValueError("Voxel carbon volume changed by more than 0.5%")
    grid = {}
    time = {}
    extended = {}
    for size in SIZES:
        a,b = (observed[size+"_"+name]["observed"] for name in ("primary","half_step"))
        time[size] = profile_check(a,b,.01)
        for name in ("primary","half_step"):
            key = size+"_"+name
            coarse = n64[key]
            fine = observed[key]
            np.testing.assert_allclose(fine["Pe"],coarse["Pe"],rtol=1e-14,atol=0.)
            if saved["cases"][key]["dt_s"]!=plan64["cases"][key]["dt_s"]:
                raise ValueError("The grid comparison uses different tracer time steps")
            grid[key] = profile_check(coarse["observed"],fine["observed"],.02)
            grid[key]["outlet_survival_difference"] = fine["observed"]["survival"][-1]-coarse["observed"]["survival"][-1]
            grid[key]["outlet_fast_survival_difference"] = fine["observed"]["fast_survival"][-1]-coarse["observed"]["fast_survival"][-1]
        extra = folder/(size+"_combined_observations.json")
        if extra.exists():
            combined = combined_observations(folder,size)
            if combined!=read_json(extra):
                raise ValueError("The additional tracer ensemble differs")
            extended[size] = profile_check(combined[size+"_primary"],combined[size+"_half_step"],.01)
    return dict(geometry_sha256=sha(base/"disordered_centres_m.txt"),dry_carbon_mass_kg=older["carbon_mass_kg"],
                voxel_solid_change_relative=voxel_change,native=native,projection=flow["checks"],
                projection_parity=read_json(base/"projection_N64_parity/summary.json")["max_face_difference_relative"],
                compact=saved["compact"],N64=n64,N96=observed,grid=grid,time_step=time,time_step_combined=extended,
                n96_tracer_plan_sha256=sha(plan_path),
                primary_grid_passed=all(grid[size+"_primary"]["passed"] for size in SIZES),
                grid_passed=all(c["passed"] for c in grid.values()),
                recovery_passed=all(c["inert_recovery_passed"] for c in observed.values()),
                time_step_passed=all(extended.get(size,c)["passed"] for size,c in time.items()))


def seal(base):
    result = evaluate(base)
    files = [HERE/name for name in ("native.py","projection.py","balance.jl","run.py","protocol.md")]
    files += [PACKING/name for name in ("openlb/capillary3d.cpp","openlb/config.mk","packing.py",
                                       "mesh.py","tracer.py","track.jl","extend.py","compare.py")]
    files += [HERE.parent/name for name in ("repeated_geometry/track.jl","repeated_geometry/predict.py",
                                           "repeated_geometry/diagnostics.py","carbon_contact/track_contacts.jl",
                                           "path_exchange/track.jl")]
    files += [ROOT/"julia/adsorption/PorousFaceTransport.jl"]
    manifest = source_manifest(files)
    for name in manifest:
        p = base/"source_snapshot"/name
        p.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(ROOT/name,p)
    summary = base/"summary.json"
    write_json(summary,dict(schema="pore2process.disordered-n96.v1",run_id=base.name,
        analysis=result,parent_summary_sha256=PARENT_SHA,source_sha256=manifest,
        artifacts_sha256={str(p.relative_to(base)):sha(p) for p in base.rglob("*") if p.is_file() and p!=summary},
        new_adsorption_runs=0,figure7_fits=0))
    write_json(OUTPUT/"latest.json",dict(summary=str(summary.relative_to(ROOT))))
    print(summary,flush=True)


def check(path):
    record = validate_record(path)
    if evaluate(path.parent)!=record["analysis"]:
        raise ValueError("The N96 comparison differs on replay")
    print(dict(sources=len(record["source_sha256"]),artifacts=len(record["artifacts_sha256"]),
               grid=record["analysis"]["grid_passed"],time=record["analysis"]["time_step_passed"],
               recovery=record["analysis"]["recovery_passed"]))


if __name__=="__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("path",type=Path)
    p.add_argument("action",choices=("handoff","trace","seal","check"))
    a = p.parse_args()
    {"handoff":handoff,"trace":trace,"seal":seal,"check":check}[a.action](a.path.resolve())
