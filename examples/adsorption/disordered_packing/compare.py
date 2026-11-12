"""Verify and compare ordered and disordered eight-cell tracer archives."""
import argparse
from itertools import combinations
from pathlib import Path
import shutil

import numpy as np

from examples.adsorption.repeated_geometry.predict import LIMITS
from examples.adsorption.transport_accuracy.run import validate_record
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT
from pore2process_demo.provenance import source_manifest, verified_source_path

from .extend import combined_observations
from .mesh import choose_cut, compact, load_field
from .native import HERE, ORDERED, ORDERED_SHA, OUTPUT
from .tracer import observe_cases, tracer_plan

SETUPS = (("ordered",64),("disordered",48),("disordered",64))
SIZES = ("d800","d450")


def profile_check(first,second,limit):
    errors = {key:float(max(abs(np.array(first[key])[2:]-np.array(second[key])[2:])))
              for key in LIMITS}
    return dict(limit=limit,errors=errors,passed=all(v<=limit for v in errors.values()))


def difference(ordered,disordered):
    return {key:(np.array(disordered[key])-np.array(ordered[key])).tolist()
            for key in (*LIMITS,"fast_conditional_capture")}


def evaluate(base):
    if sha(ORDERED/"summary.json") != ORDERED_SHA:
        raise ValueError("Pinned ordered eight-cell reference changed")
    previous = validate_record(ORDERED/"summary.json")["analysis"]["cases"]
    frozen_ordered = read_json(ORDERED/"frozen_prediction.json")["predictions"]
    common_cut = choose_cut(base)
    data = {}
    hydraulic = {}
    for label,n in SETUPS:
        name = f"{label}_N{n}"
        directory = base/name
        face = read_json(directory/"faces/summary.json")
        _, field_coords, _, native = load_field(base,label,n)
        if native != face["info"]:
            raise ValueError("Native field or geometry changed after conservative projection")
        cut = 0 if label=="ordered" else common_cut*(n//16)
        saved_plan = read_json(directory/"tracer_plan.json")
        new_plan = tracer_plan(base,label,n,cut)
        recorded_sources = saved_plan["source_sha256"]
        if any(sha(verified_source_path(name,digest))!=digest for name,digest in recorded_sources.items()):
            raise ValueError("A source used by the frozen tracer changed")
        # Catalogue metadata can be added after the tracer files are frozen.
        new_plan["source_sha256"] = recorded_sources
        if new_plan != saved_plan or compact(directory,cut,write=False) != saved_plan["compact"]:
            raise ValueError("Tracer field or matched Péclet plan differs on replay")
        replay = observe_cases(directory,saved_plan)
        archived = read_json(directory/"tracer_summary.json")
        if replay != archived["cases"] or sha(directory/"tracer_plan.json") != archived["tracer_plan_sha256"]:
            raise ValueError("Tracer evaluation differs on replay")
        data[name] = replay
        coords = np.fromfile(directory/"faces/coords.i32",dtype="<i4").reshape(-1,3)
        np.testing.assert_array_equal(coords,field_coords)
        sections = ((coords[:,0]-cut)%face["info"]["shape"][0])//n
        voids = np.bincount(sections)/n**3
        mask = np.zeros((len(voids),n,n,n),dtype=bool)
        mask[sections,(coords[:,0]-cut)%n,coords[:,1],coords[:,2]] = True
        minimum_section_difference = (min(float(np.mean(a!=b)) for a,b in combinations(mask,2))
                                      if label=="disordered" else 0.)
        if label=="disordered" and minimum_section_difference<=0.:
            raise ValueError("The resolved disordered geometry repeats an axial section")
        hydraulic[name] = dict(info=face["info"],balance=face["checks"],compact=saved_plan["compact"],
                               section_porosity=voids.tolist(),
                               minimum_section_voxel_difference=minimum_section_difference,
                               tracer_plan_sha256=archived["tracer_plan_sha256"])
    ordered_volume = hydraulic["ordered_N64"]["info"]["voxel_solid_volume_m3"]*8
    disordered_volume = hydraulic["disordered_N64"]["info"]["voxel_solid_volume_m3"]
    voxel_difference = disordered_volume/ordered_volume-1
    if abs(voxel_difference)>.005:
        raise ValueError("Voxelized solids differ by more than 0.5%")
    time = {}
    time_combined = {}
    grid = {}
    deltas = {}
    transfer = {}
    for size in SIZES:
        p,s = size+"_primary",size+"_half_step"
        time[size] = {label:profile_check(data[label][p]["observed"],data[label][s]["observed"],.01)
                      for label in ("ordered_N64","disordered_N64")}
        time_combined[size] = {}
        for label in time[size]:
            folder = base/label
            saved = folder/(size+"_combined_observations.json")
            if saved.exists():
                combined = combined_observations(folder,size)
                if combined != read_json(saved):
                    raise ValueError("Additional time-step observations differ on replay")
                time_combined[size][label] = profile_check(combined[p],combined[s],.01)
        grid[size] = dict(disordered=profile_check(data["disordered_N48"][p]["observed"],
                                                   data["disordered_N64"][p]["observed"],.02),
                          ordered=profile_check(previous[p]["observed"],data["ordered_N64"][p]["observed"],.02))
        deltas[size] = {variant:difference(data["ordered_N64"][variant]["observed"],
                                           data["disordered_N64"][variant]["observed"])
                        for variant in (p,s)}
        transfer[size] = {
            label:{variant:profile_check(frozen_ordered[variant]["repeated"],
                                        data[label][variant]["observed"],.01)
                   for variant in (p,s)} for label in ("ordered_N64","disordered_N64")}
    gradients = [hydraulic[f"disordered_N{n}"]["info"]["raw"]["gradient_Pa_m"][0] for n in (48,64)]
    return dict(tracer=data,hydraulics=hydraulic,matched_grid_geometry_effect=deltas,
                time_step=time,time_step_combined=time_combined,spatial=grid,ordered_closure_transfer=transfer,
                axial_gradient_relative_N48_to_N64=gradients[0]/gradients[1]-1,
                common_disordered_cut_sixteenths=common_cut,
                voxel_solid_difference_relative=voxel_difference,
                all_recovery_passed=all(r["inert_recovery_passed"] for c in data.values() for r in c.values()),
                time_step_passed=all(time_combined[size].get(label,score)["passed"]
                                     for size,c in time.items() for label,score in c.items()),
                spatial_passed=all(r["passed"] for c in grid.values() for r in c.values()))


def seal(base):
    result = evaluate(base)
    manifest = source_manifest([HERE/name for name in
        ("packing.py","native.py","mesh.py","tracer.py","track.jl","balance.jl","compare.py","extend.py",
         "protocol.md","test_mesh.py","test_track.jl")]
        +list((HERE/"openlb").iterdir())+[HERE.parent/"repeated_geometry/track.jl"])
    for name in manifest:
        p = base/"source_snapshot"/name
        p.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(ROOT/name,p)
    summary = base/"summary.json"
    write_json(summary,dict(schema="pore2process.disordered-packing.v1",run_id=base.name,
        analysis=result, source_sha256=manifest, ordered_reference_sha256=ORDERED_SHA,
        artifacts_sha256={str(p.relative_to(base)):sha(p) for p in base.rglob("*")
                          if p.is_file() and p!=summary},
        carbon_mass_kg=read_json(base/"plan.json")["dry_carbon_mass_kg"],
        new_adsorption_runs=0,figure7_fits=0))
    write_json(OUTPUT/"latest.json",dict(summary=str(summary.relative_to(ROOT))))
    print(summary,flush=True)


def check(path):
    saved = validate_record(path)
    if evaluate(path.parent)!=saved["analysis"]:
        raise ValueError("Ordered-disordered comparison changed on replay")
    print(dict(sources=len(saved["source_sha256"]),artifacts=len(saved["artifacts_sha256"]),
               recovery=saved["analysis"]["all_recovery_passed"],
               time=saved["analysis"]["time_step_passed"],space=saved["analysis"]["spatial_passed"]))


if __name__=="__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("path",type=Path)
    p.add_argument("--check",action="store_true")
    a = p.parse_args()
    check(a.path.resolve()) if a.check else seal(a.path.resolve())
