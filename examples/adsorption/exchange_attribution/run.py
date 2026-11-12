"""Run and replay the fixed-inventory factorial comparison."""
import argparse
from pathlib import Path
import shutil

import numpy as np

from examples.adsorption.alhamed.study import compare, crossing
from examples.adsorption.exchange_column.closure import closure
from examples.adsorption.exchange_column.run import CASES, PIXELS, source_inputs, validate_solution
from examples.adsorption.transport_accuracy.run import validate_record
from pore2process_demo.bridge import runtime_versions
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT
from pore2process_demo.provenance import source_manifest

from . import kernel
from .markers import COLUMN_SHA, HERE, MARKERS, OUTPUT, inputs, marker_specs, paired_laws

ALLOCATIONS = {"void":1,"flow":0,"slow":2}
LAWS = ("none","synthetic","matched_Pe")
LEVELS = ((1,24,8.,.08),(2,48,4.,.04),(4,96,2.,.02))


def fractions(cl, power):
    w = np.array(cl["flow_fractions"])
    s = np.array(cl["normalized_transit_times"])
    mass = w*s**power
    return mass/mass.sum()


def receiving(counts):
    counts = np.array(counts,dtype=float)
    return counts.T / counts.sum(axis=0)[:,None]


def model(base, cl, points, sp):
    plan = read_json(base/"marker_plan.json")
    specs = marker_specs(base,points,sp)
    if (plan["cases"]!=specs or plan["protocol_sha256"]!=sha(HERE/"protocol.md")
            or plan["preparation_source_sha256"]!=sha(HERE/"markers.py")):
        raise ValueError("Marker request differs from the frozen protocol")
    laws = paired_laws(base,specs,cl)
    if laws!=read_json(base/"marker_laws.json"):
        raise ValueError("Frozen matched-Peclet exchange laws differ on replay")
    n = len(cl["flow_fractions"])
    mix = {"none":np.eye(n).tolist(),"synthetic":receiving(cl["transition_counts"]).tolist(),
           "matched_Pe":{size:receiving(case["primary"]["counts"]).tolist() for size,case in laws.items()}}
    return dict(closure=cl, allocations={k:fractions(cl,p).tolist() for k,p in ALLOCATIONS.items()},
                mix=mix,marker_laws=laws,column_reference_sha256=COLUMN_SHA)


def t50(r, feed):
    return crossing(r["time_s"],np.array(r["outlet_concentration_kg_m3"])/feed,.5)


def resolution(a,b,feed):
    t = np.array(a["time_s"])
    mask = t>=24.
    old, new = t50(a,feed),t50(b,feed)
    changes = dict(t50_s=abs(new-old),
        max_fraction_after_startup=float(max(abs(np.array(a["outlet_concentration_kg_m3"])[mask]-
            np.interp(t[mask],b["time_s"],b["outlet_concentration_kg_m3"])))/feed))
    return dict(changes=changes,passed=changes["t50_s"]<=30 and changes["max_fraction_after_startup"]<=.005)


def allocation_check(result, p, trace, model, allocation):
    validate_solution(result,p,trace)
    mass = trace["dry_inventory_kg"]
    expected = mass*np.array(model["allocations"][allocation])
    w = np.array(model["closure"]["flow_fractions"])
    liquid = trace["liquid_volume_m3"]*w*np.array(model["closure"]["normalized_transit_times"])
    np.testing.assert_allclose(result["class_carbon_mass_kg"],expected,rtol=1e-12,atol=0)
    np.testing.assert_allclose(result["class_liquid_volume_m3"],liquid,rtol=1e-12,atol=0)
    np.testing.assert_allclose(sum(result["class_carbon_mass_kg"]),mass,rtol=1e-12,atol=0)


def allocation_summary(model, trace, sp):
    w = np.array(model["closure"]["flow_fractions"])
    edge = np.r_[0.,np.cumsum(w)]
    w_liquid = w*np.array(model["closure"]["normalized_transit_times"])
    volume = trace["liquid_volume_m3"]*w_liquid
    solid = trace["dry_inventory_kg"]/sp["particle_density_kg_m3"]
    first = np.diff(np.minimum(edge,.1))/w
    last = np.diff(np.maximum(edge,.9))/w
    out = {}
    for name, allocation in model["allocations"].items():
        shares = np.array(allocation)
        void = volume/(volume+solid*shares)
        out[name] = dict(fastest_tenth_carbon_fraction=float(first@shares),
                         slowest_tenth_carbon_fraction=float(last@shares),
                         local_voidage_range=[float(min(void)),float(max(void))])
    return out


def score(base, frozen, expected_pixels):
    if sha(PIXELS)!=expected_pixels:
        raise ValueError("Figure 7 marker archive changed")
    output = {}
    for size,(series,_) in CASES.items():
        scenarios = {}
        for allocation in ALLOCATIONS:
            scenarios[allocation] = {}
            for law in LAWS:
                item = frozen["cases"][f"{size}_{allocation}_{law}"]
                file = ROOT/item["source"] if item["archived"] else base/item["source"]
                case = read_json(file)
                r = compare(case["result"],series)
                scenarios[allocation][law] = dict(t50_s=r["t50"]["predicted_s"],RMSE=r["RMSE"],
                    passes_midpoint=r["t50_screen_passed"],passes_RMSE=r["RMSE_screen_passed"],
                    conservation=case["result"]["normalized_mass_residual"])
        t = {a:{k:v["t50_s"] for k,v in law.items()} for a,law in scenarios.items()}
        obs = r["t50"]["observed_interpolated_s"]
        allocation_effect = {a:{k:t[a][k]-t["void"][k] for k in LAWS} for a in ("flow","slow")}
        law_effect = {a:{k:t[a][k]-t[a]["synthetic"] for k in ("none","matched_Pe")} for a in ALLOCATIONS}
        interaction = {a:{k:law_effect[a][k]-law_effect["void"][k] for k in ("none","matched_Pe")}
                       for a in ("flow","slow")}
        candidates = [(abs(value-obs),a,k,value-obs) for a,group in t.items() for k,value in group.items()]
        _, near_allocation, near_law, near_gap = min(candidates)
        output[size] = dict(observed_t50_s=obs,scenarios=scenarios,allocation_effect_s=allocation_effect,
                            law_effect_s=law_effect,interaction_s=interaction,
                            original_gap_s=t["void"]["synthetic"]-obs,
                            closest_midpoint=dict(allocation=near_allocation,law=near_law,gap_s=near_gap))
    return output


def run(base):
    baselines, points, sp = inputs()
    _, f, _, _ = source_inputs()
    cl = closure()
    physical = model(base,cl,points,sp)
    if (base/"frozen_predictions.json").exists():
        raise ValueError("Completed cohort exists. Use --check.")
    files = [HERE/name for name in ("run.py","markers.py","kernel.py","protocol.md")]
    files += [kernel.ORIGINAL,HERE.parent/"path_exchange/track.jl",HERE.parent/"exchange_column/run.py",
              HERE.parent/"exchange_column/closure.py",HERE.parent/"alhamed/study.py",
              kernel.JULIA/"FilmParticle.jl"]
    hashes = source_manifest(files)
    for name in hashes:
        target = base/"source_snapshot"/name
        target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(ROOT/name,target)
    generated = base/"AllocationColumn.jl"
    generated.write_text(kernel.source())
    verified = kernel.verify(generated)
    write_json(base/"model.json",physical)
    write_json(base/"allocation_summary.json",allocation_summary(physical,next(iter(points.values())),sp))
    frozen = dict(model_sha256=sha(base/"model.json"),cases={})
    print("Attribution simulations",base,verified,flush=True)
    for size in CASES:
        reference = read_json(baselines[size]["synthetic"])
        for allocation, carbon in physical["allocations"].items():
            for law in LAWS:
                key = f"{size}_{allocation}_{law}"
                if allocation=="void" and law in baselines[size]:
                    old = baselines[size][law]
                    frozen["cases"][key] = dict(source=str(old.relative_to(ROOT)),archived=True,sha256=sha(old))
                    continue
                matrix = physical["mix"][law]
                if law=="matched_Pe":
                    matrix = matrix[size]
                previous = None
                refinements = []
                for index,(subcells,radial,dt,startup) in enumerate(LEVELS):
                    p = {**reference["inputs"],"subcells":subcells,"radial_cells":radial,"dt_s":dt,"startup_dt_s":startup}
                    r = kernel.solve(generated,cl["flow_fractions"],cl["normalized_transit_times"],matrix,carbon,**p)
                    allocation_check(r,p,points[size],physical,allocation)
                    file = base/f"{key}_level{index}.json"
                    write_json(file,dict(inputs=p,trace=points[size],result=r,allocation=allocation,law=law))
                    print(key,index,"midpoint min",t50(r,p["inlet_concentration_kg_m3"])/60,"runtime",r["total_s"],flush=True)
                    if previous is not None:
                        refinements.append(resolution(previous,r,p["inlet_concentration_kg_m3"]))
                        if refinements[-1]["passed"]:
                            break
                    previous = r
                if not refinements[-1]["passed"]:
                    raise ValueError(f"Numerical refinement failed for {key}: {refinements}")
                frozen["cases"][key] = dict(source=file.name,archived=False,sha256=sha(file),refinements=refinements)
    write_json(base/"frozen_predictions.json",frozen)
    print("Predictions frozen before Figure 7 scoring",sha(base/"frozen_predictions.json"),flush=True)
    scored = score(base,frozen,f["pixel_csv_sha256"])
    record = dict(schema="pore2process.exchange-attribution.v1",run_id=base.name,
        source_sha256=hashes,artifacts_sha256={str(p.relative_to(base)):sha(p) for p in base.rglob("*") if p.is_file()},
        source_pixel_sha256=f["pixel_csv_sha256"],frozen_prediction_sha256=sha(base/"frozen_predictions.json"),
        allocation_checks=verified,analysis=scored,versions=runtime_versions(),new_column_fitted_parameters=0,
        new_OpenLB_runs=0,new_DWSIM_runs=0)
    path = base/"summary.json"
    write_json(path,record)
    write_json(OUTPUT/"latest.json",dict(summary=str(path.relative_to(ROOT))))
    print("Completed",path,flush=True)


def check(path):
    record = validate_record(path)
    _, points, sp = inputs()
    base = path.parent
    physical = model(base,closure(),points,sp)
    frozen = read_json(base/"frozen_predictions.json")
    if (read_json(base/"model.json")!=physical or (base/"AllocationColumn.jl").read_text()!=kernel.source()
            or sha(base/"frozen_predictions.json")!=record["frozen_prediction_sha256"]
            or sha(base/"model.json")!=frozen["model_sha256"]):
        raise ValueError("Frozen case setup differs on replay")
    for case in frozen["cases"].values():
        file = ROOT/case["source"] if case["archived"] else base/case["source"]
        if sha(file)!=case["sha256"]:
            raise ValueError("A frozen result differs on replay")
        if not case["archived"]:
            r = read_json(file)
            allocation_check(r["result"],r["inputs"],r["trace"],physical,r["allocation"])
            levels = [read_json(p)["result"] for p in sorted(base.glob(file.name.split("_level")[0]+"_level*.json"))]
            if [resolution(a,b,r["inputs"]["inlet_concentration_kg_m3"]) for a,b in zip(levels,levels[1:])]!=case["refinements"]:
                raise ValueError("Refinement comparison differs on replay")
    if score(base,frozen,record["source_pixel_sha256"])!=record["analysis"]:
        raise ValueError("Figure 7 sensitivity differs on replay")
    print(dict(sources=len(record["source_sha256"]),artifacts=len(record["artifacts_sha256"]),scenarios=len(frozen["cases"])))


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--prepared",type=Path)
    modes.add_argument("--check",type=Path)
    args = parser.parse_args()
    check(args.check.resolve()) if args.check else run(args.prepared.resolve())
