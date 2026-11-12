"""Run and replay the shared-carbon Figure 7 study."""
import argparse
from datetime import datetime, timezone
from functools import lru_cache
from pathlib import Path
import shutil
from uuid import uuid4

import numpy as np

from examples.adsorption.alhamed.study import compare
from examples.adsorption.contact_column.inputs import model as wall_model
from examples.adsorption.contact_column.run import (SHORT_S, crossings, join_difference, refine,
                                                     scoring_curve, startup_refinement)
from examples.adsorption.exchange_attribution.markers import inputs as baseline_inputs
from examples.adsorption.exchange_column.run import (CASES, JULIA, PIXELS, solver as original_solver,
                                                     source_inputs, validate_solution)
from examples.adsorption.film_particle.metrics import column_passed
from examples.adsorption.transport_accuracy.run import validate_record
from pore2process_demo.bridge import runtime, runtime_versions
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT
from pore2process_demo.provenance import source_manifest

HERE = Path(__file__).parent
OUTPUT = ROOT/"results/adsorption/shared_contact"
PREVIOUS = ROOT/"results/adsorption/contact_column/runs/20261112T112916Z-e6cbd3c6"
PREVIOUS_SHA = "cae97ed6130b5aeb23ecbc772bfa211a86604801e153b5b87059f2e7a83c80b8"
LEVELS = ((1,32,8.,.08),(2,64,4.,.04),(4,96,2.,.02),(6,128,1.5,.015),
          (12,128,1.,.01),(18,128,.75,.0075))
STARTUP = ((6,96,.015),(8,128,.01),(12,128,.0075))


@lru_cache(maxsize=1)
def solver():
    original_solver()
    jl = runtime()
    jl.include(str(JULIA/"SharedContactColumn.jl"))
    return jl.SharedContactColumn.solve_column


def solve(cl,mix,weights,p):
    jl = runtime()
    result = solver()(np.array(cl["flow_fractions"]),np.array(cl["normalized_transit_times"]),
                      np.array(mix),np.array(weights),**p)
    out = {}
    for field in jl.fieldnames(jl.typeof(result)):
        key = str(field)
        value = getattr(result,key)
        out[key] = np.array(value,dtype=float,copy=True).tolist() if bool(jl.isa(value,jl.AbstractArray)) else float(value)
    return out


def model():
    path = PREVIOUS/"summary.json"
    if sha(path)!=PREVIOUS_SHA:
        raise ValueError("Separate-carbon archive changed")
    old = validate_record(path)
    m = wall_model()
    if m!=read_json(PREVIOUS/"model.json"):
        raise ValueError("Contact weighting or carrier exchange changed")
    baselines,points,_ = baseline_inputs()
    cases = {}
    for size in CASES:
        setup = m["sizes"][size]
        p = read_json(baselines[size]["synthetic"])["inputs"]
        cases[size] = dict(parameters=p,trace=points[size],mix=setup["mix"],
                          weights=dict(wall=setup["allocations"]["advection_rank"],
                                       flow=m["closure"]["flow_fractions"]),
                          wall_exposure_summary=setup["allocation_summary"]["advection_rank"],
                          separate_carbon_control=old["figure7"][size]["advection_rank"])
    return dict(closure=m["closure"],dry_carbon_mass_kg=m["dry_carbon_mass_kg"],sizes=cases,
                archived_contact_column_sha256=PREVIOUS_SHA)


def conservation(r,p,trace,weights):
    validate_solution(r,p,trace)
    area = 6*trace["dry_inventory_kg"]/(p["particle_density_kg_m3"]*p["particle_diameter_m"])
    np.testing.assert_allclose(r["external_particle_area_m2"],area,rtol=1e-12,atol=0)
    np.testing.assert_allclose(r["contact_area_by_path_m2"],area*np.array(weights),rtol=1e-12,atol=0)
    np.testing.assert_allclose(r["total_film_conductance_m3_s"],area*p["film_m_s"],rtol=1e-12,atol=0)
    np.testing.assert_allclose(r["cumulative_film_uptake_kg"]-r["cumulative_film_release_kg"],
                               r["adsorbed_inventory_kg"][-1],rtol=1e-8,atol=1e-16)


def full_refinement(coarse,fine,feed):
    result = refine(coarse,fine,feed)
    a,b = (crossings(r,feed) for r in (coarse,fine))
    # A crossing can move between the two windows near a flat plateau.
    # Keep it unresolved until full runs agree on the timing.
    for key in a:
        x,y = a[key],b[key]
        if x is not None and y is not None and min(x,y)<=SHORT_S<max(x,y):
            result["errors"][key+"_cross_window_relative"] = abs(x/y-1)
            result["passed"] &= result["errors"][key+"_cross_window_relative"]<=.02
    return result


def review(base,item):
    full = [read_json(base/name)["result"] for name in item["full"]]
    short = [read_json(base/name)["result"] for name in item["short"]]
    p = read_json(base/item["full"][-1])["inputs"]
    feed = p["inlet_concentration_kg_m3"]
    checks = [full_refinement(a,b,feed) for a,b in zip(full,full[1:])]
    startup = [checks[-1]["startup"]]
    path = [full[-1],*short]
    startup.extend(startup_refinement(a,b,feed) for a,b in zip(path,path[1:]))
    join = join_difference(full[-1],short[-1],feed) if short else None
    passed = checks[-1]["passed"] and startup[-1]["passed"] and (join is None or join<=.005)
    return dict(full_refinement=checks,startup_refinement=startup,join=join,numerically_qualified=passed)


def score(base,frozen,prepared,pixels_hash):
    if sha(PIXELS)!=pixels_hash:
        raise ValueError("Figure 7 markers changed")
    result = {}
    for size,(series,_) in CASES.items():
        report = {}
        for case,item in frozen["cases"][size].items():
            full = read_json(base/item["full"][-1])["result"]
            short = read_json(base/item["short"][-1])["result"] if item["short"] else None
            m = compare(scoring_curve(full,short),series)
            report[case] = dict(metrics=m,physical_screen_passed=column_passed(m),
                                numerical_screen_passed=item["review"]["numerically_qualified"],
                                total_mass_kg=full["physical_adsorbent_mass_kg"],
                                total_area_m2=full["external_particle_area_m2"],
                                film_flux_residual=full["normalized_interface_flux_mismatch"],
                                gross_film_release_kg=full["cumulative_film_release_kg"],
                                gross_film_uptake_kg=full["cumulative_film_uptake_kg"])
        report["separate_carbon"] = prepared["sizes"][size]["separate_carbon_control"]
        result[size] = report
    return result


def run():
    prepared = model()
    _,source,_,_ = source_inputs()
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")+"-"+uuid4().hex[:8]
    base = OUTPUT/"runs"/stamp
    base.mkdir(parents=True)
    files = [HERE/name for name in ("run.py","protocol.md")]
    files += [JULIA/name for name in ("SharedContactColumn.jl","test_shared_contact_column.jl","FilmParticle.jl","ExchangeColumn.jl")]
    files += [HERE.parent/name for name in ("contact_column/run.py","contact_column/inputs.py",
                                           "exchange_column/run.py","alhamed/study.py")]
    manifest = source_manifest(files)
    for name in manifest:
        target = base/"source_snapshot"/name
        target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(ROOT/name,target)
    write_json(base/"model.json",prepared)
    frozen = dict(model_sha256=sha(base/"model.json"),cases={})
    cl = prepared["closure"]
    print("Shared-contact predictions",base,flush=True)
    for size,setup in prepared["sizes"].items():
        frozen["cases"][size] = {}
        p0 = setup["parameters"]
        for case,weights in setup["weights"].items():
            item = dict(full=[],short=[])
            feed = p0["inlet_concentration_kg_m3"]
            def save(p,group,index):
                r = solve(cl,setup["mix"],weights,p)
                conservation(r,p,setup["trace"],weights)
                name = f"{size}_{case}_{group}{index}.json"
                write_json(base/name,dict(inputs=p,result=r))
                item[group].append(name)
                print(size,case,group,index,crossings(r,feed),"seconds",r["total_s"],flush=True)
                return r
            previous = None
            for index,(sub,nr,dt,early) in enumerate(LEVELS):
                p = {**p0,"subcells":sub,"radial_cells":nr,"dt_s":dt,"startup_dt_s":early}
                full = save(p,"full",index)
                if previous is not None:
                    changes = full_refinement(previous,full,feed)
                    print(size,case,changes,flush=True)
                    if changes["passed"]:
                        break
                previous = full
            if not changes["startup"]["passed"]:
                previous = full
                for index,(sub,nr,dt) in enumerate(STARTUP):
                    early = {**p,"subcells":sub,"radial_cells":nr,"startup_dt_s":dt,
                             "t_end_s":SHORT_S,"startup_s":SHORT_S}
                    short = save(early,"short",index)
                    if startup_refinement(previous,short,feed)["passed"]:
                        break
                    previous = short
            item["review"] = review(base,item)
            item["hashes"] = {name:sha(base/name) for name in item["full"]+item["short"]}
            frozen["cases"][size][case] = item
    write_json(base/"frozen_predictions.json",frozen)
    frozen_hash = sha(base/"frozen_predictions.json")
    print("Predictions frozen before Figure 7 scoring",frozen_hash,flush=True)
    scored = score(base,frozen,prepared,source["pixel_csv_sha256"])
    record = dict(schema="pore2process.shared-carbon-contact.v1",run_id=stamp,source_sha256=manifest,
        artifacts_sha256={str(p.relative_to(base)):sha(p) for p in base.rglob("*") if p.is_file()},
        frozen_prediction_sha256=frozen_hash,source_pixel_sha256=source["pixel_csv_sha256"],
        figure7=scored,versions=runtime_versions(),new_column_fitted_parameters=0,new_OpenLB_runs=0,new_DWSIM_runs=0)
    path = base/"summary.json"
    write_json(path,record)
    write_json(OUTPUT/"latest.json",dict(summary=str(path.relative_to(ROOT))))
    print("Complete",path,flush=True)


def check(path):
    record = validate_record(path)
    base = path.parent
    prepared = model()
    frozen = read_json(base/"frozen_predictions.json")
    if (read_json(base/"model.json")!=prepared or sha(base/"model.json")!=frozen["model_sha256"]
            or sha(base/"frozen_predictions.json")!=record["frozen_prediction_sha256"]):
        raise ValueError("Frozen model or prediction manifest changed")
    for size,groups in frozen["cases"].items():
        setup = prepared["sizes"][size]
        for case,item in groups.items():
            for name,digest in item["hashes"].items():
                if sha(base/name)!=digest:
                    raise ValueError("Frozen prediction changed")
                r = read_json(base/name)
                conservation(r["result"],r["inputs"],setup["trace"],setup["weights"][case])
            if review(base,item)!=item["review"]:
                raise ValueError("Numerical refinement differs on replay")
    if score(base,frozen,prepared,record["source_pixel_sha256"])!=record["figure7"]:
        raise ValueError("Figure 7 score differs on replay")
    print(dict(sources=len(record["source_sha256"]),artifacts=len(record["artifacts_sha256"])))


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check",type=Path)
    args = parser.parse_args()
    check(args.check.resolve()) if args.check else run()
