"""Compare isolated and slowly accessible carbon at fixed particle physics."""
import argparse
from datetime import datetime, timezone
from pathlib import Path
from uuid import uuid4

import numpy as np

from examples.adsorption.alhamed.study import compare
from examples.adsorption.contact_column.run import crossings
from examples.adsorption.disordered_memory_column import run as column
from examples.adsorption.disordered_repeat_join.run import archive, archived_parity, PARENT_SHA
from examples.adsorption.exchange_column.run import CASES, PIXELS
from examples.adsorption.figure7_retention.run import REFERENCE as RETENTION, checked_model, model_snapshot
from examples.adsorption.film_particle.metrics import column_passed
from examples.adsorption.shared_contact.run import full_refinement
from examples.adsorption.transport_accuracy.run import validate_record
from pore2process_demo.bridge import runtime_versions
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT
from pore2process_demo.provenance import source_manifest

from . import kernel

HERE = Path(__file__).resolve().parent
OUTPUT = ROOT/"results/adsorption/carbon_access"
SCENARIOS = dict(full=dict(direct_fraction=1.,slow_ratio=0.))
SCENARIOS.update({f"f{round(100*f)}_{name}":dict(direct_fraction=f,slow_ratio=ratio)
                  for f in (.5,.75) for name,ratio in (("isolated",0.),("slow01",.01),("slow10",.1))})


def verify(r,p,s):
    column.verify(r,p)
    f,ratio = s["direct_fraction"],s["slow_ratio"]
    np.testing.assert_allclose(r["carbon_by_cohort_kg"],np.array([f,1-f])*p["dry_mass_kg"],rtol=1e-12,atol=0.)
    np.testing.assert_allclose(r["equivalent_film_by_cohort_m_s"],p["film_m_s"]*np.array([1.,ratio]),
                               rtol=1e-12,atol=0.)
    np.testing.assert_allclose(r["effective_conductance_m3_s"],
        r["total_film_conductance_m3_s"]*(f+(1-f)*ratio),rtol=1e-12,atol=0.)
    groups = np.array(r["adsorbed_inventory_by_cohort_kg"])
    np.testing.assert_allclose(groups.sum(axis=0),r["adsorbed_inventory_kg"],rtol=1e-12,atol=1e-16)
    bc = p["affinity_m3_kg"]*p["inlet_concentration_kg_m3"]
    ceiling = p["dry_mass_kg"]*p["capacity_kg_kg"]*bc/(1+bc)
    np.testing.assert_allclose(r["equilibrium_adsorbed_ceiling_kg"],ceiling*(f if ratio==0 else 1.),
                               rtol=1e-12,atol=0.)
    if np.any(groups<0) or np.any(groups>np.array([f,1-f])[:,None]*ceiling+1e-14):
        raise ValueError("Particle cohort exceeded its physical capacity")
    if (ratio==0 or f==1.) and np.any(groups[1]!=0.):
        raise ValueError("A carbon cohort with no access adsorbed phenol")
    if ratio>0 and f<1. and groups[1,-1]<=0.:
        raise ValueError("A communicating slow cohort did not adsorb phenol")


def simulate(base,bed,baseline):
    cases = {}
    for size,settings in bed["cases"].items():
        p0 = settings["parameters"]
        cases[size] = {}
        for name,s in SCENARIOS.items():
            paths = []
            checks = []
            prior = None
            for i,(radial,dt,early) in enumerate(column.LEVELS):
                p = dict(p0,radial_cells=radial,dt_s=dt,startup_dt_s=early)
                r = kernel.solve(size,p,s)
                verify(r,p,s)
                if name=="full" and i<len(baseline["cases"][size]["refinement_files"]):
                    old = baseline["cases"][size]["refinement_files"][i]
                    archived_parity(size,old,r)
                path = base/f"{size}_{name}_{i}.json"
                write_json(path,dict(inputs=p,scenario=s,result=r))
                paths.append(path.name)
                feed = p["inlet_concentration_kg_m3"]
                print(size,name,i,crossings(r,feed),round(r["wall_s"],2),flush=True)
                if prior is not None:
                    report = full_refinement(prior,r,feed)
                    checks.append(report)
                    if report["passed"] and report["startup"]["passed"]:
                        break
                prior = r
            cases[size][name] = dict(files=paths,selected=paths[-1],refinement=checks,
                numerically_qualified=checks[-1]["passed"] and checks[-1]["startup"]["passed"])
    return cases


def evaluate(base,frozen):
    observed = read_json(RETENTION)
    if sha(PIXELS)!=observed["source"]["pixel_sha256"]:
        raise ValueError("Figure 7 marker digitization changed")
    basis = observed["basis"]
    scores = {}
    for size,(series,_) in CASES.items():
        source = observed["cases"][size]
        end = source["last_marker"]
        figure_mass = end["figure7_retained_kg"]
        readings = source["reading_scenarios_retention_kg"]
        cases = {}
        for name,item in frozen["cases"][size].items():
            run = read_json(base/item["selected"])
            r,s,p = run["result"],run["scenario"],run["inputs"]
            model = checked_model(run,basis)
            t = model["time_s"]
            state = model_snapshot(model,end["time_s"],basis["feed_rate_kg_s"])
            slow = float(np.interp(end["time_s"],t,r["adsorbed_inventory_by_cohort_kg"][1]))
            m = compare(r,series)
            cases[name] = dict(scenario=s,midpoint_s=m["t50"]["predicted_s"],
                midpoint_passed=m["t50_screen_passed"],early_crossing_s=m["t05"]["predicted_s"],
                rmse=m["RMSE"],full_figure7_screen_passed=column_passed(m),
                numerical_screen_passed=item["numerically_qualified"],snapshot=state,
                solid_ceiling_kg=r["equilibrium_adsorbed_ceiling_kg"],
                slow_adsorbed_kg=slow,excess_retention_kg=state["retained_kg"]-figure_mass,
                within_marker_reading_range=readings[0]<=state["retained_kg"]<=readings[1],
                total_dry_mass_kg=r["physical_adsorbent_mass_kg"],mass_residual=r["normalized_mass_residual"],
                film_flux_residual=r["normalized_interface_flux_mismatch"],
                clean_feed_minimum_slow_loading_time_s=(None if s["slow_ratio"]==0 else
                    basis["equilibrium_loading_at_feed_kg_kg"]*p["particle_density_kg_m3"]*p["particle_diameter_m"]/
                    (6*p["inlet_concentration_kg_m3"]*r["equivalent_film_by_cohort_m_s"][1])))
        scores[size] = dict(figure7_midpoint_s=source["experimental_midpoint"]["time_s"],
            figure7_retained_kg=figure_mass,nominal_window_s=end["time_s"],reading_scenarios_kg=readings,cases=cases)
    return scores


def run():
    baseline = archive()
    bed,_ = column.inputs()
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")+"-"+uuid4().hex[:8]
    base = OUTPUT/"runs"/stamp
    base.mkdir(parents=True)
    sources = source_manifest([HERE/name for name in ("run.py","kernel.py","cohorts.jl","protocol.md")]+[
        kernel.ORIGINAL,column.HERE/"run.py",HERE.parent/"figure7_retention/run.py"])
    plan = base/"scenario_plan.json"
    write_json(plan,dict(scenarios=SCENARIOS,bed=bed,parent_sha256=PARENT_SHA))
    generated = base/"CarbonAccessColumn.jl"
    generated.write_text(kernel.source())
    print(base,"scenario plan",sha(plan),flush=True)
    cases = simulate(base,bed,baseline)
    frozen = dict(plan_sha256=sha(plan),cases=cases,
                  curve_sha256={p.name:sha(p) for p in sorted(base.glob("d*.json"))})
    path = base/"frozen_predictions.json"
    write_json(path,frozen)
    print("Predictions frozen before reading Figure 7",sha(path),flush=True)
    scores = evaluate(base,read_json(path))
    summary = base/"summary.json"
    write_json(summary,dict(schema="pore2process.carbon_access.v1",run_id=stamp,
        source_sha256=sources,artifacts_sha256={p.name:sha(p) for p in sorted(base.iterdir())},
        parent_sha256=PARENT_SHA,retention_audit_sha256=sha(RETENTION),frozen_predictions_sha256=sha(path),
        figure7=scores,figure7_fitted_parameters=0,versions=runtime_versions()))
    write_json(OUTPUT/"latest.json",dict(summary=summary.relative_to(ROOT).as_posix()))
    print(summary,flush=True)


def check(summary):
    archive()
    record = validate_record(summary)
    base = summary.parent
    path = base/"frozen_predictions.json"
    frozen = read_json(path)
    if (record["parent_sha256"]!=PARENT_SHA or sha(RETENTION)!=record["retention_audit_sha256"]
            or sha(path)!=record["frozen_predictions_sha256"]
            or sha(base/"scenario_plan.json")!=frozen["plan_sha256"]
            or read_json(base/"scenario_plan.json")["scenarios"]!=SCENARIOS
            or (base/"CarbonAccessColumn.jl").read_text()!=kernel.source()
            or any(sha(base/name)!=digest for name,digest in frozen["curve_sha256"].items())
            or record["figure7"]!=evaluate(base,frozen)):
        raise ValueError("Frozen scenario comparison changed")
    for size,cases in frozen["cases"].items():
        for name,item in cases.items():
            runs = [read_json(base/file) for file in item["files"]]
            for r in runs:
                if r["scenario"]!=SCENARIOS[name]:
                    raise ValueError(f"Scenario inputs changed: {size}/{name}")
                verify(r["result"],r["inputs"],r["scenario"])
            feed = runs[0]["inputs"]["inlet_concentration_kg_m3"]
            if [full_refinement(a["result"],b["result"],feed)
                    for a,b in zip(runs,runs[1:])]!=item["refinement"]:
                raise ValueError(f"Refinement changed: {size}/{name}")
    print(dict(sources=len(record["source_sha256"]),artifacts=len(record["artifacts_sha256"])))


if __name__=="__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--check",type=Path)
    args = p.parse_args()
    check(args.check.resolve()) if args.check else run()
