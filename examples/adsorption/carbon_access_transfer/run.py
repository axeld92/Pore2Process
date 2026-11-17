"""Freeze predictions at new operating points and score Alhamed Figures 4–5."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
import subprocess
import sys
from uuid import uuid4

import numpy as np

from examples.adsorption.alhamed.study import crossing
from examples.adsorption.carbon_access import run as access
from examples.adsorption.disordered_memory_column import run as column
from examples.adsorption.dispersion_holdout.run import execute
from examples.adsorption.exchange_column.run import source_inputs
from examples.adsorption.figure7_retention.run import checked_model, model_snapshot
from examples.adsorption.film_particle.bridge import call
from examples.adsorption.shared_contact.run import full_refinement
from examples.adsorption.transport_accuracy.run import validate_record
from pore2process_demo.bridge import runtime_versions
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT
from pore2process_demo.provenance import source_manifest

from . import kernel, source
from .probe import OUTPUT, REFERENCE_FLOW

HERE = Path(__file__).resolve().parent
ACCESS = access.OUTPUT/"runs/20261117T090234Z-234da0ca/summary.json"
ACCESS_SHA = "37480ecb42f74a66a02785ac2221f5099822a2fdfea8ca86fa4202fd1223ca5f"


def transferred_scenarios(reference_film,new_film):
    cases = {}
    for name,s in access.SCENARIOS.items():
        ratio = s["slow_ratio"]
        resistance = (1/ratio-1)/reference_film if ratio>0 else None
        new_ratio = 1/(1+resistance*new_film) if resistance is not None else 0.
        cases[name] = dict(inputs=dict(direct_fraction=s["direct_fraction"],slow_ratio=new_ratio),
                           added_resistance_s_m=resistance,reference_slow_ratio=ratio)
    return cases


def plan(probes):
    if sha(ACCESS)!=ACCESS_SHA:
        raise ValueError("The original unfitted access study changed")
    validate_record(ACCESS)
    probe_record = validate_record(probes)
    bed,_ = column.inputs()
    ref = bed["cases"]["d800"]
    _,_,sp,_ = source_inputs()
    cases = {}
    for name,condition in source.operating_cases().items():
        Q = condition["flow_mL_min"]*1e-6/60
        L = condition["bed_cm"]*.01
        A = bed["bed_area_m2"]
        cells = round(L/ref["column_section_length_m"])
        blocks,terminal = divmod(cells,8)
        film = call("film_coefficient",diameter_m=.0008,superficial_velocity_m_s=Q/A,
            viscosity_Pa_s=sp["viscosity_Pa_s"],density_kg_m3=sp["density_kg_m3"],
            molecular_diffusivity_m2_s=sp["molecular_diffusivity_m2_s"],
            porosity=sp["porosity"],interstitial=False,multiplier=1.)
        p = dict(ref["parameters"],blocks=blocks,terminal_sections=terminal,
            dry_mass_kg=sp["bulk_density_kg_m3"]*A*L,liquid_volume_m3=sp["porosity"]*A*L,
            flow_m3_s=Q,film_m_s=film["film_m_s"],t_end_s=36000.)
        cases[name] = dict(**condition,parameters=p,sections=cells,
            Pe=Q/A/sp["porosity"]*(L/cells)/sp["molecular_diffusivity_m2_s"],
            scenarios=transferred_scenarios(ref["parameters"]["film_m_s"],p["film_m_s"]))
    return dict(cases=cases,probe_summary=probes.relative_to(ROOT).as_posix(),probe_sha256=sha(probes),
                tracer_holdout=probe_record["cases"],parent_access_sha256=ACCESS_SHA)


def simulate_case(base,name):
    setup = read_json(base/"plan.json")
    probes = ROOT/setup["probe_summary"]
    condition = setup["cases"][name]
    flow = condition["flow_mL_min"]
    feed = condition["parameters"]["inlet_concentration_kg_m3"]
    cases = {}
    for scenario,spec in condition["scenarios"].items():
        s = spec["inputs"]
        files,checks = [],[]
        previous = None
        for i,(radial,dt,early) in enumerate(column.LEVELS):
            p = dict(condition["parameters"],radial_cells=radial,dt_s=dt,startup_dt_s=early)
            r = kernel.solve(flow,probes.parent,p,s)
            verify_run(r,p,s,condition)
            path = f"{name}_{scenario}_{i}.json"
            write_json(base/path,dict(inputs=p,scenario=s,result=r))
            files.append(path)
            print(path,"t50 min",(crossing(r["time_s"],np.array(r["outlet_concentration_kg_m3"])/feed,.5) or 0)/60,
                  "wall seconds",r["wall_s"],flush=True)
            if previous is not None:
                check = full_refinement(previous,r,feed)
                checks.append(check)
                if check["passed"] and check["startup"]["passed"]:
                    break
            previous = r
        cases[scenario] = dict(files=files,selected=files[-1],refinement=checks,
            numerically_qualified=checks[-1]["passed"] and checks[-1]["startup"]["passed"])
    write_json(base/f"{name}_cases.json",cases)


def verify_run(r,p,s,condition):
    if r["sections"]!=condition["sections"]:
        raise ValueError("Column has the wrong number of tracer sections")
    access.verify(r,{key:value for key,value in p.items() if key!="terminal_sections"},s)


def crossing_metrics(pred,obs):
    scores = {}
    for key,limit in (("t05_s",.25),("t50_s",.20)):
        bounds = obs["reading_"+key]
        relative = abs(pred[key]/obs[key]-1) if pred[key] is not None and obs[key] else None
        passed = None if bounds is None else relative is not None and relative<=limit
        scores[key] = dict(predicted_s=pred[key],relative_error=relative,passed=passed)
    return scores


def compare(base,setup,frozen):
    markers = source.markers()
    metadata = read_json(source.EXTRACTION)
    dose = {}
    for key,m in markers.items():
        p = setup["cases"][m["case"]]["parameters"]
        dose[key] = m["time"][-1]*p["flow_m3_s"]*p["inlet_concentration_kg_m3"]/p["dry_mass_kg"]
    common = {i:float(min(dose[key] for key,m in markers.items() if m["figure"]==i)) for i in source.SERIES}
    comparisons = {}
    for key,m in markers.items():
        condition = setup["cases"][m["case"]]
        p = condition["parameters"]
        feed_rate = p["flow_m3_s"]*p["inlet_concentration_kg_m3"]
        common_time = common[m["figure"]]*p["dry_mass_kg"]/feed_rate
        common_time = min(common_time,m["time"][-1])
        obs = source.observation(m,feed_rate,p["dry_mass_kg"],common_time)
        scenarios = {}
        for name,item in frozen["cases"][m["case"]].items():
            saved = read_json(base/item["selected"])
            r = saved["result"]
            model = checked_model(saved,dict(p,feed_rate_kg_s=feed_rate))
            t,c = model["time_s"],model["outlet_fraction"]
            predicted = {f"t{i:02d}_s":crossing(t,c,i/100) for i in (5,50)}
            error = np.interp(m["time"],t,c)-m["fraction"]
            timing = crossing_metrics(predicted,obs)
            rmse = float(np.sqrt(np.mean(error**2)))
            state = model_snapshot(model,obs["end_s"],feed_rate)
            retained = state["retained_kg"]/p["dry_mass_kg"]
            scenarios[name] = dict(**timing,rmse=rmse,retained_kg_kg=retained,
                retained_in_reading_range=obs["reading_retention_kg_kg"][0]<=retained<=obs["reading_retention_kg_kg"][1],
                common_dose_retained_kg_kg=float(np.interp(common_time,t,model["retained_kg"])/p["dry_mass_kg"]),
                normalized_mass_residual=r["normalized_mass_residual"],
                passed=rmse<=.1 and timing["t50_s"]["passed"] is True and timing["t05_s"]["passed"] is not False)
        comparisons[key] = dict(operating_case=m["case"],figure=m["figure"],last_dose_kg_kg=float(dose[key]),
                                observed=obs,scenarios=scenarios)
    a,b = (comparisons[f"f{figure}_q233_l25"]["observed"] for figure in (4,5))
    duplicate = dict(midpoint_difference_s=abs(a["t50_s"]-b["t50_s"]),
                     last_marker_loading_difference_kg_kg=abs(a["retained_kg_kg"]-b["retained_kg_kg"]))
    return dict(source=metadata,common_dose_kg_kg={str(k):v for k,v in common.items()},
                shared_condition_digitization=duplicate,comparisons=comparisons)


def run(probes):
    setup = plan(probes)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")+"-"+uuid4().hex[:8]
    base = OUTPUT/"runs"/stamp
    base.mkdir(parents=True)
    copies = [HERE/name for name in ("run.py","kernel.py","source.py","probe.py","protocol.md")]
    sources = source_manifest(copies+[access.HERE/"cohorts.jl",access.HERE/"kernel.py",kernel.access.ORIGINAL])
    write_json(base/"plan.json",setup)
    (base/"solver.jl").write_text(kernel.source())
    print("Frozen operating and access inputs",base,sha(base/"plan.json"),flush=True)
    def process(name):
        cmd = [sys.executable,"-m","examples.adsorption.carbon_access_transfer.run","--case",name,"--base",str(base)]
        elapsed = execute(cmd,base/f"{name}.log",timeout=2400)
        print(name,"prediction complete",elapsed,flush=True)
    with ThreadPoolExecutor(max_workers=3) as pool:
        list(pool.map(process,setup["cases"]))
    frozen = dict(plan_sha256=sha(base/"plan.json"),
                  cases={name:read_json(base/f"{name}_cases.json") for name in setup["cases"]})
    frozen["curves_sha256"] = {file:sha(base/file) for scenarios in frozen["cases"].values()
                              for item in scenarios.values() for file in item["files"]}
    path = base/"frozen.json"
    write_json(path,frozen)
    print("Predictions frozen before reading marker fractions",sha(path),flush=True)
    scored = compare(base,setup,frozen)
    summary = base/"summary.json"
    write_json(summary,dict(schema="pore2process.carbon-access-transfer.v1",
        run_id=stamp,frozen_sha256=sha(path),source_sha256=sources,figure_parameters_fitted=0,
        artifacts_sha256={p.name:sha(p) for p in sorted(base.iterdir())},
        versions=runtime_versions(),score=scored))
    write_json(OUTPUT/"latest.json",dict(summary=summary.relative_to(ROOT).as_posix()))
    print(summary,flush=True)


def check(path):
    record = validate_record(path)
    base = path.parent
    setup,frozen = [read_json(base/name) for name in ("plan.json","frozen.json")]
    probe = ROOT/setup["probe_summary"]
    validate_record(probe)
    if (sha(probe)!=setup["probe_sha256"] or sha(base/"frozen.json")!=record["frozen_sha256"]
            or sha(base/"plan.json")!=frozen["plan_sha256"] or (base/"solver.jl").read_text()!=kernel.source()
            or record["score"]!=compare(base,setup,frozen)):
        raise ValueError("Frozen holdout scores differ on replay")
    for name,scenarios in frozen["cases"].items():
        for scenario,item in scenarios.items():
            condition = setup["cases"][name]
            runs = [read_json(base/file) for file in item["files"]]
            for r in runs:
                if r["scenario"]!=condition["scenarios"][scenario]["inputs"]:
                    raise ValueError("Saved scenario parameters differ")
                verify_run(r["result"],r["inputs"],r["scenario"],condition)
            feed = condition["parameters"]["inlet_concentration_kg_m3"]
            if [full_refinement(a["result"],b["result"],feed) for a,b in zip(runs,runs[1:])]!=item["refinement"]:
                raise ValueError("Refinement differs on replay")
    print(dict(sources=len(record["source_sha256"]),artifacts=len(record["artifacts_sha256"])))


if __name__=="__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--probes",type=Path)
    p.add_argument("--check",type=Path)
    p.add_argument("--case")
    p.add_argument("--base",type=Path)
    args = p.parse_args()
    if args.case:
        simulate_case(args.base.resolve(),args.case)
    elif args.check:
        check(args.check.resolve())
    else:
        run(args.probes.resolve() if args.probes else ROOT/read_json(OUTPUT/"latest_probes.json")["summary"])
