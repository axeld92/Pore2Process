"""Freeze age-conditioned N160 tracer predictions and check the archived holdout."""
import argparse
from datetime import datetime, timezone
from pathlib import Path
from tempfile import TemporaryDirectory
from uuid import uuid4

import numpy as np

from examples.adsorption.disordered_exchange.run import PARENT_SHA, archive_case, parent_inputs
from examples.adsorption.disordered_section_exchange.run import (
    LIMITS, errors, train_columns, verify_saved_field,
)
from examples.adsorption.dispersion_holdout.run import execute
from examples.adsorption.transport_accuracy.run import validate_record
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT
from pore2process_demo.provenance import source_manifest

HERE = Path(__file__).resolve().parent
OUTPUT = ROOT/"results/adsorption/disordered_travel_memory"
PROBES = ROOT/"results/adsorption/disordered_section_exchange/runs/20261113T173037Z-5d62c7c1"
PROBES_SHA = "cad8ffcd4934850da3cb95744d0ab9ca3c6251c7871ee7e2ba5c7132552100ce"
MODELS = ("baseline","two","four","eight","four_first","four_second")
PRODUCTS = ("rules.f64",*(f"{name}.csv" for name in MODELS))


def inputs():
    parent,plan = parent_inputs()
    verify_saved_field(parent,plan)
    path = PROBES/"summary.json"
    if sha(path)!=PROBES_SHA:
        raise ValueError("Section-specific run summary changed")
    probe_run = validate_record(path)
    probe_plan = read_json(PROBES/"probe_plan.json")
    count = probe_plan["count"]
    n = plan["marker_count"]
    if count!=2*n or probe_plan["seed"] < plan["seed"]+n:
        raise ValueError("Probe ensemble must contain two halves independent of the archived holdout")
    if probe_run["parent_sha256"]!=PARENT_SHA or probe_run["geometry_sha256"]!=parent["analysis"]["geometry_sha256"]:
        raise ValueError("Probe run refers to different N160 geometry or trajectories")
    return parent,plan,count


def probe_columns(name,count):
    data = np.fromfile(PROBES/(name+".f64"),dtype="<f8").reshape(count,16)
    return data[:,:8].copy(order="C").tobytes()


def solve(case, spec, training, probes):
    command = [*spec["command"][:4],str(HERE/"solve.jl"),str(training),str(probes),str(case)]
    timing = execute(command,case/"solve.log",timeout=120)
    write_json(case/"execution.json",dict(command=command,**timing))


def predictions(base, plan):
    cases = {}
    for name,spec in plan["cases"].items():
        case = base/name
        model = {key:np.loadtxt(case/(key+".csv"),delimiter=",") for key in MODELS}
        previous = np.loadtxt(PROBES/name/"sectional.csv",delimiter=",")
        np.testing.assert_allclose(model["baseline"],previous,rtol=1e-12,atol=1e-15)
        cases[name] = dict(Pe=spec["Pe"],dt_s=spec["dt_s"],
            training_sha256=sha(case/"training.f64"),probe_passages_sha256=sha(case/"probe_passages.f64"),
            rules_sha256=sha(case/"rules.f64"),**{key:value.tolist() for key,value in model.items()})
    return cases


def score(base,parent,plan,predicted):
    cases = {}
    for name,model in predicted.items():
        data = archive_case(parent,plan,name)
        if (base/name/"training.f64").read_bytes()!=data[:,:2].tobytes():
            raise ValueError("Training input differs from archived sections 1–2")
        arrivals = data[:,2:8]
        if not np.isfinite(arrivals).all() or np.any(np.diff(data[:,:8],axis=1)<=0):
            raise ValueError("Incomplete or unordered holdout passages")
        observed = np.column_stack((arrivals.mean(axis=0),arrivals.var(axis=0)))
        if np.any(observed<=0):
            raise ValueError("Expected nonzero holdout arrival means and variances")
        cases[name] = dict(observed=observed.tolist(),**{key:errors(observed,model[key]) for key in MODELS})
    return dict(limits=LIMITS,cases=cases,passed=all(v["four"]["variance_passed"] for v in cases.values()))


def frozen(base,plan,sources):
    return dict(parent_sha256=PARENT_SHA,probe_summary_sha256=PROBES_SHA,
                source_sha256=sources,cases=predictions(base,plan))


def run():
    parent,plan,count = inputs()
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")+"-"+uuid4().hex[:8]
    base = OUTPUT/"runs"/stamp
    base.mkdir(parents=True)
    sources = source_manifest([HERE/name for name in ("run.py","model.jl","solve.jl","test_model.jl","protocol.md")]+[
        HERE.parent/name for name in ("disordered_section_exchange/model.jl","disordered_section_exchange/run.py",
                                      "disordered_exchange/model.jl","disordered_exchange/run.py",
                                      "dispersion_holdout/run.py","transport_accuracy/run.py")])
    for name,spec in plan["cases"].items():
        case = base/name
        case.mkdir()
        training,probes = case/"training.f64",case/"probe_passages.f64"
        training.write_bytes(train_columns(name,plan["marker_count"]))
        probes.write_bytes(probe_columns(name,count))
        solve(case,spec,training,probes)
    prediction = frozen(base,plan,sources)
    path = base/"frozen_prediction.json"
    write_json(path,prediction)
    digest = sha(path)
    print("Predictions frozen before scoring sections 3–8",digest,flush=True)
    analysis = score(base,parent,plan,prediction["cases"])
    record = base/"summary.json"
    write_json(record,dict(schema="pore2process.disordered_travel_memory.v1",
        run_id=stamp,parent_sha256=PARENT_SHA,probe_summary_sha256=PROBES_SHA,source_sha256=sources,
        frozen_prediction_sha256=digest,analysis=analysis,new_native_flow_runs=0,
        new_tracer_probes=0,new_adsorption_runs=0,figure7_fits=0,
        artifacts_sha256={p.relative_to(base).as_posix():sha(p) for p in sorted(base.rglob("*")) if p.is_file()}))
    write_json(OUTPUT/"latest.json",dict(summary=str(record.relative_to(ROOT))))
    for name,case in analysis["cases"].items():
        print(name,{key:round(100*case[key]["maximum_variance_error"],4) for key in MODELS})
    print(record)


def check(path):
    record = validate_record(path)
    parent,plan,count = inputs()
    base = path.parent
    prediction = frozen(base,plan,record["source_sha256"])
    saved = base/"frozen_prediction.json"
    if (record["parent_sha256"]!=PARENT_SHA or record["probe_summary_sha256"]!=PROBES_SHA
            or record["frozen_prediction_sha256"]!=sha(saved) or read_json(saved)!=prediction
            or record["analysis"]!=score(base,parent,plan,prediction["cases"])):
        raise ValueError("Frozen prediction or score differs on replay")
    for name,spec in plan["cases"].items():
        case = base/name
        probes = case/"probe_passages.f64"
        if probes.read_bytes()!=probe_columns(name,count):
            raise ValueError("Saved age ranks must come from the independent probes")
        with TemporaryDirectory(prefix="pore2process-age-") as tmp:
            replay = Path(tmp)
            solve(replay,spec,case/"training.f64",probes)
            for product in PRODUCTS:
                if (replay/product).read_bytes()!=(case/product).read_bytes():
                    raise ValueError(f"Julia replay differs for {name}/{product}")
    print(dict(sources=len(record["source_sha256"]),artifacts=len(record["artifacts_sha256"]),
               passed=record["analysis"]["passed"]))


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check",type=Path)
    args = parser.parse_args()
    check(args.check.resolve()) if args.check else run()
