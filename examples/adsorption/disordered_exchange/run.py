"""Freeze and score the Julia two-section disordered tracer reduction."""
import argparse
from datetime import datetime, timezone
from pathlib import Path
import shutil
from uuid import uuid4

import numpy as np

from examples.adsorption.dispersion_holdout.run import execute
from examples.adsorption.disordered_packing.compare import profile_check
from examples.adsorption.disordered_n160.native import OUTPUT as N160
from examples.adsorption.path_exchange.predict import ranks
from examples.adsorption.repeated_geometry.predict import observe
from examples.adsorption.transport_accuracy.run import validate_record
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT
from pore2process_demo.provenance import source_manifest

HERE = Path(__file__).resolve().parent
OUTPUT = ROOT/"results/adsorption/disordered_exchange"
PARENT = N160/"runs/20261113T130058Z-5ba9e799"
PARENT_SHA = "c352eae8ef39f34c9bf5fdccb3185a029a522d584decf2b1ab03e305b6e59755"
FOLDER = "disordered_N160"
METRICS = ("survival", "fast_survival", "conditional_capture", "fast_conditional_capture",
           "mean_arrival_s", "variance_s2", "fast_persistence")
MODEL_FILES = ("counts","clear","delay","delay2","transition","wallfree","factorized",
               "flux","states","histories","fast_share","pair","receiver")


def parent_inputs():
    if sha(PARENT/"summary.json")!=PARENT_SHA:
        raise ValueError("N160 parent summary changed")
    summary = read_json(PARENT/"summary.json")
    path = PARENT/FOLDER/"tracer_plan.json"
    if sha(path)!=summary["analysis"]["tracer_plan_sha256"]:
        raise ValueError("N160 plan changed")
    plan = read_json(path)
    return summary,plan


def archive_case(summary,plan,name):
    relative = f"{FOLDER}/{name}.f64"
    path = PARENT/relative
    if sha(path)!=summary["artifacts_sha256"][relative]:
        raise ValueError("Pinned N160 marker data changed")
    return np.memmap(path,mode="r",dtype="<f8",shape=(plan["marker_count"],16))


def training(data):
    return np.asarray(data[:,[0,1,8,9]],dtype="<f8",order="C")


def run_models(base,summary,plan):
    for name,case in plan["cases"].items():
        directory = base/name
        directory.mkdir()
        source = directory/"training.f64"
        training(archive_case(summary,plan,name)).tofile(source)
        julia,_,*flags = case["command"][:4]
        cmd = [julia,*flags,str(HERE/"solve.jl"),str(source),str(directory)]
        timing = execute(cmd,directory/"solve.log",timeout=120)
        write_json(directory/"execution.json",dict(command=cmd,**timing))


def freeze(base,summary,plan):
    cases = {}
    for name,spec in plan["cases"].items():
        directory = base/name
        frozen = training(archive_case(summary,plan,name)).tobytes()
        if (directory/"training.f64").read_bytes()!=frozen:
            raise ValueError("Training data differ from N160 sections 1–2")
        models = {key:np.loadtxt(directory/(key+".csv"),delimiter=",",ndmin=1).tolist()
                  for key in MODEL_FILES}
        cases[name] = dict(Pe=spec["Pe"],dt_s=spec["dt_s"],training_sha256=sha(directory/"training.f64"),
                           **models)
    return dict(parent_sha256=PARENT_SHA,dry_carbon_mass_kg=summary["analysis"]["dry_carbon_mass_kg"],
                marker_count=plan["marker_count"],cases=cases)


def observations(data,case):
    actual = observe(data)
    times = data[:,:8]
    if actual["recovery"][-1]!=1.:
        raise ValueError("Inert moments require complete trajectories")
    dt = np.diff(times,axis=1,prepend=0.)
    limits = np.rint(len(data)*np.arange(1,11)/10).astype(int)
    labels = np.column_stack([np.searchsorted(limits,ranks(d),side="right") for d in dt.T])
    persistent = np.logical_and.accumulate(labels<9,axis=1)[labels[:,0]<9]
    training_joint = np.array(case["counts"])/len(data)
    variation = []
    for donor,receiver in zip(labels.T[:-1],labels.T[1:],strict=True):
        table = np.bincount(donor*10+receiver,minlength=100).reshape(10,10)/len(data)
        variation.append(float(.5*np.abs(table-training_joint).sum()))
    return dict(**actual,mean_arrival_s=times.mean(axis=0).tolist(),
                variance_s2=times.var(axis=0).tolist(),fast_persistence=persistent.mean(axis=0).tolist(),
                seam_exchange_total_variation=variation)


def compare(actual,predicted):
    prior = {key:values[:2] for key,values in actual.items() if key in METRICS}
    profiles = {key:prior[key]+list(values) for key,values in zip(METRICS,np.array(predicted).T,strict=True)}
    wall = profile_check(actual,profiles,.02)
    limits = dict(mean_arrival_s=.05,variance_s2=.10)
    timing_errors = {key:float(max(abs(np.asarray(profiles[key][2:])/actual[key][2:]-1))) for key in limits}
    timing = dict(errors=timing_errors,limits=limits,passed=all(timing_errors[k]<=limits[k] for k in limits))
    return dict(wall=wall,timing=timing,profiles=profiles,passed=wall["passed"] and timing["passed"])


def evaluate(summary,plan,frozen):
    cases = {}
    for name,case in frozen["cases"].items():
        actual = observations(archive_case(summary,plan,name),case)
        observed_training = np.array([actual["survival"][:2],actual["fast_survival"][:2]]).T
        np.testing.assert_allclose(case["histories"],observed_training,rtol=0.,atol=1e-15)
        cases[name] = dict(observed=actual,**{key:compare(actual,case[key]) for key in ("pair","receiver")})
    return dict(cases=cases,passed=all(case["pair"]["passed"] for case in cases.values()))


def run():
    summary,plan = parent_inputs()
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")+"-"+uuid4().hex[:8]
    base = OUTPUT/"runs"/stamp
    base.mkdir(parents=True)
    sources = source_manifest([HERE/name for name in ("run.py","model.jl","solve.jl","protocol.md")]+[
        HERE.parent/name for name in ("path_exchange/predict.py","repeated_geometry/predict.py",
                                      "disordered_packing/compare.py")])
    for name in sources:
        target = base/"source_snapshot"/name
        target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(ROOT/name,target)
    run_models(base,summary,plan)
    frozen = freeze(base,summary,plan)
    path = base/"frozen_prediction.json"
    write_json(path,frozen)
    digest = sha(path)
    print("Predictions frozen before evaluating sections 3–8:",digest,flush=True)
    result = evaluate(summary,plan,frozen)
    record = base/"summary.json"
    write_json(record,dict(schema="pore2process.disordered-exchange.v1",run_id=stamp,
        frozen_prediction_sha256=digest,source_sha256=sources,parent_sha256=PARENT_SHA,
        artifacts_sha256={str(p.relative_to(base)):sha(p) for p in base.rglob("*") if p.is_file()},
        analysis=result,new_adsorption_runs=0,figure7_fits=0))
    write_json(OUTPUT/"latest.json",dict(summary=str(record.relative_to(ROOT))))
    for name,case in result["cases"].items():
        print(name,case["pair"]["wall"],case["pair"]["timing"])
    print(record)


def check(path):
    record = validate_record(path)
    summary,plan = parent_inputs()
    frozen_path = path.parent/"frozen_prediction.json"
    frozen = read_json(frozen_path)
    if (record["frozen_prediction_sha256"]!=sha(frozen_path)
            or frozen!=freeze(path.parent,summary,plan)
            or record["analysis"]!=evaluate(summary,plan,frozen)):
        raise ValueError("The frozen disordered tracer prediction differs on replay")
    print(dict(sources=len(record["source_sha256"]),artifacts=len(record["artifacts_sha256"]),
               passed=record["analysis"]["passed"]))


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check",type=Path)
    args = parser.parse_args()
    check(args.check.resolve()) if args.check else run()
