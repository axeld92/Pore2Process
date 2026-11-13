"""Compose local Julia transport rules and score the frozen N160 prediction."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from uuid import uuid4

import numpy as np

from examples.adsorption.disordered_exchange.run import FOLDER, PARENT, PARENT_SHA, archive_case, parent_inputs
from examples.adsorption.dispersion_holdout.run import execute
from examples.adsorption.transport_accuracy.run import validate_record
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT
from pore2process_demo.provenance import source_manifest, verified_source_path

HERE = Path(__file__).resolve().parent
OUTPUT = ROOT/"results/adsorption/disordered_section_exchange"
PREDICTORS = ("stationary", "sectional", "probe_first", "probe_second")
RULE_TABLES = tuple(f"section_{section}_{name}.csv" for section in range(3,9)
                    for name in ("counts","transition","delay","delay2"))
LIMITS = {"mean": .05, "variance": .10}
TRACKER = ("disordered_packing/track.jl", "repeated_geometry/track.jl",
           "carbon_contact/track_contacts.jl", "path_exchange/track.jl")


def probe_plan(base, plan):
    count = 2*plan["marker_count"]
    seed = plan["seed"] + 4*plan["marker_count"]
    commands = {}
    for name, case in plan["cases"].items():
        command = case["command"]
        dx, diffusion, dt, stop, n, old_seed, grid, ends = command[7:]
        if (int(n),int(old_seed),int(grid)) != (plan["marker_count"],plan["seed"],160):
            raise ValueError("Unexpected parent tracer configuration")
        commands[name] = [*command[:4],str(HERE.parent/TRACKER[0]),str(PARENT/FOLDER/"compact"),
                          str(base/(name+".f64")),dx,diffusion,dt,stop,str(count),str(seed),grid,ends]
    return dict(count=count,seed=seed,commands=commands)


def verify_saved_field(summary, plan):
    for name in TRACKER:
        file = HERE.parent/name
        relative = file.relative_to(ROOT).as_posix()
        if verified_source_path(relative,summary["source_sha256"][relative]) != file:
            raise ValueError(f"Replay needs the archived tracker: {relative}")
    for name,digest in plan["compact_sha256"].items():
        if sha(PARENT/FOLDER/"compact"/name) != digest:
            raise ValueError(f"Saved compact field changed: {name}")
    for relative in ("plan.json",f"{FOLDER}/raw.json",f"{FOLDER}/compact/summary.json"):
        if sha(PARENT/relative)!=summary["artifacts_sha256"][relative]:
            raise ValueError(f"Saved geometry or flow metadata changed: {relative}")


def run_probes(base, probes):
    write_json(base/"probe_plan.json",probes)

    def trace(job):
        name, command = job
        result = execute(command,base/(name+".log"),timeout=1200)
        write_json(base/(name+"_execution.json"),result)
        print(name,"fresh probes complete",result,flush=True)

    with ThreadPoolExecutor(max_workers=2) as pool:
        list(pool.map(trace,probes["commands"].items()))


def train_columns(name, count):
    path = PARENT/FOLDER/(name+".f64")
    # Read only the training columns before the prediction freeze. The full
    # archive digest and equality of this extract are checked during scoring.
    with path.open("rb",buffering=0) as file:
        if path.stat().st_size != count*16*8:
            raise ValueError("Wrong archived marker count or format")
        training = bytearray()
        for row in range(count):
            file.seek(row*16*8)
            training.extend(file.read(2*8))
    return bytes(training)


def probe_durations(base, name, count):
    data = np.fromfile(base/(name+".f64"),dtype="<f8").reshape(count,16)
    if not np.isfinite(data[:,:8]).all() or np.any(data[:,0]<=0):
        raise ValueError("Probe did not complete all first passages")
    return np.diff(data[:,:8],axis=1).astype("<f8")


def solve(base, plan, probes):
    frozen = {}
    for name,spec in plan["cases"].items():
        directory = base/name
        directory.mkdir()
        train = directory/"training.f64"
        train.write_bytes(train_columns(name,plan["marker_count"]))
        local = directory/"probe_durations.f64"
        probe_durations(base,name,probes["count"]).tofile(local)
        cmd = [*spec["command"][:4],str(HERE/"solve.jl"),str(train),str(local),str(directory)]
        timing = execute(cmd,directory/"solve.log",timeout=120)
        write_json(directory/"execution.json",dict(command=cmd,**timing))
        frozen[name] = read_prediction(directory,spec)
    return frozen


def read_prediction(directory, spec):
    return dict(Pe=spec["Pe"],dt_s=spec["dt_s"],
                training_sha256=sha(directory/"training.f64"),
                probe_durations_sha256=sha(directory/"probe_durations.f64"),
                rule_sha256={name:sha(directory/name) for name in RULE_TABLES},
                **{name:np.loadtxt(directory/(name+".csv"),delimiter=",").tolist() for name in PREDICTORS})


def errors(observed, predicted):
    predicted = np.asarray(predicted)
    if predicted.shape != (6,3) or observed.shape != (6,2) or not np.isfinite(predicted).all():
        raise ValueError("Expected finite predictions for all six held-out sections")
    np.testing.assert_allclose(predicted[:,0],1.,rtol=0.,atol=1e-12)
    relative = (predicted[:,1:]-observed)/observed
    mean,variance = relative.T
    peak = int(np.argmax(abs(variance)))
    return dict(mean_relative_error=mean.tolist(),variance_relative_error=variance.tolist(),
                maximum_mean_error=float(max(abs(mean))),maximum_variance_error=float(abs(variance[peak])),
                variance_peak_section=peak+3,mean_passed=bool(max(abs(mean))<=LIMITS["mean"]),
                variance_passed=bool(abs(variance[peak])<=LIMITS["variance"]))


def score(base, summary, plan, frozen):
    cases = {}
    for name,predictions in frozen.items():
        data = archive_case(summary,plan,name)
        if (base/name/"training.f64").read_bytes() != data[:,:2].tobytes():
            raise ValueError("Training input differs from original sections 1–2")
        arrivals = data[:,2:8]
        if not np.isfinite(arrivals).all() or np.any(np.diff(data[:,:8],axis=1)<=0):
            raise ValueError("Incomplete or unordered holdout passages")
        observed = np.column_stack((arrivals.mean(axis=0),arrivals.var(axis=0)))
        if np.any(observed<=0):
            raise ValueError("Expected nonzero holdout arrival means and variances")
        cases[name] = dict(observed=observed.tolist(),
                           **{key:errors(observed,predictions[key]) for key in PREDICTORS})
    return dict(limits=LIMITS,cases=cases,passed=all(v["sectional"]["variance_passed"] for v in cases.values()))


def freeze_and_score(base, summary, plan, probes):
    if read_json(base/"probe_plan.json") != probes:
        raise ValueError("Fresh probes do not match the preregistered plan")
    sources = source_manifest([HERE/name for name in ("run.py","model.jl","solve.jl","test_model.jl","protocol.md")]+[
        HERE.parent/name for name in (*TRACKER,"disordered_exchange/model.jl","disordered_exchange/run.py",
                                      "dispersion_holdout/run.py","transport_accuracy/run.py")])
    predictions = solve(base,plan,probes)
    path = base/"frozen_prediction.json"
    write_json(path,dict(parent_sha256=PARENT_SHA,probe_plan_sha256=sha(base/"probe_plan.json"),
                         source_sha256=sources,cases=predictions))
    digest = sha(path)
    print("Predictions frozen before scoring sections 3–8",digest,flush=True)
    analysis = score(base,summary,plan,predictions)
    record = base/"summary.json"
    write_json(record,dict(schema="pore2process.disordered-section-exchange.v1",
        run_id=base.name,source_sha256=sources,parent_sha256=PARENT_SHA,frozen_prediction_sha256=digest,
        artifacts_sha256={p.relative_to(base).as_posix():sha(p) for p in sorted(base.rglob("*")) if p.is_file()},
        dry_carbon_mass_kg=summary["analysis"]["dry_carbon_mass_kg"],
        geometry_sha256=summary["analysis"]["geometry_sha256"],
        throughput_m3_s=plan["compact"]["throughput_m3_s"],analysis=analysis,
        new_native_flow_runs=0,new_adsorption_runs=0,figure7_fits=0))
    write_json(OUTPUT/"latest.json",dict(summary=str(record.relative_to(ROOT))))
    for name,case in analysis["cases"].items():
        print(name,{key:case[key]["maximum_variance_error"] for key in PREDICTORS})
    print(record)


def run():
    summary,plan = parent_inputs()
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")+"-"+uuid4().hex[:8]
    base = OUTPUT/"runs"/stamp
    base.mkdir(parents=True)
    probes = probe_plan(base,plan)
    print("Probing saved N160 flow",base,flush=True)
    verify_saved_field(summary,plan)
    run_probes(base,probes)
    freeze_and_score(base,summary,plan,probes)


def check(path):
    record = validate_record(path)
    summary,plan = parent_inputs()
    verify_saved_field(summary,plan)
    base = path.parent
    probes = probe_plan(base,plan)
    probe_path = base/"probe_plan.json"
    frozen_path = base/"frozen_prediction.json"
    predictions = {name:read_prediction(base/name,spec) for name,spec in plan["cases"].items()}
    for name in predictions:
        if (base/name/"probe_durations.f64").read_bytes() != probe_durations(base,name,probes["count"]).tobytes():
            raise ValueError(f"Probe duration inputs differ from the fresh tracer output: {name}")
    expected = dict(parent_sha256=PARENT_SHA,probe_plan_sha256=sha(probe_path),
                    source_sha256=record["source_sha256"],cases=predictions)
    if (record["parent_sha256"]!=PARENT_SHA or probes!=read_json(probe_path)
            or record["frozen_prediction_sha256"]!=sha(frozen_path) or expected!=read_json(frozen_path)
            or record["analysis"]!=score(base,summary,plan,predictions)):
        raise ValueError("Frozen section-specific prediction differs on replay")
    print(dict(sources=len(record["source_sha256"]),artifacts=len(record["artifacts_sha256"]),
               passed=record["analysis"]["passed"]))


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check",type=Path)
    args = parser.parse_args()
    check(args.check.resolve()) if args.check else run()
