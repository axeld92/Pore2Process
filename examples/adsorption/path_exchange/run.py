"""Measure path-time exchange from the archived TRT flow and score a frozen prediction."""
import argparse
from datetime import datetime, timezone
from pathlib import Path
import shutil
from uuid import uuid4

import numpy as np

from examples.adsorption.dispersion_holdout.run import DT, PARENT, PARENT_SHA, execute, inputs, repeat_mesh
from examples.adsorption.porous_outlet_tracer.run import PULSE_S, endpoint, history
from examples.adsorption.residence_time_holdout.run import metrics
from examples.adsorption.stationary_curved_wall import reference as ref
from examples.adsorption.transport_accuracy.run import validate_record
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT
from pore2process_demo.provenance import source_manifest

from .predict import compare, pair_prediction, summary

HERE = Path(__file__).parent
OUTPUT = ROOT / "results/adsorption/path_exchange"
ARCHIVE = ROOT / "results/adsorption/residence_time_holdout/runs/20261112T044839Z-2c87f73c"
PATHS_SHA = "8158092b17c178d9621d981bb3f8fc64335267b85a232cb05189acaad1680398"
BENCHMARK_SHA = "c42826da68b23b15300e5cb2bf74d7506096b8672130808114b297639062522d"
HISTORY_SHA = "4413ca03a19bbae95b3237d1bee946ffd74fc43901c3b5facc20d8617504e78f"
COUNT = 65536
STOP_S = 8.
SEED = 3526139
SETTINGS = dict(primary=(ref.D, DT, SEED), half_step=(ref.D, DT/2, SEED),
                replicate=(ref.D, DT, SEED + COUNT), advection=(0., DT, SEED))


def commands(base, runtime, mesh):
    dx = mesh["dx_m"] * 3
    return {name: [runtime[0], "--threads=10", *runtime[1:], str(HERE / "track.jl"),
                   *map(str, (base / "compact", base / (name + ".f64"), dx, d, dt, STOP_S, COUNT, seed))]
            for name, (d, dt, seed) in SETTINGS.items()}


def one_cell_case(times, raw_times, expected, stats):
    first = np.sort(raw_times[np.isfinite(raw_times[:, 0]), 0])
    weight = np.full(len(first), 1 / len(first))
    # Use the existing full comparison window to compute 99% arrival,
    # while compare() uses the measured one-cell window for CDF error.
    curve, values = summary(times, first, weight)
    report = compare(expected, stats, curve, values)
    names = ("t05", "rise_10_to_75", "mean", "variance", "t99")
    return dict(metrics=values, comparison=report, qualified=all(report["passes"][k] for k in names))


def prediction(base, mesh, short_steps):
    if sha(ARCHIVE / "paths.csv") != PATHS_SHA:
        raise ValueError("Archived one-cell path distribution changed")
    delays, weights = np.loadtxt(ARCHIVE / "paths.csv", delimiter=",", skiprows=1, unpack=True)
    times = DT * np.arange(2 * short_steps + 1)
    dose = mesh["flow_m3_s"] * PULSE_S
    transit = mesh["volume_m3"] / (2 * mesh["flow_m3_s"])
    measured = history(PARENT.parent / "N384/history.csv")
    cdf = measured["recovered_kg"] / dose
    stats = endpoint(measured, transit)
    one = dict(**metrics(measured["time_s"], cdf), mean=stats["mean_exit_s"], variance=stats["variance_s2"])
    baseline, baseline_stats = summary(times, 2 * delays, weights)
    curves = {"persistent": baseline}
    result = {"persistent": baseline_stats}
    for name in SETTINGS:
        trajectories = np.fromfile(base / (name + ".f64"), dtype="<f8").reshape(-1, 2)
        if len(trajectories) != COUNT:
            raise ValueError("Marker output does not match the requested count")
        curves[name], result[name] = pair_prediction(trajectories, delays, weights, times)
        result[name]["one_cell_raw_check"] = one_cell_case(measured["time_s"], trajectories, cdf, one)
    for name in ("half_step", "replicate"):
        a, b = result["primary"], result[name]
        limits = dict(cdf=.005, variance=.05, t99=.05)
        errors = dict(cdf=float(max(abs(curves[name] - curves["primary"]))),
                      variance=abs(b["variance"] / a["variance"] - 1),
                      t99=abs(b["t99"] / a["t99"] - 1))
        result[name]["sensitivity"] = dict(errors=errors, limits=limits,
                                           passed=all(v <= limits[k] for k, v in errors.items()))
    return times, curves, result


def score(base, mesh, result):
    if sha(ARCHIVE / "summary.json") != BENCHMARK_SHA or sha(ARCHIVE / "holdout/history.csv") != HISTORY_SHA:
        raise ValueError("Archived two-cell benchmark changed")
    data = history(ARCHIVE / "holdout/history.csv")
    saved = history(base / "prediction.csv")
    np.testing.assert_array_equal(data["time_s"], saved["time_s"])
    cdf = data["recovered_kg"] / (mesh["flow_m3_s"] * PULSE_S)
    stats = endpoint(data, mesh["volume_m3"] / mesh["flow_m3_s"])
    if not stats["tail_passed"]:
        raise ValueError("Incomplete benchmark tail")
    measured = dict(**metrics(data["time_s"], cdf), mean=stats["mean_exit_s"], variance=stats["variance_s2"])
    comparisons = {name: compare(cdf, measured, saved[name], stats) for name, stats in result.items()}
    primary = comparisons["primary"]
    numerical = all(result[k]["sensitivity"]["passed"] for k in ("half_step", "replicate"))
    one_cell = result["primary"]["one_cell_raw_check"]["qualified"]
    return dict(measured=measured, comparisons=comparisons, primary_numerically_stable=numerical,
                primary_one_cell_qualified=one_cell, primary_prediction_passed=all(primary["passes"].values()),
                qualified=numerical and one_cell and all(primary["passes"].values()))


def run():
    _, runtime, info, coarse, steps = inputs()
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ") + "-" + uuid4().hex[:8]
    base = OUTPUT / "runs" / stamp
    base.mkdir(parents=True)
    files = [HERE / name for name in ("run.py", "predict.py", "track.jl", "test_track.jl", "protocol.md")]
    files += [HERE.parent / p for p in ("dispersion_holdout/run.py", "residence_time_holdout/run.py",
                                       "residence_time_holdout/paths.py")]
    hashes = source_manifest(files)
    for name in hashes:
        dest = base / "source_snapshot" / name
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT / name, dest)
    mesh = repeat_mesh(base / "compact", info, coarse, write=True)
    req = commands(base, runtime, mesh)
    record = dict(schema="pore2process.path-exchange.v1", run_id=stamp, source_sha256=hashes,
                  mesh=mesh, requests=req, one_cell_sha256=PARENT_SHA, paths_sha256=PATHS_SHA,
                  benchmark_sha256=BENCHMARK_SHA, new_openlb_runs=0, new_full_bed_runs=0, new_dwsim_runs=0)
    write_json(base / "plan.json", record)
    print("Marker study", base, flush=True)
    for name, command in req.items():
        duration = execute(command, base / (name + ".log"), timeout=3600)
        write_json(base / (name + ".json"), duration)
        print(name, duration, flush=True)
    times, curves, result = prediction(base, mesh, steps)
    np.savetxt(base / "prediction.csv", np.column_stack([times, *curves.values()]), delimiter=",",
               header=",".join(("time_s", *curves)), comments="")
    write_json(base / "frozen_prediction.json", dict(curve_sha256=sha(base / "prediction.csv"), cases=result))
    record["frozen_prediction_sha256"] = sha(base / "frozen_prediction.json")
    print("Predictions frozen", record["frozen_prediction_sha256"], flush=True)
    record["analysis"] = score(base, mesh, result)
    record["artifacts_sha256"] = {str(p.relative_to(base)): sha(p) for p in base.rglob("*") if p.is_file()}
    path = base / "summary.json"
    write_json(path, record)
    write_json(OUTPUT / "latest.json", dict(summary=str(path.relative_to(ROOT))))
    print("Complete", path, record["analysis"], flush=True)


def check(path):
    record = validate_record(path)
    _, runtime, info, coarse, steps = inputs()
    base = path.parent
    mesh = repeat_mesh(base / "compact", info, coarse, write=False)
    if mesh != record["mesh"] or commands(base, runtime, mesh) != record["requests"]:
        raise ValueError("Marker case input or command differs")
    times, curves, result = prediction(base, mesh, steps)
    frozen = read_json(base / "frozen_prediction.json")
    saved = history(base / "prediction.csv")
    np.testing.assert_array_equal(saved["time_s"], times)
    for name, curve in curves.items():
        np.testing.assert_array_equal(saved[name], curve)
    if (sha(base / "frozen_prediction.json") != record["frozen_prediction_sha256"]
            or frozen != dict(curve_sha256=sha(base / "prediction.csv"), cases=result)
            or score(base, mesh, result) != record["analysis"]):
        raise ValueError("Frozen predictions or held-out comparison differ")
    print(dict(sources=len(record["source_sha256"]), artifacts=len(record["artifacts_sha256"]),
               qualified=record["analysis"]["qualified"]))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", type=Path)
    args = parser.parse_args()
    check(args.check.resolve()) if args.check else run()
