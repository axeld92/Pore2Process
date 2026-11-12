"""Freeze a one-cell path model, then test a complete two-cell pore pulse."""
import argparse
from datetime import datetime, timezone
from pathlib import Path
import shutil
from uuid import uuid4

import numpy as np

from examples.adsorption.dispersion_holdout.run import DT, PARENT, PARENT_SHA, execute, inputs, repeat_mesh
from examples.adsorption.porous_outlet_tracer.run import PULSE_S, endpoint, history
from examples.adsorption.stationary_curved_wall import reference as ref
from examples.adsorption.transport_accuracy.run import validate_record
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT
from pore2process_demo.provenance import source_manifest

from .paths import moments, paths, pulse_cdf

HERE = Path(__file__).parent
OUTPUT = ROOT / "results/adsorption/residence_time_holdout"
BRICK = HERE.parent / "porous_outlet_tracer/solve_brick.jl"
PREFIX = ROOT / "results/adsorption/dispersion_holdout/runs/20261112T033234Z-32e27f5e"
PREFIX_SHA = "6b6ea7f8c84a3106c40711e11f76895f533ce495274d3ea8a3153d67c31399a7"
TAIL_S = 0.60
LIMITS = dict(profile_rmse=.030, mean=.05, variance=.10, t05=.05,
              rise_10_to_75=.10, t99=.10, tail_mass=.25)


def metrics(times, cumulative):
    quantiles = {f"t{int(100*f):02}": float(np.interp(f, cumulative, times)) if cumulative[-1] >= f
                 else None for f in (.05, .10, .50, .75, .90, .95, .99)}
    return dict(**quantiles, rise_10_to_75=quantiles["t75"] - quantiles["t10"],
                tail_mass=float(1 - np.interp(TAIL_S, times, cumulative)))


def calibration(mesh, steps):
    data = history(PARENT.parent / "N384/history.csv")
    dose = mesh["flow_m3_s"] * PULSE_S
    cumulative = data["recovered_kg"] / dose
    transit = mesh["volume_m3"] / (2 * mesh["flow_m3_s"])
    delays, weights, adjustment = paths(cumulative, DT, PULSE_S, transit)
    times = DT * np.arange(steps + 1)
    short = pulse_cdf(data["time_s"], delays, weights, PULSE_S)
    long = pulse_cdf(times, 2 * delays, weights, PULSE_S)
    short_error = float(max(abs(short - cumulative)))
    if short_error > .005 or abs(weights.sum() - 1) > 1e-12 or abs(np.dot(delays, weights) / transit - 1) > 1e-12:
        raise ValueError(f"Calibration failed its one-cell mass, mean, or shape check (max CDF error {short_error})")
    mean, variance = moments(2 * delays, weights, PULSE_S)
    if 2 * delays[-1] + PULSE_S > times[-1]:
        raise ValueError("Frozen path distribution extends beyond the planned holdout")
    metadata = dict(calibration_parent_sha256=PARENT_SHA, one_cell_max_cdf_error=short_error,
                    one_cell_recovered=float(cumulative[-1]), slow_path_adjustment_s=adjustment,
                    one_cell_transit_s=transit, holdout_predicted_mean_s=mean,
                    holdout_predicted_variance_s2=variance, predicted=metrics(times, long))
    return metadata, {"paths.csv": (delays, weights), "one_cell.csv": (data["time_s"], short),
                      "prediction.csv": (times, long)}


def request(base, runtime, mesh, steps):
    return dict(command=[runtime[0], "--threads=10", *runtime[1:], str(BRICK), "solve",
                         *map(str, (base / "compact", base / "holdout", 3, mesh["dx_m"], ref.D, DT, steps))],
                steps=steps, dose_kg=mesh["flow_m3_s"] * PULSE_S)


def verify_frozen(base, mesh, req, frozen_hash):
    manifest = base / "frozen_prediction.json"
    if sha(manifest) != frozen_hash:
        raise ValueError("Frozen manifest changed")
    frozen = read_json(manifest)
    model, curves = calibration(mesh, req["steps"])
    if model != frozen["model"]:
        raise ValueError("Frozen model differs from the one-cell-only calibration")
    for name, expected in curves.items():
        if sha(base / name) != frozen["files"][name]:
            raise ValueError(f"Frozen curve changed: {name}")
        np.testing.assert_array_equal(np.loadtxt(base / name, delimiter=",", skiprows=1, unpack=True), expected)
    return model


def analyze(base, mesh, req, model):
    h = history(base / "holdout/history.csv")
    times, predicted = np.loadtxt(base / "prediction.csv", delimiter=",", skiprows=1, unpack=True)
    dose = req["dose_kg"]
    np.testing.assert_array_equal(h["time_s"], times)
    np.testing.assert_allclose(h["injected_kg"], mesh["flow_m3_s"] * np.minimum(times, PULSE_S),
                               rtol=1e-12, atol=0)
    if not np.isfinite(h.view(np.float64)).all():
        raise ValueError("Non-finite tracer history")
    transfer = np.diff(h["recovered_kg"], prepend=0.) / dose
    np.testing.assert_allclose(transfer, h["interval_exit_kg"] / dose, rtol=0, atol=1e-14)
    actual = h["recovered_kg"] / dose
    mass = float(max(abs(h["liquid_kg"] + h["recovered_kg"] - h["injected_kg"])) / dose)
    lo, hi = float(min(h["min_kg_m3"])), float(max(h["max_kg_m3"]))
    transit = mesh["volume_m3"] / mesh["flow_m3_s"]
    stats = endpoint(h, transit)
    if sha(PREFIX / "summary.json") != PREFIX_SHA:
        raise ValueError("Archived early comparison changed")
    prefix = read_json(PREFIX / "summary.json")
    if sha(PREFIX / "holdout/history.csv") != prefix["artifacts_sha256"]["holdout/history.csv"]:
        raise ValueError("Archived early history changed")
    old = history(PREFIX / "holdout/history.csv")
    np.testing.assert_array_equal(h[:len(old)], old)
    if (base / "holdout/final.f64").stat().st_size != mesh["fine_cells"] * 8:
        raise ValueError("Complete final tracer field required")
    measured = metrics(times, actual)
    estimate = model["predicted"]
    errors = dict(profile_rmse=float(np.sqrt(np.mean((actual - predicted)**2))),
                  mean=abs(model["holdout_predicted_mean_s"] / stats["mean_exit_s"] - 1),
                  variance=abs(model["holdout_predicted_variance_s2"] / stats["variance_s2"] - 1))
    for name in ("t05", "rise_10_to_75", "t99", "tail_mass"):
        value = measured[name]
        errors[name] = abs(estimate[name] / value - 1) if value is not None and value > 0 else None
    checks = {name: value is not None and value <= LIMITS[name] for name, value in errors.items()}
    checks.update(mass_conservation=mass <= 1e-8, bounds=lo >= -1e-8 and hi <= 1+1e-8,
                  nonnegative_outlet=bool(min(transfer) >= -1e-14), complete_recovery=stats["tail_passed"])
    return dict(samples=len(h), reproduced_prefix_samples=len(old), pore_stats=stats,
                measured=measured, predicted=estimate, errors=errors, checks=checks,
                limits=LIMITS, qualified=all(checks.values()), max_mass_residual=mass,
                concentration_range_kg_m3=[lo, hi])


def run():
    _, runtime, info, coarse, short_steps = inputs()
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ") + "-" + uuid4().hex[:8]
    base = OUTPUT / "runs" / stamp
    base.mkdir(parents=True)
    dependencies = [HERE / name for name in ("run.py", "paths.py", "protocol.md")]
    dependencies += [BRICK, BRICK.with_name("solve.jl"), HERE.parent / "dispersion_holdout/run.py"]
    dependencies += [ROOT / "julia/adsorption" / name for name in
                     ("BrickOpenTracer.jl", "BoundedFaceTransport.jl", "ThreadedBoundedFaceTransport.jl",
                      "OpenBoundedFaceTransport.jl")]
    hashes = source_manifest(dependencies)
    for name in hashes:
        target = base / "source_snapshot" / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT / name, target)
    mesh = repeat_mesh(base / "compact", info, coarse, write=True)
    req = request(base, runtime, mesh, 2 * short_steps)
    record = dict(schema="pore2process.residence-time-holdout.v1", run_id=stamp, source_sha256=hashes,
                  mesh=mesh, request=req, new_OpenLB_runs=0, new_full_bed_runs=0, new_DWSIM_runs=0)
    write_json(base / "plan.json", record)
    model, curves = calibration(mesh, req["steps"])
    headers = dict(zip(curves, ("delay_s,mass_fraction", "time_s,fraction_recovered", "time_s,fraction_recovered")))
    for name, data in curves.items():
        np.savetxt(base / name, np.column_stack(data), delimiter=",", header=headers[name], comments="")
    frozen = base / "frozen_prediction.json"
    write_json(frozen, dict(model=model, files={name: sha(base / name) for name in curves}))
    frozen_hash = sha(frozen)
    write_json(base / "launch.json", dict(frozen_prediction_sha256=frozen_hash))
    print("Frozen before pore run", base, model, flush=True)
    (base / "holdout").mkdir()
    timing = execute(req["command"], base / "holdout/stdout.log", timeout=28800)
    write_json(base / "holdout/execution.json", timing)
    verify_frozen(base, mesh, req, frozen_hash)
    record.update(frozen_prediction_sha256=frozen_hash, analysis=analyze(base, mesh, req, model))
    record["artifacts_sha256"] = {str(p.relative_to(base)): sha(p) for p in base.rglob("*") if p.is_file()}
    summary = base / "summary.json"
    write_json(summary, record)
    write_json(OUTPUT / "latest.json", dict(summary=str(summary.relative_to(ROOT))))
    print("Completed residence time holdout", summary, record["analysis"], flush=True)


def check(path):
    record = validate_record(path)
    _, runtime, info, coarse, short_steps = inputs()
    base = path.parent
    mesh = repeat_mesh(base / "compact", info, coarse, write=False)
    req = request(base, runtime, mesh, 2 * short_steps)
    model = verify_frozen(base, mesh, req, record["frozen_prediction_sha256"])
    if mesh != record["mesh"] or req != record["request"] or analyze(base, mesh, req, model) != record["analysis"]:
        raise ValueError("Holdout comparison differs on replay")
    print(dict(summary=str(path), sources=len(record["source_sha256"]),
               artifacts=len(record["artifacts_sha256"]), qualified=record["analysis"]["qualified"]))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", type=Path)
    args = parser.parse_args()
    check(args.check.resolve()) if args.check else run()
