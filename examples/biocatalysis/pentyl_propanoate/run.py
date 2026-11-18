"""Fit batch-only kinetics, then predict and compare the withheld Lorenzo bed runs."""

import argparse
import csv
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from pore2process_demo.bridge import runtime
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT


HERE = Path(__file__).resolve().parent
SOURCE = ROOT / "literature/sources"
BATCH = SOURCE / "lorenzo2026/batch_2to1.csv"
EQUILIBRIUM = SOURCE / "lorenzo2026/apparent_equilibrium_2to1.csv"
HOLDOUT = SOURCE / "lorenzo2025/conversion_2to1.csv"
PROTOCOL = HERE / "protocol.json"
MIXING = HERE / "mixing/reference.json"
OUTPUT = ROOT / "results/biocatalysis/pentyl_propanoate/runs"


def read_rows(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def write_rows(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def array(rows, key):
    return np.array([float(row[key]) for row in rows])


def error_metrics(errors):
    return {
        "count": len(errors),
        "rmse_percentage_points": float(100 * np.sqrt(np.mean(errors**2))),
        "mae_percentage_points": float(100 * np.mean(abs(errors))),
        "bias_percentage_points": float(100 * np.mean(errors)),
        "max_absolute_error_percentage_points": float(100 * np.max(abs(errors))),
    }


def fit_record(model, fit, temperature):
    return {
        "coefficients": np.array(model.coefficients(fit.kinetics)).tolist(),
        "mass_specific_k_mol_kg_s": float(model.rate_constant(fit.kinetics, temperature)),
        "apparent_Kx": float(model.equilibrium_constant(fit.kinetics, temperature)),
        "equilibrium_conversion": float(model.equilibrium_conversion(fit.kinetics, temperature)),
        "iterations": int(fit.iterations),
        "gradient_inf": float(fit.gradient_inf),
    }


def calibrate(model, plan, output):
    data = read_rows(BATCH)
    assays = array(data, "assay_id").astype(int)
    temperatures = array(data, "temperature_C")
    times = array(data, "time_min")
    loading = array(data, "catalyst_mass_pct_of_initial_acid")
    observed = array(data, "acid_conversion")
    apparent_k = {int(row["temperature_C"]): float(row["apparent_equilibrium_constant_Kx"])
                  for row in read_rows(EQUILIBRIUM)}
    if set(assays) != set(plan["batch_assays"]) or set(temperatures) != set(apparent_k):
        raise ValueError("The batch assay or apparent equilibrium series differs from the protocol.")

    pooled = model.fit_batch(temperatures, loading, times, observed)
    predicted = np.empty(len(data))
    local_kinetics = {}
    per_temperature = {}
    for temperature, k in apparent_k.items():
        mask = temperatures == temperature
        fit = model.fit_at_temperature(temperature, loading[mask], times[mask], observed[mask], k)
        local_kinetics[temperature] = fit.kinetics
        predicted[mask] = np.array(fit.predicted)
        per_temperature[str(temperature)] = fit_record(model, fit, temperature)

    pooled_prediction = np.array(pooled.predicted)
    residual = predicted - observed
    pooled_residual = pooled_prediction - observed
    after_start = times > 0
    report = {
        "sources": {str(path.relative_to(ROOT)): sha(path) for path in (PROTOCOL, BATCH, EQUILIBRIUM)},
        "primary_fit_by_temperature": per_temperature,
        "primary_fit": error_metrics(residual[after_start]),
        "primary_fit_by_assay": {
            str(assay): error_metrics(residual[(assays == assay) & after_start])
            for assay in plan["batch_assays"]
        },
        "pooled_fit_at_60_C": fit_record(model, pooled, 60),
        "pooled_fit": error_metrics(pooled_residual[after_start]),
    }
    write_json(output / "batch_fit.json", report)
    write_rows(output / "batch_fit.csv", [
        dict(row, predicted_conversion=primary, pooled_predicted_conversion=smooth)
        for row, primary, smooth in zip(data, predicted, pooled_prediction, strict=True)
    ])
    return local_kinetics, pooled.kinetics, report


def predict(model, plan, local_kinetics, pooled, mixing_path):
    mixing = read_json(mixing_path)
    variance = mixing["water_tracer"]["pulse_free_relative_variance"]
    peclet = float(model.axial_peclet(variance))
    primary_grid, fine_grid = plan["axial_dispersion_cell_counts"]
    catalyst_kg = plan["bed_catalyst_mass_kg"]
    predictions = []
    for temperature in plan["bed_controlled_temperatures_C"]:
        kinetics = local_kinetics[temperature]
        for flow in plan["bed_total_feed_flow_ml_min"]:
            acid_mol_s = plan["bed_feed_acid_mol_L"] * flow / 60_000
            plug = model.bed_prediction(kinetics, temperature, catalyst_kg, acid_mol_s)
            dispersed = model.bed_prediction(
                kinetics, temperature, catalyst_kg, acid_mol_s, peclet=peclet, cells=primary_grid,
            )
            fine = model.bed_prediction(
                kinetics, temperature, catalyst_kg, acid_mol_s, peclet=peclet, cells=fine_grid,
            )
            smooth = model.bed_prediction(pooled, temperature, catalyst_kg, acid_mol_s)
            predictions.append({
                "temperature_C": temperature,
                "flow_ml_min": flow,
                "acid_feed_mol_s": acid_mol_s,
                "plug_flow_conversion": float(plug.conversion),
                "water_mixing_scenario_conversion": float(fine.conversion),
                "pooled_kinetics_plug_flow_conversion": float(smooth.conversion),
                "axial_grid_change": float(abs(dispersed.conversion - fine.conversion)),
                "axial_balance_error": float(fine.balance_error),
            })
    maximum_grid_change = max(row["axial_grid_change"] for row in predictions)
    maximum_balance_error = max(row["axial_balance_error"] for row in predictions)
    if (maximum_grid_change > plan["max_grid_change_absolute_conversion"]
            or maximum_balance_error > plan["max_stoichiometric_balance_error"]):
        raise RuntimeError("The reacting-bed grid or acid balance check failed.")
    return predictions, {
        "mixing_source": str(mixing_path.relative_to(ROOT)),
        "mixing_source_sha256": sha(mixing_path),
        "water_tracer_pulse_free_relative_variance": variance,
        "aqueous_mixing_check": mixing["finite_pulse_comparison"],
        "aqueous_mixing_benchmark_passed": mixing["water_tracer_benchmark_passed"],
        "scenario_axial_peclet": peclet,
        "max_grid_change_absolute_conversion": maximum_grid_change,
        "max_stoichiometric_balance_error": maximum_balance_error,
    }


def compare(predictions):
    def condition(row):
        return float(row["temperature_C"]), float(row["flow_ml_min"])

    held_out = read_rows(HOLDOUT)
    observations = {condition(row): float(row["acid_conversion"]) for row in held_out}
    if len(observations) != len(held_out) or set(observations) != set(map(condition, predictions)):
        raise ValueError("The held-out figures contain duplicate or unmatched temperature/flow conditions.")
    comparison = [dict(row, observed_acid_conversion=observations[condition(row)])
                  for row in predictions]
    observed = array(comparison, "observed_acid_conversion")
    if not np.all(np.isfinite(observed) & (observed >= 0) & (observed <= 1)):
        raise ValueError("Invalid experimental conversion.")
    metrics = {
        name: error_metrics(array(comparison, name) - observed)
        for name in ("plug_flow_conversion", "water_mixing_scenario_conversion",
                     "pooled_kinetics_plug_flow_conversion")
    }
    return comparison, metrics


def plot_comparison(rows, path):
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    figure = Figure(figsize=(9, 5), layout="constrained")
    FigureCanvasAgg(figure)
    axes = figure.subplots(2, 3, sharex=True, sharey=True).ravel()
    for axis, temperature in zip(axes[:5], sorted({r["temperature_C"] for r in rows}), strict=True):
        cases = [row for row in rows if row["temperature_C"] == temperature]
        flow = array(cases, "flow_ml_min")
        axis.plot(flow, 100 * array(cases, "plug_flow_conversion"), color="#1776a8",
                  linewidth=2, label="Batch kinetics, plug flow")
        axis.plot(flow, 100 * array(cases, "water_mixing_scenario_conversion"), color="#a84b17",
                  linewidth=1.2, linestyle="--", label="OpenLB water-mixing scenario")
        axis.scatter(flow, 100 * array(cases, "observed_acid_conversion"), color="#222222",
                     s=21, zorder=3, label="Held-out measurements")
        axis.set_title(f"{temperature} °C")
        axis.tick_params(labelbottom=True)
        axis.set_xticks([0.7, 1, 2, 3], ["0.7", "1", "2", "3"])
        axis.set_ylim(30, 100)
        axis.grid(axis="y", alpha=0.23)
    handles, labels = axes[0].get_legend_handles_labels()
    axes[-1].set_axis_off()
    axes[-1].legend(handles, labels, loc="center left", frameon=False, fontsize=9)
    figure.supxlabel("Liquid feed flow (mL/min)")
    figure.supylabel("Propanoic acid conversion (%)")
    figure.suptitle("Pentyl propanoate at a 2:1 alcohol-to-acid feed ratio")
    figure.savefig(path, dpi=150)


def run(output, mixing_summary=MIXING, *, batch_only=False):
    output = output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    plan = read_json(PROTOCOL)
    julia = runtime()
    julia.include(str(HERE / "model.jl"))
    model = julia.PentylPropanoate
    local_kinetics, pooled, batch = calibrate(model, plan, output)
    if batch_only:
        return batch
    predictions, numerical = predict(model, plan, local_kinetics, pooled, mixing_summary.resolve())
    write_rows(output / "bed_predictions_before_comparison.csv", predictions)
    comparison, metrics = compare(predictions)
    write_rows(output / "bed_comparison.csv", comparison)
    plot_comparison(comparison, output / "bed_comparison.png")
    result = {
        "julia_version": str(julia.seval("VERSION")),
        "model_sha256": sha(HERE / "model.jl"),
        "runner_sha256": sha(HERE / "run.py"),
        "held_out_conversion_source": str(HOLDOUT.relative_to(ROOT)),
        "held_out_conversion_sha256": sha(HOLDOUT),
        "bed_comparison_sha256": sha(output / "bed_comparison.csv"),
        "batch": batch,
        "bed": numerical | {"comparison": metrics},
    }
    write_json(output / "summary.json", result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    parser.add_argument("--output", type=Path, default=OUTPUT / stamp)
    parser.add_argument("--mixing-summary", type=Path, default=MIXING)
    parser.add_argument("--batch-only", action="store_true")
    args = parser.parse_args()
    report = run(args.output, args.mixing_summary, batch_only=args.batch_only)
    print("Report:", args.output / ("batch_fit.json" if args.batch_only else "summary.json"))
    print("Batch primary fit:" if args.batch_only else "Held-out conversion:",
          report["primary_fit"] if args.batch_only else report["bed"]["comparison"])


if __name__ == "__main__":
    main()
