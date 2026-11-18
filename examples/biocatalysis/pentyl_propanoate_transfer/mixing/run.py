"""Test the measured water-tracer spread using frozen batch kinetics."""

import argparse
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from examples.biocatalysis.pentyl_propanoate.run import (
    HERE as BASE,
    HOLDOUT,
    PROTOCOL,
    array,
    error_metrics,
    read_rows,
    write_rows,
)
from pore2process_demo.bridge import runtime
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT


TRACER = HOLDOUT.with_name("water_tracer.csv")
EXTRACTION = HOLDOUT.with_name("extraction.json")
MODEL = BASE / "model.jl"
FROZEN = BASE / "reference/summary.json"
OUTPUT = ROOT / "results/biocatalysis/pentyl_propanoate_transfer/mixing"


def water_outlet_moments():
    rows = [row for row in read_rows(TRACER) if row["tracer"] == "blue_dextran"]
    time = array(rows, "time_min")
    concentration = array(rows, "concentration_mg_ml")
    if (len(time) < 2 or not np.all(np.isfinite(time)) or not np.all(np.diff(time) > 0)
            or not np.all(np.isfinite(concentration) & (concentration >= 0))):
        raise ValueError("Expected increasing tracer times and nonnegative concentrations.")
    area = float(np.trapezoid(concentration, x=time))
    if not area > 0:
        raise ValueError("The water-tracer curve must have positive area.")
    mean = float(np.trapezoid(time * concentration, x=time) / area)
    variance = float(np.trapezoid((time - mean)**2 * concentration, x=time) / area)
    return {"marker_count": len(rows), "area_mg_min_ml": area,
            "mean_min": mean, "variance_min2": variance}


def remove_pulse(model, outlet, inlet_mean, inlet_variance):
    mean = outlet["mean_min"] - inlet_mean
    variance = outlet["variance_min2"] - inlet_variance
    if not mean > 0 or not variance > 0:
        raise ValueError("The assumed inlet pulse leaves no positive bed mean or variance.")
    relative = variance / mean**2
    return {
        "inlet_mean_min": inlet_mean,
        "inlet_variance_min2": inlet_variance,
        "bed_mean_min": mean,
        "bed_variance_min2": variance,
        "bed_relative_variance": relative,
        "danckwerts_peclet": float(model.axial_peclet(relative)),
    }


def predict(model, protocol, frozen_kinetics, scenarios):
    mass = protocol["bed_catalyst_mass_kg"]
    concentration = protocol["bed_feed_acid_mol_L"]
    coarse, fine = protocol["axial_dispersion_cell_counts"]
    rows = []
    for temperature in protocol["bed_controlled_temperatures_C"]:
        kinetics = model.Kinetics(*frozen_kinetics[str(temperature)]["coefficients"])
        for flow in protocol["bed_total_feed_flow_ml_min"]:
            args = (kinetics, temperature, mass, concentration * flow / 60_000)
            plug = model.bed_prediction(*args)
            row = {"temperature_C": temperature, "flow_ml_min": flow,
                   "plug_flow": float(plug.conversion)}
            grid_errors, balance_errors = [], []
            for name, scenario in scenarios.items():
                peclet = scenario["danckwerts_peclet"]
                lower_grid = model.bed_prediction(*args, peclet=peclet, cells=coarse)
                upper_grid = model.bed_prediction(*args, peclet=peclet, cells=fine)
                row[name] = float(upper_grid.conversion)
                grid_errors.append(float(abs(lower_grid.conversion - upper_grid.conversion)))
                balance_errors.extend((float(lower_grid.balance_error), float(upper_grid.balance_error)))
            row["max_grid_change_absolute"] = max(grid_errors)
            row["max_acid_balance_error"] = max(balance_errors)
            if (not all(0 <= row[name] <= 1 for name in ("plug_flow", *scenarios))
                    or row["max_grid_change_absolute"] > protocol["max_grid_change_absolute_conversion"]
                    or row["max_acid_balance_error"] > protocol["max_stoichiometric_balance_error"]):
                raise RuntimeError(f"Invalid or unresolved dispersion prediction at {temperature} C, {flow} mL/min.")
            rows.append(row)
    return rows


def compare(predictions, scenarios):
    def condition(row):
        return float(row["temperature_C"]), float(row["flow_ml_min"])

    held_out = read_rows(HOLDOUT)
    measured = {condition(row): float(row["acid_conversion"]) for row in held_out}
    if len(measured) != len(held_out) or set(measured) != set(map(condition, predictions)):
        raise ValueError("The held-out conversion conditions do not match the frozen protocol.")
    rows = [dict(row, observed_acid_conversion=measured[condition(row)]) for row in predictions]
    observed = array(rows, "observed_acid_conversion")
    if not np.all((observed >= 0) & (observed <= 1)):
        raise ValueError("Expected finite conversion fractions in the held-out dataset.")
    plug = array(rows, "plug_flow")
    metrics = {}
    for name in ("plug_flow", *scenarios):
        conversion = array(rows, name)
        error = conversion - observed
        metrics[name] = error_metrics(error) | {
            "overpredicted_count": int(np.count_nonzero(error > 0)),
            "minimum_signed_error_percentage_points": float(100 * np.min(error)),
            "maximum_reduction_from_plug_flow_percentage_points": float(100 * np.max(plug - conversion)),
        }
    return rows, metrics


def run(output):
    frozen = read_json(FROZEN)
    protocol_hash = frozen["batch"]["sources"][str(PROTOCOL.relative_to(ROOT))]
    if sha(MODEL) != frozen["model_sha256"] or sha(PROTOCOL) != protocol_hash:
        raise ValueError("The batch model or protocol differs from the frozen benchmark.")
    protocol = read_json(PROTOCOL)
    experiment = read_json(EXTRACTION)["water_tracer"]
    pulse = experiment["pulse_volume_ml"] / experiment["water_flow_ml_min"]
    outlet = water_outlet_moments()

    julia = runtime()
    if not julia.seval("isdefined(@__MODULE__, :PentylPropanoate)"):
        julia.include(str(MODEL))
    model = julia.PentylPropanoate
    scenarios = {
        "nominal_rectangular_pulse": remove_pulse(model, outlet, pulse / 2, pulse**2 / 12),
        "impulse_at_start": remove_pulse(model, outlet, 0.0, 0.0),
        # The latest possible time and zero inlet variance maximize the residual
        # variance divided by squared mean for any inlet supported on [0, pulse].
        "maximum_spreading_bound": remove_pulse(model, outlet, pulse, 0.0),
    }
    table = experiment["reported_rtd_moments"]["blue_dextran"]
    table_check = remove_pulse(model, table, pulse / 2, pulse**2 / 12)
    predictions = predict(model, protocol, frozen["batch"]["primary_fit_by_temperature"], scenarios)
    output.mkdir(parents=True, exist_ok=False)
    write_rows(output / "bed_predictions_before_comparison.csv", predictions)
    rows, metrics = compare(predictions, scenarios)
    write_rows(output / "bed_comparison.csv", rows)
    record = {
        "julia_version": str(julia.seval("VERSION")),
        "source_hashes": {str(path.relative_to(ROOT)): sha(path)
                          for path in (TRACER, EXTRACTION, MODEL, FROZEN, PROTOCOL, HOLDOUT)},
        "bed_comparison_sha256": sha(output / "bed_comparison.csv"),
        "fixed_catalyst_kg": protocol["bed_catalyst_mass_kg"],
        "fixed_inlet_acid_mol_l": protocol["bed_feed_acid_mol_L"],
        "water_flow_ml_min": experiment["water_flow_ml_min"],
        "nominal_pulse_duration_min": pulse,
        "figure_a1a_blue_dextran": outlet,
        "table_a1_blue_dextran": table,
        "table_a1_rectangular_pulse_check": table_check,
        "conditional_danckwerts_transfer": scenarios,
        "grid_cells": protocol["axial_dispersion_cell_counts"],
        "max_grid_change_absolute_conversion": max(row["max_grid_change_absolute"] for row in rows),
        "max_stoichiometric_balance_error": max(row["max_acid_balance_error"] for row in rows),
        "comparison": metrics,
    }
    write_json(output / "summary.json", record)
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    parser.add_argument("--output", type=Path, default=OUTPUT / stamp)
    args = parser.parse_args()
    record = run(args.output.resolve())
    print("Report:", args.output / "summary.json")
    for name, metrics in record["comparison"].items():
        print(f"{name}: RMSE {metrics['rmse_percentage_points']:.2f} pp, "
              f"{metrics['overpredicted_count']}/{metrics['count']} overpredictions")


if __name__ == "__main__":
    main()
