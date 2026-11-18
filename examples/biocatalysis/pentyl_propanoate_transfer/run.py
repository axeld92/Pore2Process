"""Transfer measured batch progress at fixed catalyst exposure into the bed."""

from argparse import ArgumentParser
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from examples.biocatalysis.pentyl_propanoate.run import (
    BATCH, HOLDOUT, PROTOCOL, error_metrics, read_rows, write_rows,
)
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT


HERE = Path(__file__).resolve().parent
FIGURE_S3 = ROOT / "literature/sources/lorenzo2025/batch_2to1_markers.csv"
ORIGINAL_COMPARISON = PROTOCOL.parent / "reference/bed_comparison.csv"
OUTPUT = ROOT / "results/biocatalysis/pentyl_propanoate_transfer/rate"
MOLAR_MASS_ACID_KG_MOL = 0.0740785
CATALYST_FRACTION_2026 = 0.10
CATALYST_CONCENTRATION_2025_KG_L = 0.0253


def batch_series(rows, temperature, *, from_acid=False):
    """Use one assay and preserve the reported time order and conversion fluctuations."""
    selected = sorted((row for row in rows if int(row["temperature_C"]) == temperature),
                      key=lambda row: float(row["time_min"]))
    times = np.array([float(row["time_min"]) for row in selected])
    if (len(times) < 2 or times[0] != 0 or not np.all(np.isfinite(times))
            or np.any(np.diff(times) <= 0)):
        raise ValueError(f"Expected a strictly increasing batch time series at {temperature} °C.")
    if from_acid:
        acid = np.array([float(row["x_acid"]) for row in selected])
        conversion = 1 - acid / acid[0]
    else:
        conversion = np.array([float(row["acid_conversion"]) for row in selected])
    if conversion[0] != 0 or not np.all((conversion >= 0) & (conversion <= 1)):
        raise ValueError(f"Invalid batch conversion at {temperature} °C.")
    return times, conversion


def at_exposure(loading_kg_per_mol, mass_kg, acid_feed_mol_min, series):
    """Return the available batch conversion at the equivalent plug-flow exposure."""
    time_min = (mass_kg / acid_feed_mol_min) / loading_kg_per_mol
    times, conversion = series
    if not times[0] <= time_min <= times[-1]:
        return time_min, None
    return time_min, float(np.interp(time_min, times, conversion))


def make_predictions(plan):
    batch_2026 = [
        row for row in read_rows(BATCH)
        if float(row["catalyst_mass_pct_of_initial_acid"]) == 100 * CATALYST_FRACTION_2026
    ]
    batch_2025 = read_rows(FIGURE_S3)
    feed_mol_L = plan["bed_feed_acid_mol_L"]
    mass_kg = plan["bed_catalyst_mass_kg"]
    # Table S2 concentration basis for 2025, weighed acid-mass ratio for 2026.
    loading_2025 = CATALYST_CONCENTRATION_2025_KG_L / feed_mol_L
    loading_2026 = CATALYST_FRACTION_2026 * MOLAR_MASS_ACID_KG_MOL
    predictions = []
    for temperature in plan["bed_controlled_temperatures_C"]:
        printed = batch_series(batch_2026, temperature)
        acid_based = batch_series(batch_2026, temperature, from_acid=True)
        figure_s3 = batch_series(batch_2025, temperature)
        for flow in plan["bed_total_feed_flow_ml_min"]:
            acid_mol_min = feed_mol_L * flow / 1_000
            exposure = mass_kg / acid_mol_min
            time_2026, conversion_2026 = at_exposure(loading_2026, mass_kg, acid_mol_min, printed)
            _, depletion = at_exposure(loading_2026, mass_kg, acid_mol_min, acid_based)
            time_2025, conversion_2025 = at_exposure(loading_2025, mass_kg, acid_mol_min, figure_s3)
            predictions.append({
                "temperature_C": temperature,
                "flow_ml_min": flow,
                "catalyst_mass_kg": mass_kg,
                "acid_feed_mol_min": acid_mol_min,
                "exposure_kg_min_per_mol_acid": exposure,
                "equivalent_2025_batch_time_min": time_2025,
                "last_digitized_2025_batch_time_min": float(figure_s3[0][-1]),
                "predicted_from_2025_batch": conversion_2025,
                "equivalent_2026_batch_time_min": time_2026,
                "last_2026_batch_time_min": float(printed[0][-1]),
                "predicted_from_2026_batch": conversion_2026,
                "predicted_from_2026_acid_depletion": depletion,
            })
    return predictions, {"2025": loading_2025, "2026": loading_2026}


def compare(predictions):
    def condition(row):
        return int(row["temperature_C"]), float(row["flow_ml_min"])

    held_out = read_rows(HOLDOUT)
    measured = {condition(row): float(row["acid_conversion"]) for row in held_out}
    baseline = {condition(row): float(row["plug_flow_conversion"])
                for row in read_rows(ORIGINAL_COMPARISON)}
    conditions = set(map(condition, predictions))
    if len(measured) != len(held_out) or set(measured) != conditions or set(baseline) != conditions:
        raise ValueError("Batch transfer, holdout data and the original grid must match.")
    results = [
        dict(row, original_plug_flow_prediction=baseline[condition(row)],
             measured_acid_conversion=measured[condition(row)]) for row in predictions
    ]
    alternatives = ("predicted_from_2025_batch", "predicted_from_2026_batch",
                    "predicted_from_2026_acid_depletion")

    def score(subset, name):
        errors = np.array([row[name] - row["measured_acid_conversion"] for row in subset])
        return error_metrics(errors) | {"overpredicted_cases": int(np.count_nonzero(errors > 0))}

    scores_available = {
        name: score([row for row in results if row[name] is not None], name)
        for name in alternatives
    }
    common = [row for row in results if all(row[name] is not None for name in alternatives)]
    scores_common = {name: score(common, name) for name in (*alternatives, "original_plug_flow_prediction")}
    return results, scores_available, scores_common


def run(output):
    output.mkdir(parents=True, exist_ok=False)
    plan = read_json(PROTOCOL)
    predictions, loadings = make_predictions(plan)
    write_rows(output / "predictions_before_comparison.csv", predictions)
    results, all_available, common = compare(predictions)
    write_rows(output / "bed_comparison.csv", results)
    summary = {
        "model": "Interpolated measured batch progress at fixed catalyst exposure. No kinetic fit.",
        "inputs_sha256": {
            str(path.relative_to(ROOT)): sha(path)
            for path in (BATCH, FIGURE_S3, PROTOCOL, HOLDOUT, ORIGINAL_COMPARISON, HERE / "run.py")
        },
        "catalyst_kg_per_batch_mol_acid": loadings,
        "all_available_cases": all_available,
        "common_cases": common,
        "unavailable": [
            {"temperature_C": row["temperature_C"], "flow_ml_min": row["flow_ml_min"],
             "batch_series": year}
            for row in predictions for year in ("2025", "2026")
            if row[f"predicted_from_{year}_batch"] is None
        ],
        "bed_comparison_sha256": sha(output / "bed_comparison.csv"),
    }
    write_json(output / "summary.json", summary)
    return summary


def main():
    parser = ArgumentParser(description=__doc__)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    parser.add_argument("--output", type=Path, default=OUTPUT / stamp)
    args = parser.parse_args()
    result = run(args.output.resolve())
    print("Common-case comparison", result["common_cases"])
    print("Summary", args.output / "summary.json")


if __name__ == "__main__":
    main()
