"""Test initial-rate responses to water activity using independent batch data."""

from argparse import ArgumentParser
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from math import exp, isfinite, log
from pathlib import Path

from examples.biocatalysis.pentyl_propanoate.run import read_rows, write_rows
from pore2process_demo.io import sha, write_json
from pore2process_demo.paths import ROOT


HERE = Path(__file__).resolve().parent
BATCH_RATES = HERE / "batch_initial_rates.csv"
OUTPUT = ROOT / "results/biocatalysis/pentyl_propanoate_film_water/water_activity"


@dataclass(frozen=True)
class BatchRate:
    water_activity: float
    rate: float


def read_batch_series(path):
    series = defaultdict(list)
    for row in read_rows(path):
        source, solvent = row["source_id"], row["solvent"]
        rate = float(row["specific_activity_as_printed"])
        coordinates = {"nominal": row["initial_water_activity"]}
        if row["figure_initial_water_activity_approx"]:
            coordinates["figure_approximate"] = row["figure_initial_water_activity_approx"]
        for basis, coordinate in coordinates.items():
            key = source, solvent, basis
            water_activity = float(coordinate)
            if not (0 <= water_activity <= 1 and isfinite(rate) and rate > 0):
                raise ValueError(f"Invalid water activity or rate in {key}.")
            series[key].append(BatchRate(water_activity, rate))
    for key, rates in series.items():
        rates.sort(key=lambda point: point.water_activity)
        if len(rates) < 3 or len({point.water_activity for point in rates}) != len(rates):
            raise ValueError(f"Need distinct water activities and an interior check in {key}.")
    return series


def endpoint_slope(low, high):
    return log(high.rate / low.rate) / (high.water_activity - low.water_activity)


def compare(scenario, key, slope, anchor, test_points):
    source, solvent, activity_basis = key
    rows = []
    for point in test_points:
        measured = point.rate / anchor.rate
        predicted = exp(slope * (point.water_activity - anchor.water_activity))
        rows.append({
            "comparison": scenario,
            "source_id": source,
            "solvent": solvent,
            "activity_basis": activity_basis,
            "anchor_water_activity": anchor.water_activity,
            "test_water_activity": point.water_activity,
            "observed_relative_rate": measured,
            "predicted_relative_rate": predicted,
            "signed_error_in_relative_rate": predicted - measured,
        })
    return rows


def analyse(series):
    comparisons, within_study, transfer = [], [], []
    for key, rates in series.items():
        low, high = rates[0], rates[-1]
        slope = endpoint_slope(low, high)
        scored = compare("same_series_endpoints", key, slope, low, rates[1:-1])
        within_study.append({
            "series": key,
            "low_water_activity": low.water_activity,
            "high_water_activity": high.water_activity,
            "high_rate_relative_to_low": high.rate / low.rate,
            "slope_per_unit_water_activity": slope,
            "scored_count": len(scored),
            "max_absolute_error_in_relative_rate": max(abs(row["signed_error_in_relative_rate"]) for row in scored),
        })
        comparisons.extend(scored)

    donor_key = ("Tjornelund2023", "MTBE", "nominal")
    low, high = series[donor_key][0], series[donor_key][-1]
    slope = endpoint_slope(low, high)
    for basis in ("nominal", "figure_approximate"):
        key = ("DuttaBanik2017", "MTBE", basis)
        receiver = series[key]
        # Anchor both readings of the second paper at the same measured rate.
        anchor = next(point for point in receiver if point.water_activity == 0.23)
        test_points = [
            point for point in receiver
            if anchor.water_activity < point.water_activity <= high.water_activity
        ]
        if not low.water_activity <= anchor.water_activity or not test_points:
            raise ValueError("The external check must stay inside the donor's activity range.")
        scored = compare("separate_MTBE_study", key, slope, anchor, test_points)
        transfer.append({
            "donor_series": donor_key,
            "recipient_series": key,
            "recipient_anchor_water_activity": anchor.water_activity,
            "slope_per_unit_water_activity": slope,
            "scored_count": len(scored),
            "max_absolute_error_in_relative_rate": max(abs(row["signed_error_in_relative_rate"]) for row in scored),
            "highest_activity_comparison": scored[-1],
        })
        comparisons.extend(scored)
    return comparisons, {"within_study": within_study, "between_studies": transfer}


def run(output):
    comparisons, summary = analyse(read_batch_series(BATCH_RATES))
    output.mkdir(parents=True, exist_ok=False)
    write_rows(output / "comparison.csv", comparisons)
    summary.update({
        "fit_scope": "Published batch initial rates only. Endpoints supply the slope. No reactor data enter the analysis.",
        "target_reaction_water_activity": None,
        "source_sha256": {
            str(path.relative_to(ROOT)): sha(path) for path in (BATCH_RATES, HERE / "run.py")
        },
        "comparison_sha256": sha(output / "comparison.csv"),
    })
    write_json(output / "summary.json", summary)
    return summary


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    parser.add_argument("--output", type=Path, default=OUTPUT / stamp)
    args = parser.parse_args()
    report = run(args.output.resolve())
    for series in report["within_study"]:
        print(*series["series"], f'high/low rate {series["high_rate_relative_to_low"]:.3f}')
    for transfer in report["between_studies"]:
        at_high = transfer["highest_activity_comparison"]
        print(*transfer["recipient_series"], f'a_w {at_high["test_water_activity"]:.2f}',
              f'predicted {at_high["predicted_relative_rate"]:.3f}, measured {at_high["observed_relative_rate"]:.3f}')
