"""Audit the conversion readout using independently published liquid data."""

import argparse
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from examples.biocatalysis.pentyl_propanoate.run import error_metrics, read_rows, write_rows
from examples.biocatalysis.pentyl_propanoate_transfer.published import run as published
from pore2process_demo.bridge import runtime
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT


HERE = Path(__file__).resolve().parent
PROPERTIES = HERE / "properties.json"
OUTPUT = ROOT / "results/biocatalysis/pentyl_propanoate_transport/phase"


def acid_conversion_readout(remaining_acid_mol, sample_volume_ml, feed_volume_ml):
    """Read 1 - C_sample/C_feed on a basis of one mole of acid entering."""
    if sample_volume_ml <= 0 or feed_volume_ml <= 0 or remaining_acid_mol < 0:
        raise ValueError("The sample and feed require positive volumes and nonnegative acid.")
    return 1 - remaining_acid_mol * feed_volume_ml / sample_volume_ml


def pair_excess_volume(acid, cosolvent, coefficients):
    total = acid + cosolvent
    if total == 0:
        return 0.0
    fraction = acid / total
    return total * fraction * (1 - fraction) * np.polynomial.polynomial.polyval(
        2 * fraction - 1, coefficients
    )


def ternary_split(moles, endpoints):
    """Find a measured or interpolated tie line crossing the ternary composition."""
    total = float(sum(moles))
    composition = moles / total
    endpoints = np.array(endpoints, dtype=float).reshape(-1, 2, 3)
    endpoints /= endpoints.sum(axis=2, keepdims=True)
    for first, last in zip(endpoints[:-1], endpoints[1:], strict=True):
        def residual(weight):
            organic, aqueous = first + weight * (last - first)
            line = aqueous - organic
            offset = composition - organic
            return line[0] * offset[1] - line[1] * offset[0]

        left, right = 0.0, 1.0
        sign_left = residual(left)
        if sign_left * residual(right) > 0:
            continue
        for _ in range(56):
            center = (left + right) / 2
            if sign_left * residual(center) > 0:
                left = center
            else:
                right = center
        organic, aqueous = first + center * (last - first)
        fraction = (composition[2] - organic[2]) / (aqueous[2] - organic[2])
        if 1e-12 < fraction < 1 - 1e-12:
            amounts = total * (1 - fraction) * organic, total * fraction * aqueous
            if not np.allclose(amounts[0] + amounts[1], moles, rtol=0, atol=1e-12):
                raise RuntimeError("The interpolated tie line does not close the ternary balance.")
            return amounts
    return None


def scenario_names(properties):
    return (
        "fixed_6g_molar_conversion", "whole_liquid_additive", "whole_liquid_binary_excess",
        "omit_generated_water", "omit_water_printed_inlet_basis",
        *(f"{phase}_{temperature.replace('.', '_')}_K" for temperature in
          properties["ternary_liquid_liquid"]["mole_fraction_endpoints_by_temperature_K"]
          for phase in ("pentanol_rich", "water_rich")),
    )


def readout_scenarios(conversion, source, properties):
    component_data = properties["pure_liquids"]
    molar_volume_ml = 1000 * np.array(component_data["molar_mass_g_mol"]) / np.array(
        component_data["density_kg_m3"]
    )
    acid, alcohol, water, ester = moles = np.array([1 - conversion, 2 - conversion,
                                                    conversion, conversion])
    feed_add = float(np.dot([1, 2, 0, 0], molar_volume_ml))
    feed_printed = 1000 / source["bed"]["feed_acid_mol_per_L"]
    outlet_add = float(moles @ molar_volume_ml)
    water_omitted = outlet_add - water * molar_volume_ml[2]
    binary = properties["binary_excess_volume"]
    acid_pentanol = binary["acid_pentanol"]["coefficients_cm3_mol"]
    acid_water = binary["acid_water"]["coefficients_cm3_mol"]
    feed_binary = feed_add + pair_excess_volume(1, 2, acid_pentanol)
    outlet_binary = (outlet_add + pair_excess_volume(acid, alcohol, acid_pentanol)
                     + pair_excess_volume(acid, water, acid_water))
    scenarios = {
        "fixed_6g_molar_conversion": conversion,
        "whole_liquid_additive": acid_conversion_readout(acid, outlet_add, feed_add),
        "whole_liquid_binary_excess": acid_conversion_readout(acid, outlet_binary, feed_binary),
        "omit_generated_water": acid_conversion_readout(acid, water_omitted, feed_add),
        "omit_water_printed_inlet_basis": acid_conversion_readout(acid, water_omitted, feed_printed),
        "whole_liquid_additive_volume_ratio": outlet_add / feed_add,
        "whole_liquid_binary_excess_volume_ratio": float(outlet_binary / feed_binary),
        "omit_generated_water_volume_ratio": float(water_omitted / feed_add),
    }
    for temperature, endpoints in properties["ternary_liquid_liquid"]["mole_fraction_endpoints_by_temperature_K"].items():
        label = temperature.replace(".", "_") + "_K"
        split = ternary_split(moles[:3], endpoints)
        if split is None:
            scenarios.update({f"pentanol_rich_{label}": None, f"water_rich_{label}": None})
            continue
        organic_moles, aqueous_moles = split
        # Ester is parked in the pentanol-rich phase without affecting the ternary tie lines.
        organic_volume = float(organic_moles @ molar_volume_ml[:3] + ester * molar_volume_ml[3])
        aqueous_volume = float(aqueous_moles @ molar_volume_ml[:3])
        scenarios[f"pentanol_rich_{label}"] = acid_conversion_readout(float(organic_moles[0]), organic_volume, feed_add)
        scenarios[f"water_rich_{label}"] = acid_conversion_readout(float(aqueous_moles[0]), aqueous_volume, feed_add)
    return scenarios


def predict(source, properties):
    jl = runtime()
    if not jl.seval("isdefined(@__MODULE__, :LorenzoPackedBed)"):
        jl.include(str(published.MODEL))
    model = jl.LorenzoPackedBed
    params = model.PingPong(np.array([source["rate_parameters"][name] for name in published.PARAMETERS]))
    base_rows, _, numerics = published.predict(model, params, source)
    if numerics["nominal_catalyst_mass_g"] != 6.0:
        raise ValueError("This audit requires the previously used, full 6 g catalyst loading.")
    predictions = []
    for row in base_rows:
        predictions.append(dict(
            temperature_C=row["temperature_C"], flow_ml_min=row["flow_ml_min"],
            **readout_scenarios(row["void_liquid_eq4"], source, properties),
        ))
    return predictions, numerics


def compare(predictions, scenarios, observations):
    observed = {published.condition(row): float(row["acid_conversion"]) for row in observations}
    if (len(observed) != 20 or len(observations) != 20 or len(predictions) != 20
            or set(observed) != set(map(published.condition, predictions))):
        raise ValueError("Expected the 20 held-out reactor measurements on the frozen grid.")
    if any(not 0 <= value < 1 for value in observed.values()):
        raise ValueError("The volume diagnostic requires observed conversion in [0, 1).")
    comparisons = []
    for row in predictions:
        measurement = observed[published.condition(row)]
        baseline = row["fixed_6g_molar_conversion"]
        comparisons.append(dict(
            row,
            observed_acid_conversion=measurement,
            volume_ratio_needed_if_all_unreacted_acid_is_sampled=(1 - baseline) / (1 - measurement),
        ))

    results = {}
    for name in scenarios:
        available = [row for row in comparisons if row[name] is not None]
        if not available:
            continue
        measured = np.array([row["observed_acid_conversion"] for row in available])
        baseline = np.array([row["fixed_6g_molar_conversion"] for row in available])
        predicted = np.array([row[name] for row in available])
        error = predicted - measured
        shift_pp = 100 * (predicted - baseline)
        results[name] = dict(
            **error_metrics(error),
            overpredicted_cases=int(np.count_nonzero(error > 0)),
            baseline_rmse_on_same_cases_pp=error_metrics(baseline - measured)["rmse_percentage_points"],
            min_shift_from_molar_conversion_pp=float(min(shift_pp)),
            max_shift_from_molar_conversion_pp=float(max(shift_pp)),
        )
    return comparisons, results


def run(output):
    source = read_json(published.MODEL_INPUTS)
    properties = read_json(PROPERTIES)
    predictions, numerics = predict(source, properties)
    output.mkdir(parents=True, exist_ok=False)
    write_rows(output / "predictions_before_comparison.csv", predictions)
    # All readout projections are fixed before opening the measured conversion file.
    comparisons, results = compare(predictions, scenario_names(properties), read_rows(published.HOLDOUT))
    write_rows(output / "comparison.csv", comparisons)
    anchors = [row for row in comparisons if (row["temperature_C"], row["flow_ml_min"]) == (60, 3.0)]
    summary = {
        "nominal_catalyst_mass_g": numerics["nominal_catalyst_mass_g"],
        "fixed_6g_rate_law": "Original Eq. S1, mass-consistent 6 g per void liquid basis, printed Table B1 internal effectiveness",
        "forecast_type": "Readout postprocessing under specified liquid scenarios",
        "pure_volumes_and_binary_excess_reference_K": properties["pure_liquids"]["reference_temperature_K"],
        "max_volume_ratio_change_by_scenario": {
            key: max(abs(row[key] - 1) for row in comparisons)
            for key in ("whole_liquid_additive_volume_ratio", "whole_liquid_binary_excess_volume_ratio",
                        "omit_generated_water_volume_ratio")
        },
        "comparison": results,
        "sixty_C_three_ml_min": anchors[0],
        "required_volume_ratio_if_no_acid_removed_min_max": [
            reduce(row["volume_ratio_needed_if_all_unreacted_acid_is_sampled"] for row in comparisons)
            for reduce in (min, max)
        ],
        "source_sha256": {
            str(path.relative_to(ROOT)): sha(path) for path in (
                HERE / "run.py", PROPERTIES, published.MODEL_INPUTS, published.MODEL,
                Path(published.__file__).resolve(), published.HOLDOUT,
            )
        },
        "comparison_sha256": sha(output / "comparison.csv"),
    }
    write_json(output / "summary.json", summary)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    default = OUTPUT / datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    parser.add_argument("--output", type=Path, default=default)
    args = parser.parse_args()
    summary = run(args.output.resolve())
    print(args.output.resolve() / "summary.json")
    for label, stats in summary["comparison"].items():
        print(f"{label}: {stats['rmse_percentage_points']:.3f} pp RMSE, "
              f"{stats['overpredicted_cases']}/{stats['count']} high")


if __name__ == "__main__":
    main()
