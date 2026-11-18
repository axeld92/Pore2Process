"""Compare water inventory and bead-surface activity at matched conversion."""

from argparse import ArgumentParser
from datetime import datetime, timezone
from functools import partial
from math import exp, log, pi
from pathlib import Path

import numpy as np

from examples.biocatalysis.pentyl_propanoate.run import read_rows, write_rows
from examples.biocatalysis.pentyl_propanoate_film_water.film import run as prior_film
from examples.biocatalysis.pentyl_propanoate_film_water.water_activity import run as prior_batch
from examples.biocatalysis.pentyl_propanoate_transfer.published import run as published
from examples.biocatalysis.pentyl_propanoate_transport.phase import run as prior_phase
from pore2process_demo.bridge import runtime
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT

from .activity import (
    TernaryProjection,
    binary_activity,
    liquid_moles,
    measured_pentanol_curve,
    projected_activity,
)


HERE = Path(__file__).resolve().parent
INPUT = HERE / "inputs.json"
BATCH_COMPOSITION = ROOT / "literature/sources/lorenzo2026/batch_2to1.csv"
OUTPUT = ROOT / "results/biocatalysis/pentyl_propanoate_water_balance"
GAS_CONSTANT_J_MOL_K = 8.31446261815324
CONVERSIONS = tuple(step / 200 for step in range(151))
TRANSFER_CASES = (
    ("packed_bed", "wilson_geankoplis", 0),
    ("packed_bed_slow_water", "wilson_geankoplis", -1),
    ("quiescent_sphere", "isolated_sphere", 0),
)


def water_saturation_pressure(temperature_K, correlation):
    critical = correlation["critical_temperature_K"]
    if not 273.16 <= temperature_K <= critical:
        raise ValueError("Temperature is outside the reference water saturation curve.")
    reduced = 1 - temperature_K / critical
    polynomial = sum(factor * reduced**power for factor, power in
                     zip(correlation["factors"], correlation["exponents"], strict=True))
    return correlation["critical_pressure_Pa"] * exp(critical / temperature_K * polynomial)


def gas_water_at_saturation(temperature_K, gas_volume_ml, correlation):
    if not np.isfinite(gas_volume_ml) or gas_volume_ml < 0:
        raise ValueError("The gas volume must be finite and nonnegative.")
    return gas_volume_ml * 1e-6 * water_saturation_pressure(temperature_K, correlation) / (
        GAS_CONSTANT_J_MOL_K * temperature_K
    )


def water_diffusivity(liquid_properties, measurement, temperature_K, viscosity_Pa_s, uncertainty_units):
    reference_temperature = measurement["temperature_C"] + 273.15
    viscosities = liquid_properties["viscosity_density"]["measurements"]
    reference_viscosity = exp(np.interp(
        reference_temperature, [item["temperature_K"] for item in viscosities],
        [log(item["viscosity_Pa_s"]) for item in viscosities],
    ))
    reference_diffusivity = (measurement["diffusivity_m2_s"]
                             + uncertainty_units * measurement["reported_plus_minus_m2_s"])
    return reference_diffusivity * (temperature_K / reference_temperature) * (
        reference_viscosity / viscosity_Pa_s
    )


def activity_models(temperature_C, binary_curve, lower_spread_curve, ternary):
    models = {
        "hot_binary_pentanol": partial(binary_activity, curve=binary_curve),
        "hot_binary_pentanol_lower_spread": partial(binary_activity, curve=lower_spread_curve),
    }
    if 30 <= temperature_C <= 50:
        weight_at_50 = (temperature_C - 30) / 20
        projections = [(weight, ternary[temp]) for temp, weight in
                       ((30, 1 - weight_at_50), (50, weight_at_50)) if weight > 0]
        for ester_label, include_ester in (("ester_excluded", False), ("ester_as_pentanol", True)):
            for curve_label, curve in (("linear", None), ("binary_shape", binary_curve)):
                models[f"{ester_label}_{curve_label}"] = partial(
                    projected_activity, weighted_projections=projections,
                    ester_as_pentanol=include_ester, curve=curve,
                )
    return models


def difference(first, second):
    return None if first is None or second is None else first - second


def conditional_batch_rate_ratio(wetter, reference, series):
    if wetter is None or reference is None or not (
        series[0].water_activity <= min(wetter, reference)
        and max(wetter, reference) <= series[-1].water_activity
    ):
        return None
    activities = [point.water_activity for point in series]
    log_rates = [log(point.rate) for point in series]
    reference_log_rate, surface_log_rate = np.interp([reference, wetter], activities, log_rates)
    return exp(surface_log_rate - reference_log_rate)


def batch_composition_audit(path):
    """Check the selected published mass-balance water values against ester."""
    pairs = [(float(row["x_water"]), float(row["x_ester"])) for row in read_rows(path)
             if int(row["assay_id"]) in (3, 9, 10, 11, 12) and float(row["acid_conversion"]) > 0]
    return {
        "post_start_records": len(pairs),
        "equal_water_and_ester_at_printed_precision": sum(water == ester for water, ester in pairs),
        "largest_absolute_water_ester_difference": max(abs(water - ester) for water, ester in pairs),
        "source_method": "Water estimated by the authors' mass balance. These values provide no independent water activity or water retention measurement.",
    }


def calculate(inputs, phase_data, liquid_properties, reactor, jl, conversions=CONVERSIONS):
    molar_masses = phase_data["pure_liquids"]["molar_mass_g_mol"]
    acid_mass_per_mol, water_mass_per_mol = molar_masses[0], molar_masses[2]
    batch, moisture = inputs["batch_charge"], inputs["fresh_catalyst_moisture"]
    acid_mass = batch["initial_acid_volume_ml"] * batch["acid_density_at_25_C_g_ml"]
    batch_acid_mol = acid_mass / acid_mass_per_mol
    batch_catalyst_g = acid_mass * batch["catalyst_mass_fraction_of_acid"]
    gas_volume_ml = batch["vessel_volume_ml"] - batch["initial_liquid_charge_ml"]
    high_moisture_fraction = moisture["water_mass_fraction_as_received"] + moisture["reported_plus_minus_mass_fraction"]
    batch_water_mol = batch_catalyst_g * high_moisture_fraction / water_mass_per_mol

    bed = reactor["bed"]
    area = pi * bed["diameter_m"]**2 / 4
    liquid_volume_m3 = area * bed["length_m"] * bed["bed_porosity"]
    catalyst_kg, diameter, porosity = (
        bed["catalyst_mass_kg"], bed["bead_mean_diameter_m"], bed["bed_porosity"]
    )
    if catalyst_kg != 0.006:
        raise ValueError("Use the printed 6 g inventory.")
    catalyst_g_L = catalyst_kg / liquid_volume_m3
    exterior_area_per_liquid_m3 = 6 * (1 - porosity) / (diameter * porosity)
    acid, alcohol = bed["feed_acid_mol_per_L"], bed["feed_alcohol_mol_per_L"]
    if alcohol != 2 * acid:
        raise ValueError("Use the frozen nominal 2:1 alcohol-to-acid ratio.")

    model = jl.LorenzoPackedBed
    kinetic_params = model.PingPong(np.array([reactor["rate_parameters"][key] for key in published.PARAMETERS]))
    vle = inputs["pentanol_water_vle"]
    binary_curve = measured_pentanol_curve(vle["observations"])
    lower_spread_curve = measured_pentanol_curve([
        row for row in vle["observations"] if row[-1] <= vle["maximum_repeat_difference_for_sensitivity"]
    ])
    ternary = {
        round(float(temp) - 273.15): TernaryProjection.from_tie_lines(rows)
        for temp, rows in phase_data["ternary_liquid_liquid"]["mole_fraction_endpoints_by_temperature_K"].items()
    }
    phase_source_endpoints = []
    for temperature, projection in ternary.items():
        organic_water = projection.endpoints[:, 0, 2]
        aqueous_water = projection.endpoints[:, 1, 2]
        phase_source_endpoints.append({
            "temperature_C": temperature,
            "organic_water_mole_fraction": [float(organic_water.min()), float(organic_water.max())],
            "aqueous_water_mole_fraction": [float(aqueous_water.min()), float(aqueous_water.max())],
        })
    analogue_points = prior_batch.read_batch_series(prior_batch.BATCH_RATES)[
        ("Tjornelund2023", "MTBE", "nominal")
    ]
    low, high = analogue_points[0], analogue_points[-1]
    analogue_log_slopes = [prior_batch.endpoint_slope(a, b)
                          for a, b in zip(analogue_points[:-1], analogue_points[1:], strict=True)]
    lookup = {round(row["temperature_K"] - 273.15): row for row in liquid_properties["conditions"]}
    headspace_rows, results = [], []
    for temperature in bed["reaction_temperature_C"]:
        properties = lookup[temperature]
        temperature_K = properties["temperature_K"]
        rho, mu = properties["density_kg_m3"], properties["viscosity_Pa_s"]
        saturation = inputs["pure_water_vapor_pressure"]
        gas_water = gas_water_at_saturation(temperature_K, gas_volume_ml, saturation)
        maximum_loss = gas_water / batch_acid_mol
        headspace_rows.append({
            "temperature_C": temperature,
            "pure_water_vapor_pressure_Pa": water_saturation_pressure(temperature_K, saturation),
            "saturated_vapor_water_mol": gas_water,
            "saturated_vapor_water_mg": 1000 * gas_water * water_mass_per_mol,
            "water_per_initial_acid_mol_mol": maximum_loss,
        })
        eta = reactor["effectiveness_by_temperature"][str(temperature)]["eta_internal"]
        rates = {
            conversion: max(0.0, eta * float(model.disappearance_rate(
                kinetic_params, temperature, conversion, catalyst_g_L, acid0=acid, alcohol0=alcohol,
            )))
            for conversion in conversions
        }
        if np.any(np.diff(list(rates.values())) > 1e-12):
            raise ValueError("The film-free rate envelope requires a decreasing reference rate along the progress path.")
        activity = activity_models(temperature, binary_curve, lower_spread_curve, ternary)
        water_mobility = {name: water_diffusivity(
            liquid_properties, inputs["water_in_pentanol_diffusivity"], temperature_K, mu, deviation,
        ) for name, _, deviation in TRANSFER_CASES}
        for flow in bed["feed_flow_ml_per_min"]:
            velocity = (flow / 60e6) / area
            coefficients = {
                name: float(jl.BeadFilm.film_coefficients(
                    velocity, diameter, rho, mu, water_mobility[name], porosity, correlation,
                ).coefficient_m_s)
                for name, correlation, _ in TRANSFER_CASES
            }
            for conversion, rate in rates.items():
                bulk = liquid_moles(conversion)
                headspace = liquid_moles(conversion, -min(conversion, maximum_loss))
                released_water = liquid_moles(conversion, batch_water_mol / batch_acid_mol)
                # mol/(L min) to mol/(m³ s), then divide by a_l to get flux.
                water_flux = rate * 1000 / (60 * exterior_area_per_liquid_m3)
                for map_name, evaluate in activity.items():
                    bulk_activity = evaluate(bulk)
                    lower_batch_activity = evaluate(headspace)
                    moisture_batch_activity = evaluate(released_water)
                    for film_name, coefficient in coefficients.items():
                        delta_water_mol_L = water_flux / coefficient / 1000
                        surface = liquid_moles(conversion, delta_water_mol_L / acid)
                        surface_activity = evaluate(surface)
                        delta_total = difference(surface_activity, lower_batch_activity)
                        results.append({
                            "temperature_C": temperature,
                            "bed_flow_ml_min": flow,
                            "matched_acid_conversion": conversion,
                            "activity_scenario": map_name,
                            "water_film_scenario": film_name,
                            "bulk_water_mole_fraction": float(bulk[2] / bulk.sum()),
                            "water_diffusivity_m2_s": water_mobility[film_name],
                            "water_film_coefficient_m_s": coefficient,
                            "generated_water_flux_mol_m2_s": water_flux,
                            "surface_water_increment_mol_L": delta_water_mol_L,
                            "bulk_water_activity": bulk_activity,
                            "batch_aw_with_saturated_headspace_loss": lower_batch_activity,
                            "batch_aw_with_fresh_bead_water_release": moisture_batch_activity,
                            "bed_surface_water_activity": surface_activity,
                            "headspace_delta_aw": difference(bulk_activity, lower_batch_activity),
                            "bead_moisture_delta_aw": difference(moisture_batch_activity, bulk_activity),
                            "water_film_delta_aw": difference(surface_activity, bulk_activity),
                            "film_plus_headspace_delta_aw": delta_total,
                            "conditional_rate_relative_to_batch": conditional_batch_rate_ratio(
                                surface_activity, lower_batch_activity, analogue_points,
                            ),
                        })
    return headspace_rows, results, {
        "batch_acid_initial_mol": batch_acid_mol,
        "batch_catalyst_g": batch_catalyst_g,
        "bed_catalyst_g": catalyst_kg * 1000,
        "batch_nominal_headspace_ml": gas_volume_ml,
        "bed_liquid_volume_ml": liquid_volume_m3 * 1e6,
        "batch_high_analogue_moisture_mg": 1000 * batch_water_mol * water_mass_per_mol,
        "bed_high_analogue_moisture_mg": 1e6 * catalyst_kg * high_moisture_fraction,
        "batch_high_analogue_moisture_mol_per_initial_acid_mol": batch_water_mol / batch_acid_mol,
        "bed_high_analogue_moisture_mol_per_initial_acid_mol_in_one_void_volume":
            (1000 * catalyst_kg * high_moisture_fraction / water_mass_per_mol) / (1000 * liquid_volume_m3 * acid),
        "water_rate_analogue_segment_log_slopes": analogue_log_slopes,
        "water_rate_analogue_interval": [low.water_activity, high.water_activity],
        "binary_activity_data": list(zip(*binary_curve, strict=True)),
        "binary_activity_data_lower_repeat_spread": list(zip(*lower_spread_curve, strict=True)),
        "binary_reported_aw_max_difference_from_ideal_vapor_estimate": max(
            abs((1 - x_pentanol) * gamma - (1 - y_pentanol) * 1000 * vle["pressure_kPa"]
                / water_saturation_pressure(temperature_K, inputs["pure_water_vapor_pressure"]))
            for temperature_K, x_pentanol, y_pentanol, gamma, _ in vle["observations"]
        ),
        "phase_source_endpoints": phase_source_endpoints,
    }


def summarize(rows):
    summary = {}
    for name, _, _ in TRANSFER_CASES:
        cases = [row for row in rows if row["water_film_scenario"] == name]
        summary[name] = {}
        for label, subset in (
            ("all_conversions", cases),
            ("conversion_at_least_0_20", [row for row in cases if row["matched_acid_conversion"] >= 0.20]),
        ):
            ranges = {}
            for key in ("headspace_delta_aw", "bead_moisture_delta_aw", "water_film_delta_aw",
                        "film_plus_headspace_delta_aw", "conditional_rate_relative_to_batch"):
                values = [row[key] for row in subset if row[key] is not None]
                ranges[key] = [min(values), max(values)]
            summary[name][label] = ranges
    return summary


def run(output):
    inputs, phase_data, liquid, source = map(read_json, (
        INPUT, prior_phase.PROPERTIES, prior_film.PROPERTIES, published.MODEL_INPUTS,
    ))
    jl = runtime()
    for name, path in (("LorenzoPackedBed", published.MODEL), ("BeadFilm", prior_film.HERE / "model.jl")):
        if not jl.seval(f"isdefined(@__MODULE__, :{name})"):
            jl.include(str(path))
    headspace, comparisons, summary = calculate(inputs, phase_data, liquid, source, jl)
    output.mkdir(parents=True, exist_ok=False)
    write_rows(output / "headspace.csv", headspace)
    write_rows(output / "matched_activity.csv", comparisons)
    representative = [row for row in comparisons if row["water_film_scenario"] == "packed_bed"
                      and row["matched_acid_conversion"] in (0, 0.2, 0.5, 0.7)
                      and (row["temperature_C"], row["bed_flow_ml_min"]) in ((50, 3), (60, 3))]
    write_rows(output / "representative.csv", representative)
    summary.update({
        "kinetic_scope": "Water production uses the frozen batch rate law and internal effectiveness. A film-free bead flux screens the water concentration rise; no new conversion is fitted or solved.",
        "observed_reactor_conversions_read": 0,
        "published_batch_composition_check": batch_composition_audit(BATCH_COMPOSITION),
        "conditional_ranges": summarize(comparisons),
        "julia": str(jl.seval("VERSION")),
        "source_sha256": {
            str(path.relative_to(ROOT)): sha(path)
            for path in (INPUT, BATCH_COMPOSITION, prior_phase.PROPERTIES, prior_phase.HERE / "run.py",
                         prior_film.PROPERTIES, prior_film.HERE / "model.jl", published.MODEL_INPUTS,
                         published.MODEL, prior_batch.BATCH_RATES, prior_batch.HERE / "run.py",
                         HERE / "run.py", HERE / "activity.py")
        },
        "output_sha256": {name: sha(output / name) for name in (
            "headspace.csv", "matched_activity.csv", "representative.csv",
        )},
    })
    write_json(output / "summary.json", summary)
    return summary


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    parser.add_argument("--output", type=Path, default=OUTPUT / stamp)
    args = parser.parse_args()
    report = run(args.output.resolve())
    for name, result in report["conditional_ranges"].items():
        maximum = result["conversion_at_least_0_20"]["water_film_delta_aw"][1]
        print(name, f"maximum activity rise above 20% conversion {maximum:.4f}")
