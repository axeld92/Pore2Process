"""Audit catalyst conditioning using independent sorption, tracer and reuse data."""

from argparse import ArgumentParser
from datetime import datetime, timezone
from itertools import product
from math import pi
from pathlib import Path

import numpy as np

from examples.biocatalysis.pentyl_propanoate.run import array, error_metrics, read_rows, write_rows
from examples.biocatalysis.pentyl_propanoate_film_water.film import run as film
from examples.biocatalysis.pentyl_propanoate_film_water.water_activity import run as batch
from examples.biocatalysis.pentyl_propanoate_transfer.published import run as published
from examples.biocatalysis.pentyl_propanoate_transfer.mixing import run as mixing
from examples.biocatalysis.pentyl_propanoate_transport.phase import run as phase
from examples.biocatalysis.pentyl_propanoate_water_balance import run as water
from examples.biocatalysis.pentyl_propanoate_water_balance.activity import (
    TernaryProjection, liquid_moles, measured_pentanol_curve,
)
from pore2process_demo.bridge import runtime
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT

HERE = Path(__file__).resolve().parent
INPUT = HERE / "inputs.json"
PARTICLES = HERE.parent / "pentyl_propanoate_transport/particle/evidence.json"
OUTPUT = ROOT / "results/biocatalysis/pentyl_propanoate_conditioning"
FILM_REFERENCE = film.HERE / "reference/bed_comparison.csv"


def bead_storage(total_kg_m3, catalyst_kg_m3, qmax, K):
    """Finite-batch equilibrium for initially dry catalyst and retained water."""
    if not (all(np.isfinite([total_kg_m3, catalyst_kg_m3, qmax, K]))
            and min(total_kg_m3, catalyst_kg_m3, qmax) >= 0 and K > 0):
        raise ValueError("Use nonnegative water and catalyst inventories and a positive K.")
    lo, hi = 0., total_kg_m3
    for _ in range(60):
        c = (lo + hi) / 2
        if c + catalyst_kg_m3 * qmax * c / (K + c) > total_kg_m3:
            hi = c
        else:
            lo = c
    concentration = (lo + hi) / 2
    return concentration, qmax * concentration / (K + concentration)


def normalized_moments(rows):
    t, c = (array(rows, name) for name in ("time_min", "concentration_mg_ml"))
    if not (len(t) > 1 and np.all(np.isfinite(t)) and np.all(np.diff(t) > 0)
            and np.all(np.isfinite(c) & (c >= 0))):
        raise ValueError("Use increasing measurement times and finite nonnegative concentrations.")
    area = float(np.trapezoid(c, x=t))
    if area <= 0:
        raise ValueError("Use a positive tracer integral.")
    mean = float(np.trapezoid(t * c, x=t) / area)
    return {"area_mg_min_ml": area, "mean_min": mean,
            "variance_min2": float(np.trapezoid((t - mean)**2 * c, x=t) / area)}


def exchange_time(blue, azide, inlet_mean, inlet_variance):
    """Conditional exponential mobile/immobile-exchange clock in minutes."""
    mu = blue["mean_min"] - inlet_mean
    var = blue["variance_min2"] - inlet_variance
    extra = azide["mean_min"] - blue["mean_min"]
    if not (mu > 0 and var >= 0 and extra > 0):
        raise ValueError("The inferred mobile residence and extra storage must be positive.")
    beta = extra / mu
    tau = (azide["variance_min2"] - inlet_variance - (1 + beta)**2 * var) / (2 * extra)
    return {"immobile_to_mobile_storage_ratio": beta, "return_time_min": tau,
            "physical_nonnegative_return_time": bool(tau >= 0)}


def tracer_check(source, rows):
    moments = {name: normalized_moments([row for row in rows if row["tracer"] == name])
               for name in ("blue_dextran", "sodium_azide")}
    flow = source["water_flow_ml_min"]
    duration = source["pulse_duration_min_if_injected_at_water_flow"]
    for name, result in moments.items():
        dose = source["pulse_volume_ml"] * source["reported_pulse_concentration_mg_ml"][name]
        result["reported_injected_dose_mg"] = dose
        result["integrated_output_dose_mg"] = flow * result["area_mg_min_ml"]
    records = {}
    for key, data in (("digitized", moments), ("printed", source["reported_rtd_moments"])):
        b, a = data["blue_dextran"], data["sodium_azide"]
        # Shared, but unmeasured pulse constrained to [0, pulse duration].
        # The mean is >= 0; its variance is <= duration^2/4. Using separate
        # extremes gives an upper envelope, not a fitted pulse.
        upper = exchange_time(b, a, 0., duration**2 / 4)
        records[key] = {
            "shared_rectangular_pulse": exchange_time(b, a, duration/2, duration**2/12),
            "conditional_return_time_upper_min": upper["return_time_min"],
            "mean_extra_storage_ml": flow * (a["mean_min"] - b["mean_min"]),
        }
    return {"normalized_moments": moments, "pulse_duration_scenario_min": duration,
            "exchange": records,
            "scope": "Shared pulse and mobile hydrodynamic clock, Poisson entry and exponential exit from immobile storage. The recovered tracer masses disagree with doses. Return time is a conditional aqueous estimate, not reactive-liquid wetting time."}


def reuse_proxies(inputs):
    digit = inputs["same_alcohol_reuse"]["time_proxy_digitization"]
    x0, x20 = digit["x_ticks_0_20_h_px"]
    y0, y100 = digit["y_ticks_0_100_pct_px"]
    y = y0 + (y100 - y0) * digit["target_conversion_pct"] / 100
    times, reading_bounds = {}, {}
    def crossing_time(xa, ya, xb, yb):
        if not xa < xb or not yb <= y <= ya:
            raise ValueError("The plotted markers must bracket the fixed conversion.")
        crossing_x = xa + (xb - xa) * (ya - y) / (ya - yb)
        return (crossing_x - x0) * 20 / (x20 - x0)

    pixels = digit["reading_half_range_px"]
    for cycle, pair in digit["bracketing_marker_xy_px"].items():
        (xa, ya), (xb, yb) = pair
        times[cycle] = crossing_time(xa, ya, xb, yb)
        perturbed = [crossing_time(xa+dx1, ya+dy1, xb+dx2, yb+dy2)
                     for dx1, dy1, dx2, dy2 in product((-pixels, pixels), repeat=4)]
        reading_bounds[cycle] = [min(perturbed), max(perturbed)]
    # Independent comparison at one fixed progress, not residual-activity fit.
    mild = min(times["1"] / times[cycle] for cycle in ("5", "10"))
    values = inputs["propanoic_acid_contact"]["activities_U_g"]
    assay = {key: val["mean"] / values["fresh"]["mean"] for key, val in values.items()}
    acid_vs_alcohol = assay["after_propanoic_acid"] / assay["after_benzyl_alcohol"]
    return {"acetate_time_to_70pct_h": times, "acetate_marker_reading_range_h": reading_bounds,
            "acetate_mild_time_scale": mild,
            "acetate_mild_time_scale_marker_reading_range": [
                reading_bounds["1"][0] / max(reading_bounds[c][1] for c in ("5", "10")),
                reading_bounds["1"][1] / max(reading_bounds[c][0] for c in ("5", "10")),
            ],
            "propanoic_assay_activity_relative_to_fresh": assay,
            "propanoic_acid_vs_benzyl_alcohol_activity_ratio": acid_vs_alcohol}


def trace_stability(inputs):
    digit = inputs["bed_time_trace"]["digitization"]
    y0, y100 = digit["y_ticks_0_100_pct_px"]
    values = (y0 - np.array(digit["marker_y_px"])) / (y0 - y100)
    steady = values[1:]
    return {
        "times_h": digit["times_h"], "conversion": values.tolist(),
        "observed_30_to_180min_span_pp": float(100*np.ptp(steady)),
        "observed_30_to_180min_change_pp": float(100*(steady[-1] - steady[0])),
        "reading_half_range_pp": 100*digit["reading_half_range_px"]/(y0-y100),
    }


def batch_water(inputs, water_inputs, phase_data, bed):
    charge = water_inputs["batch_charge"]
    acid_MW, _, water_MW, _ = phase_data["pure_liquids"]["molar_mass_g_mol"]
    acid_mass = charge["initial_acid_volume_ml"] * charge["acid_density_at_25_C_g_ml"]
    acid_mol = acid_mass / acid_MW
    cat_g = acid_mass * charge["catalyst_mass_fraction_of_acid"]
    volume_L = charge["initial_liquid_charge_ml"] / 1000
    moisture = water_inputs["fresh_catalyst_moisture"]
    initial_water_sensitivity = (moisture["water_mass_fraction_as_received"]
                                 + moisture["reported_plus_minus_mass_fraction"])
    vle = water_inputs["pentanol_water_vle"]
    curve = measured_pentanol_curve(vle["observations"])
    lower = measured_pentanol_curve([
        row for row in vle["observations"] if row[-1] <= vle["maximum_repeat_difference_for_sensitivity"]
    ])
    ternary = {
        round(float(t)-273.15): TernaryProjection.from_tie_lines(rows) for t, rows in
        phase_data["ternary_liquid_liquid"]["mole_fraction_endpoints_by_temperature_K"].items()
    }
    analogue = batch.read_batch_series(batch.BATCH_RATES)[("Tjornelund2023", "MTBE", "nominal")]
    batch_storage_rows, activity_rows = [], []
    isotherms = inputs["moisture_isotherms"]
    for name, iso in isotherms["cases"].items():
        qmax, K = iso["qmax_kg_kg"], iso["K_kg_m3"]
        for conversion in (0.01, 0.02, 0.04, 0.1, 0.2, 0.5, 0.65, 0.75):
            generated_g = acid_mol * conversion * water_MW
            max_retained_g = min(generated_g, cat_g * qmax)
            # 1 g/L = 1 kg/m³. Only report the fitted equilibrium within a
            # common concentration range of the two source graphs.
            c, q = bead_storage(generated_g / volume_L, cat_g / volume_L, qmax, K)
            in_range = c <= isotherms["common_test_concentration_max_kg_m3"]
            batch_storage_rows.append({
                "case": name, "conversion": conversion,
                "produced_water_mg": 1000*generated_g, "max_retained_water_mg": 1000*max_retained_g,
                "max_fraction_produced_retained": max_retained_g / generated_g,
                "max_net_retention_with_initial_moisture_mg": 1000*min(
                    generated_g, cat_g*max(0., qmax - initial_water_sensitivity)),
                "within_common_isotherm_concentration_range": in_range,
                "within_range_equilibrium_free_water_g_L": c if in_range else None,
                "within_range_equilibrium_retained_water_mg": q*cat_g*1000 if in_range else None,
            })
            moles = liquid_moles(conversion)
            dry_bead_capacity_moles = liquid_moles(conversion, -max_retained_g/(acid_mol*water_MW))
            for temperature in bed["reaction_temperature_C"]:
                for label, evaluate in water.activity_models(temperature, curve, lower, ternary).items():
                    bulk, retained = evaluate(moles), evaluate(dry_bead_capacity_moles)
                    activity_rows.append({
                        "case": name, "conversion": conversion, "temperature_C": temperature,
                        "activity_scenario": label, "bulk_aw_equal_retained_inventories": bulk,
                        "minimum_batch_aw_at_analogue_capacity": retained,
                        "bed_minus_batch_aw_scenario": water.difference(bulk, retained),
                        "conditional_bed_batch_rate_ratio": water.conditional_batch_rate_ratio(
                            bulk, retained, analogue),
                    })
    capacities = {name: {"batch_water_capacity_mg": 1000*cat_g*iso["qmax_kg_kg"],
                         "bed_water_capacity_mg": 1e6*bed["catalyst_mass_kg"]*iso["qmax_kg_kg"]}
                  for name, iso in isotherms["cases"].items()}
    activity_summary = {}
    for conversion in sorted({row["conversion"] for row in activity_rows}):
        rows = [row for row in activity_rows if row["conversion"] == conversion]
        shifts = [row["bed_minus_batch_aw_scenario"] for row in rows
                  if row["bed_minus_batch_aw_scenario"] is not None]
        ratios = [row["conditional_bed_batch_rate_ratio"] for row in rows
                  if row["conditional_bed_batch_rate_ratio"] is not None]
        activity_summary[str(conversion)] = {
            "max_bed_minus_batch_aw_scenario": max(shifts, default=None),
            "min_conditional_analogue_bed_batch_rate_ratio": min(ratios, default=None),
            "activity_scenarios_inside_composition_range": len(shifts),
            "rate_scenarios_inside_controlled_aw_range": len(ratios),
        }
    return batch_storage_rows, activity_rows, {
        "batch_acid_mol": acid_mol, "batch_catalyst_g": cat_g,
        "initial_water_mass_fraction_sensitivity": initial_water_sensitivity,
        "analogue_saturation_capacity": capacities,
        "activity_sensitivity_summary": activity_summary,
        "use": "Initially dry, nominal constant liquid volume. At finite conversion, compare a batch retaining up to the analogue saturation capacity with a steady bed having no net storage. This is a deliberately large conditional inventory change, not a bound for the actual unmeasured mixture.",
    }


def bed_mean_remaining(time_s, remaining, elapsed_min, displacement_min):
    """Step exposure ages in a plug-displaced bed, with pre-contact age <= 0."""
    age = elapsed_min*60 - np.linspace(0., displacement_min*60, 2049)
    return float(np.trapezoid(np.interp(age, time_s, remaining, left=1.),
                              dx=1/2048))


def moisture_screen(jl, inputs, water_inputs, liquid, reactor, particle_data):
    plan = inputs["transient_screen"]
    bed = reactor["bed"]
    area = pi * bed["diameter_m"]**2 / 4
    total_volume = area * bed["length_m"]
    void = bed["bed_porosity"]
    cat_density = bed["catalyst_mass_kg"] / ((1-void)*total_volume)
    radius = bed["bead_mean_diameter_m"]/2
    tau = reactor["dispersion_and_transport"]["particle_tortuosity_used"]
    T = plan["temperature_C"] + 273.15
    prop = next(row for row in liquid["conditions"] if row["temperature_K"] == T)
    d = water.water_diffusivity(liquid, water_inputs["water_in_pentanol_diffusivity"], T,
                                prop["viscosity_Pa_s"], plan["water_diffusivity_uncertainty_units"])
    k = float(jl.BeadFilm.film_coefficients(
        plan["film_flow_ml_min"]/60e6/area, 2*radius, prop["density_kg_m3"],
        prop["viscosity_Pa_s"], d, void, plan["film"],
    ).coefficient_m_s)
    dry_eps = particle_data["lorenzo_2025"]["particle_porosity"]["dry_nitrogen_sorption"]["value"]
    porosities = {"wet_aqueous_basis": bed["particle_porosity_wet"], "dry_nitrogen_basis": dry_eps}
    isotherms = inputs["moisture_isotherms"]["cases"]
    setups = [(name, eps_name, eps, iso["qmax_kg_kg"], iso["K_kg_m3"], c0, c1)
              for (name, iso), (eps_name, eps), concentration, direction in product(
                  isotherms.items(), porosities.items(),
                  inputs["moisture_isotherms"]["common_test_concentrations_kg_m3"], ("wetting", "drying"))
              for c0, c1 in [(0., concentration) if direction == "wetting" else (concentration, 0.)]]
    setups += [("no_sorption", "wet_aqueous_basis", porosities["wet_aqueous_basis"], 0., 1., 0., 1.)]
    curves, summary_rows = [], []
    max_grid, max_balance = 0., 0.
    fields = ("remaining", "remaining_shell", "pore_remaining")
    for name, eps_name, eps, qmax, K, c0, cb in setups:
        solutions = [jl.BeadMoisture.relax_water(
            radius, eps, tau, cat_density, d, k, qmax, K, c0, cb, cells=n, dt=dt,
            duration=plan["duration_s"], sample_step=plan["sample_step_s"],
            shell_depth=plan["outer_shell_depth_m"],
        ) for n, dt in zip(plan["grid_cells"], plan["step_s"], strict=True)]
        coarse, fine = solutions
        times = np.array(fine.times_s)
        series = {field: np.array(getattr(fine, field)) for field in fields}
        difference = max(float(np.max(abs(np.array(getattr(coarse, f)) - series[f]))) for f in fields)
        max_grid = max(max_grid, difference)
        max_balance = max(max_balance, *(float(sol.max_relative_balance_error) for sol in solutions))
        if any(not np.all((values >= -1e-8) & (values <= 1+1e-8))
               or np.any(np.diff(values) > 1e-8) for values in series.values()):
            raise RuntimeError("The passive water step must relax monotonically within its inventory range.")
        case = {"isotherm": name, "porosity_basis": eps_name, "particle_porosity": eps,
                "initial_c_kg_m3": c0, "boundary_c_kg_m3": cb}
        for i, t in enumerate(times):
            curves.append(case | {"time_min": t/60} | {f: values[i] for f, values in series.items()})
        for flow in bed["feed_flow_ml_per_min"]:
            displacement = total_volume*void*1e6/flow
            for elapsed in (30, 180):
                summary_rows.append(case | {
                    "flow_ml_min": flow, "elapsed_min": elapsed, "void_displacement_min": displacement,
                    **{f: float(np.interp(elapsed*60, times, series[f])) for f in fields},
                    **{f"plug_age_average_{f}": bed_mean_remaining(times, series[f], elapsed, displacement)
                       for f in fields},
                })
    if max_grid > plan["max_grid_difference_fraction"] or max_balance > plan["max_mass_balance_relative_error"]:
        raise RuntimeError(f"Unresolved moisture-grid or balance check: {max_grid=}, {max_balance=}")
    slow_flow = min(bed["feed_flow_ml_per_min"])
    memory = {}
    for porosity, direction, elapsed in product(porosities, ("wetting", "drying"), (30, 180)):
        selected = [r for r in summary_rows if r["porosity_basis"] == porosity
                    and r["isotherm"] != "no_sorption" and r["flow_ml_min"] == slow_flow
                    and r["elapsed_min"] == elapsed and
                    (r["initial_c_kg_m3"] < r["boundary_c_kg_m3"]) == (direction == "wetting")]
        memory[f"{porosity}_{direction}_{elapsed}min"] = {
            f: [min(r[f] for r in selected), max(r[f] for r in selected)]
            for f in ("plug_age_average_remaining", "plug_age_average_remaining_shell")
        }
    return curves, summary_rows, {
        "radius_m": radius, "catalyst_per_envelope_volume_kg_m3": cat_density,
        "tortuosity": tau, "water_diffusivity_m2_s": d, "water_film_coefficient_m_s": k,
        "maximum_grid_difference_fraction": max_grid, "maximum_mass_balance_relative_error": max_balance,
        "film_flow_ml_min": plan["film_flow_ml_min"],
        "slowest_flow_water_memory_ranges": memory,
        "scope": "Analogue local-equilibrium storage plus liquid diffusion. All sorbed water counted as extra storage; possible overlap with liquid-filled pore water prolongs the estimate. Water storage is not an enzyme-rate history. Uniform boundary step, then an illustrative plug-displacement age distribution at constant slow-flow film.",
    }


def scaling_cases(reuse, void_fraction):
    # Illustrative histories share an independently known geometric rate
    # convention (eps), not an inferred activity from the bed measurements.
    independent = {
        "central_film": 1.,
        "acetate_reuse_time_proxy": reuse["acetate_mild_time_scale"],
        "benzyl_alcohol_contact_proxy": reuse["propanoic_assay_activity_relative_to_fresh"]["after_benzyl_alcohol"],
        "acid_1_5_reuse_proxy": reuse["propanoic_assay_activity_relative_to_fresh"]["after_second_cycle_acid_alcohol_1_5"],
        "neat_acid_vs_alcohol_proxy": reuse["propanoic_acid_vs_benzyl_alcohol_activity_ratio"],
    }
    illustrations = {
        "basis_30min_illustration": void_fraction,
        "continued_decay_180min_illustration": void_fraction**6,
        "recovered_wetting_180min_illustration": 1 - (1-void_fraction)**6,
    }
    return independent, illustrations


def bed_predictions(jl, reactor, liquid, factors):
    model = jl.LorenzoPackedBed
    p = model.PingPong(np.array([reactor["rate_parameters"][name] for name in published.PARAMETERS]))
    bed = reactor["bed"]
    mass, void, diameter = (bed[k] for k in ("catalyst_mass_kg", "bed_porosity", "bead_mean_diameter_m"))
    if mass != 0.006:
        raise ValueError("Keep the full printed 6 g catalyst mass.")
    area = pi * bed["diameter_m"]**2 / 4
    liquid_volume = area * bed["length_m"] * void
    a_l = 6*(1-void)/(diameter*void)
    pe = bed["length_m"]*0.20/(void**2*diameter)
    acid, alcohol = bed["feed_acid_mol_per_L"], bed["feed_alcohol_mol_per_L"]
    by_temp = {round(row["temperature_K"]-273.15): row for row in liquid["conditions"]}
    rows = []
    grid, balance, film_balance = 0., 0., 0.
    for temp, flow in product(bed["reaction_temperature_C"], bed["feed_flow_ml_per_min"]):
        liq = by_temp[temp]
        k = jl.BeadFilm.film_coefficients(
            flow/60e6/area, diameter, liq["density_kg_m3"], liq["viscosity_Pa_s"],
            liq["diffusivity_m2_s"]["central"], void, "wilson_geankoplis",
        ).coefficient_m_s
        row = {"temperature_C": temp, "flow_ml_min": flow, "catalyst_mass_kg": mass}
        for name, factor in factors.items():
            eta = reactor["effectiveness_by_temperature"][str(temp)]["eta_internal"] * factor
            restricted = jl.BeadFilm.FilmKinetics(p, 60*float(k)*a_l*acid, eta)
            coarse, fine = [model.packed_bed(
                restricted, temp, liquid_volume*1e6/flow, mass/liquid_volume, eta, pe,
                acid0=acid, alcohol0=alcohol, cells=n,
            ) for n in (512, 1024)]
            local = jl.BeadFilm.inlet_balance(restricted, temp, mass/liquid_volume, acid0=acid, alcohol0=alcohol)
            row[name] = float(fine.conversion)
            grid = max(grid, abs(float(coarse.conversion)-row[name]))
            balance = max(balance, float(coarse.balance_error), float(fine.balance_error))
            film_balance = max(film_balance, float(local.film_balance_error_mol_L_min))
        rows.append(row)
    if grid > 1e-4 or balance > 1e-9 or film_balance > 1e-10:
        raise RuntimeError("The fixed-mass scenario grid or reaction/film balance failed.")
    baseline = {published.condition(row): float(row["wilson_geankoplis_central"]) for row in read_rows(FILM_REFERENCE)}
    replay = max(abs(row["central_film"]-baseline[published.condition(row)]) for row in rows)
    if replay > 1e-10:
        raise RuntimeError("The unscaled prediction no longer reproduces the earlier central film.")
    return rows, {"max_axial_grid_change": grid, "max_normalized_acid_balance": balance,
                  "max_film_balance_mol_L_min": film_balance,
                  "max_prior_film_replay_difference": replay}


def compare_bed(rows, cases):
    observed_rows = read_rows(published.HOLDOUT)
    measured = {published.condition(row): float(row["acid_conversion"]) for row in observed_rows}
    if len(rows) != 20 or len(observed_rows) != 20 or len(measured) != 20 or set(measured) != set(map(published.condition, rows)):
        raise ValueError("Use exactly the frozen 20-case holdout grid.")
    y = np.array([measured[published.condition(row)] for row in rows])
    if not np.all(np.isfinite(y) & (y >= 0) & (y <= 1)):
        raise ValueError("Expected finite measured fractional conversions.")
    results = [row | {"observed_acid_conversion": value} for row, value in zip(rows, y, strict=True)]
    comparison = {name: error_metrics(np.array([row[name] for row in rows])-y) for name in cases}
    return results, comparison


def run(output):
    inputs, reactor, liquid = (read_json(path) for path in (INPUT, published.MODEL_INPUTS, film.PROPERTIES))
    particle, water_inputs, phase_data = (read_json(path) for path in (PARTICLES, water.INPUT, phase.PROPERTIES))
    source = read_json(mixing.EXTRACTION)["water_tracer"]
    tracer = tracer_check(source, read_rows(mixing.TRACER))
    reuse = reuse_proxies(inputs)
    independent, illustration = scaling_cases(reuse, reactor["bed"]["bed_porosity"])
    storage, aw, capacity = batch_water(inputs, water_inputs, phase_data, reactor["bed"])
    jl = runtime()
    for name, path in (("LorenzoPackedBed", published.MODEL), ("BeadFilm", film.HERE/"model.jl"),
                       ("BeadMoisture", HERE/"moisture.jl")):
        if not jl.seval(f"isdefined(@__MODULE__, :{name})"):
            jl.include(str(path))
    curves, hydration, moisture = moisture_screen(jl, inputs, water_inputs, liquid, reactor, particle)
    predictions, diagnostics = bed_predictions(jl, reactor, liquid, independent | illustration)
    output.mkdir(parents=True, exist_ok=False)
    write_rows(output/"batch_storage.csv", storage)
    write_rows(output/"batch_activity.csv", aw)
    write_rows(output/"hydration_summary.csv", hydration)
    write_rows(output/"hydration_curves.csv", curves)
    # Only score the reactor after saving all predictions. Even the
    # independently chosen magnitude illustrations are not scored for fit.
    write_rows(output/"predictions_before_comparison.csv", predictions)
    results, scores = compare_bed(predictions, independent)
    write_rows(output/"bed_comparison.csv", results)
    drifts = {
        case: [100*(row[case]-row["basis_30min_illustration"]) for row in predictions]
        for case in illustration if case != "basis_30min_illustration"
    }
    summary = {
        "water_tracer": tracer, "reuse_proxies": reuse, "bed_time_trace": trace_stability(inputs),
        "batch_water_capacity": capacity, "bead_moisture": moisture,
        "independent_activity_scales": independent, "magnitude_illustrations": illustration,
        "illustrative_damaged_fraction_for_nominal_basis_scale":
            (1-reactor["bed"]["bed_porosity"])/(1-independent["neat_acid_vs_alcohol_proxy"]),
        "quasisteady_30_to_180min_conversion_drift_range_pp": {
            name: [min(values), max(values)] for name, values in drifts.items()
        },
        "reaction_numerics": diagnostics, "comparison": scores,
        "sixty_C_three_ml_min": next(row for row in results if published.condition(row) == (60, 3)),
        "julia": str(jl.seval("VERSION")),
        "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in (
            INPUT, PARTICLES, water.INPUT, water.HERE/"run.py", water.HERE/"activity.py",
            phase.PROPERTIES, phase.HERE/"run.py", batch.BATCH_RATES, batch.HERE/"run.py",
            film.PROPERTIES, film.HERE/"model.jl", FILM_REFERENCE,
            published.MODEL_INPUTS, published.MODEL, published.HOLDOUT,
            mixing.TRACER, mixing.EXTRACTION, HERE/"run.py", HERE/"moisture.jl",
        )},
        "output_sha256": {path.name: sha(path) for path in sorted(output.glob("*.csv"))},
    }
    write_json(output/"summary.json", summary)
    return summary


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    parser.add_argument("--output", type=Path, default=OUTPUT/stamp)
    args = parser.parse_args()
    data = run(args.output.resolve())
    print("Output:", args.output)
    print("Aqueous exchange:", data["water_tracer"]["exchange"])
    for name, score in data["comparison"].items():
        print(name, f'held-out RMSE {score["rmse_percentage_points"]:.2f} pp')
