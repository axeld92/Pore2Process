"""Test mass proportionality on whole held-out batch loadings, then screen transfer."""

from argparse import ArgumentParser
from datetime import datetime, timezone
from itertools import product
from math import pi
from pathlib import Path

import numpy as np

from examples.biocatalysis.pentyl_propanoate import run as base
from examples.biocatalysis.pentyl_propanoate.run import array, error_metrics, read_rows, write_rows
from examples.biocatalysis.pentyl_propanoate_film_water.film import run as film
from examples.biocatalysis.pentyl_propanoate_transfer.mixing import run as mixing
from examples.biocatalysis.pentyl_propanoate_transfer.published import run as published
from pore2process_demo.bridge import runtime
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT


HERE = Path(__file__).resolve().parent
PLAN = HERE / "protocol.json"
OUTPUT = ROOT / "results/biocatalysis/pentyl_propanoate_loading"
TRACER = Path(mixing.__file__).resolve().parent / "reference/summary.json"
FILM_REFERENCE = film.HERE / "reference/bed_comparison.csv"
MODELS = ("proportional", "power", "crowding")
VARIANTS = {
    "printed": ("printed", "nominal_ratio", False),
    "equal_curve_weights": ("printed", "nominal_ratio", True),
    "measured_initial_ratio": ("printed", "printed_ratio", False),
    "rederived_acid": ("from_acid", "nominal_ratio", False),
    "rederived_ester": ("from_ester", "nominal_ratio", False),
}


def batch_rows(source_rows, plan):
    source = plan["batch_source"]
    acid_g = source["nominal_acid_volume_ml"] * source["acid_density_at_25_C_g_ml"]
    liquid_L = source["nominal_liquid_volume_ml"]/1000
    groups = {}
    for assay, expected in zip(source["assays"], source["catalyst_loadings_pct_of_initial_acid_mass"], strict=True):
        rows = sorted([r for r in source_rows if int(r["assay_id"]) == assay], key=lambda r: float(r["time_min"]))
        t, loading, T = (array(rows, name) for name in ("time_min", "catalyst_mass_pct_of_initial_acid", "temperature_C"))
        if not (len(t) >= 5 and t[0] == 0.0 and np.all(np.diff(t) > 0)
                and np.all(loading == expected) and np.all(T == source["temperature_C"])):
            raise ValueError(f"Assay {assay} does not match the frozen batch design.")
        initial_acid = float(rows[0]["x_acid"])
        if not 0 < initial_acid < 1 or not 0 < float(rows[0]["x_alcohol"]) < 1:
            raise ValueError("Expected positive initial acid and alcohol mole fractions.")
        initial_ratio = float(rows[0]["x_alcohol"])/initial_acid
        mass_g = acid_g * expected / 100
        exposure_per_min = expected/100 * source["acid_molar_mass_kg_mol"] * 60
        groups[expected] = [{
            "assay": assay, "loading_pct": expected, "mass_g": mass_g,
            "nominal_concentration_g_L": mass_g/liquid_L, "time_min": float(ti),
            "exposure_kg_s_mol_acid": exposure_per_min*ti,
            "nominal_ratio": source["nominal_alcohol_acid_molar_ratio"],
            "printed_ratio": initial_ratio, "printed": float(row["acid_conversion"]),
            "from_acid": 1-float(row["x_acid"])/initial_acid,
            "from_ester": float(row["x_ester"])/initial_acid,
        } for row, ti in zip(rows, t, strict=True)]
        for response in ("printed", "from_acid", "from_ester"):
            values = array(groups[expected], response)
            if not (values[0] == 0 and np.all(np.isfinite(values) & (values >= 0) & (values <= 1))):
                raise ValueError(f"Nonphysical conversion response {response} in assay {assay}.")
    return groups


def crossing_time(time, conversion, target):
    """First upward crossing, linearly interpolated between measured values."""
    time, conversion = np.asarray(time), np.asarray(conversion)
    if not (len(time) == len(conversion) and len(time) >= 2 and np.all(np.diff(time) > 0)
            and np.all(np.isfinite(time)) and np.all(np.isfinite(conversion))
            and np.isfinite(target) and conversion[0] <= target):
        raise ValueError("Expected time-ordered physical conversion records and a covered target.")
    if target == conversion[0]:
        return float(time[0]), [float(time[0]), float(time[0])]
    crosses = np.flatnonzero(conversion >= target)
    if len(crosses) == 0:
        raise ValueError(f"No observation reaches conversion {target}.")
    upper = crosses[0]
    lower = upper-1
    fraction = (target-conversion[lower])/(conversion[upper]-conversion[lower])
    return float(time[lower] + fraction*(time[upper]-time[lower])), [float(time[lower]), float(time[upper])]


def matched_conversion(groups, plan):
    reference = plan["batch_source"]["reference_loading_pct"]
    uncertainty = plan["batch_source"]["reported_conversion_absolute_uncertainty"]
    molar_mass = plan["batch_source"]["acid_molar_mass_kg_mol"]
    levels = plan["matched_conversion_levels"]
    records = []
    for response in ("printed", "from_acid", "from_ester"):
        for loading, rows in groups.items():
            time, conversion = array(rows, "time_min"), array(rows, response)
            previous_time, previous_progress = 0., 0.
            for target in levels:
                t, bracket = crossing_time(time, conversion, target)
                exposure_per_min = loading/100*molar_mass*60
                records.append({
                    "response": response, "loading_pct": loading, "conversion": target,
                    "time_min": t, "bracket_low_min": bracket[0], "bracket_high_min": bracket[1],
                    "time_at_target_minus_reported_conversion_margin_min": (
                        crossing_time(time, conversion, target-uncertainty)[0] if response == "printed" else None),
                    "time_at_target_plus_reported_conversion_margin_min": (
                        crossing_time(time, conversion, target+uncertainty)[0] if response == "printed" else None),
                    "mass_normalized_time_pct_min": loading*t,
                    "exposure_kg_s_mol_acid": exposure_per_min*t,
                    "interval_start_conversion": previous_progress,
                    "interval_mass_specific_rate_mol_kg_s":
                        (target-previous_progress)/(exposure_per_min*(t-previous_time)),
                })
                previous_time, previous_progress = t, target
    summaries = {}
    for response in ("printed", "from_acid", "from_ester"):
        summaries[response] = {}
        for target in levels:
            rows = [r for r in records if r["conversion"] == target and r["response"] == response]
            times = {r["loading_pct"]: r["time_min"] for r in rows}
            u, y = np.log(list(times)), np.log(list(times.values()))
            alpha = -np.dot(u-u.mean(), y-y.mean())/np.dot(u-u.mean(), u-u.mean())
            anchor = next(r for r in rows if r["loading_pct"] == reference)
            for r in rows:
                r["cumulative_mass_specific_rate_relative_to_10pct"] = (
                    anchor["exposure_kg_s_mol_acid"]/r["exposure_kg_s_mol_acid"])
                r["interval_mass_specific_rate_relative_to_10pct"] = (
                    r["interval_mass_specific_rate_mol_kg_s"]/anchor["interval_mass_specific_rate_mol_kg_s"])
            summaries[response][str(target)] = {
                "time_min_by_loading": {str(w): t for w,t in times.items()},
                "cumulative_mass_specific_rate_relative_to_10pct": {
                    str(r["loading_pct"]): r["cumulative_mass_specific_rate_relative_to_10pct"] for r in rows},
                "interval_mass_specific_rate_relative_to_10pct": {
                    str(r["loading_pct"]): r["interval_mass_specific_rate_relative_to_10pct"] for r in rows},
                "descriptive_total_rate_power_exponent": float(alpha),
                "time_based_speedup_10_to_20pct": times[10]/times[20],
            }
    return records, summaries


def time_based_power(groups, levels, response="printed"):
    # Estimate a separate progress intercept at each target conversion,
    # with one common exponent. Only the supplied groups enter the fit.
    loading = np.log(list(groups))
    centered = loading-loading.mean()
    denom = np.dot(centered, centered)
    logs = []
    for X in levels:
        times = np.array([crossing_time(array(r, "time_min"), array(r, response), X)[0] for r in groups.values()])
        logs.append(np.log(times))
    estimate = -sum(np.dot(centered, y-y.mean()) for y in logs)/(denom*len(levels))
    return float(estimate)


def liquid_sample_screen(groups, plan):
    data = plan["batch_source"]
    volume, sample = data["nominal_liquid_volume_ml"], data["reported_sample_volume_ml"]
    out = {}
    for loading,rows in groups.items():
        t = array(rows,"time_min")
        # The volume before the first sample is the nominal charge. Each
        # later interval follows 0.1 mL additional withdrawal. This is a
        # scenario of no solids loss, no replenishment, and no t=0 sample.
        fraction = 1-sample*np.arange(len(t)-1)/volume
        if min(fraction) <= 0:
            raise ValueError("Nominal sampling would exhaust the liquid.")
        effective_time = float(np.sum(np.diff(t)/fraction))
        out[str(loading)] = {
            "max_volume_removed_before_final_sample_fraction": float(1-fraction[-1]),
            "maximum_concentration_increase_factor": float(1/fraction[-1]),
            "cumulative_proportional_exposure_increase_factor": effective_time/t[-1],
        }
    return out


def curve_alignment_check(groups, plan):
    levels = plan["matched_conversion_levels"]
    reference = plan["batch_source"]["reference_loading_pct"]
    folds, records = {}, []
    for held_out, rows in groups.items():
        other = {w: curve for w, curve in groups.items() if w != held_out}
        alpha = time_based_power(other, levels)
        pred, actual = {name: [] for name in ("proportional", "power")}, []
        for row in rows[1:]:
            predictions = {}
            for name, exponent in (("proportional", 1.), ("power", alpha)):
                t = row["time_min"] * (held_out/reference)**exponent
                curves = [(array(curve, "time_min")*(w/reference)**exponent,
                           array(curve, "printed")) for w, curve in other.items()]
                predictions[name] = (float(np.mean([np.interp(t, a, y) for a,y in curves]))
                                     if all(a[0] <= t <= a[-1] for a,y in curves) else None)
            # Compare only points inside every training curve for both
            # scales. This does not extrapolate the low-loading curves.
            shared_support = all(value is not None for value in predictions.values())
            records.append({
                "held_out_loading_pct": held_out, "time_min": row["time_min"],
                "observed_conversion": row["printed"], "training_power_exponent": alpha,
                "in_common_support_of_both_models": shared_support,
                **predictions,
            })
            if shared_support:
                actual.append(row["printed"])
                for name,value in predictions.items():
                    pred[name].append(value)
        if not actual:
            raise ValueError("A held-out curve has no common exposure support.")
        folds[str(held_out)] = {
            "training_power_exponent": alpha, "compared_points": len(actual),
            "held_out_post_start_points": len(rows)-1,
            "scores": {law: error_metrics(np.asarray(values)-actual) for law, values in pred.items()},
        }
    return records, {
        "all_loadings_descriptive_power_exponent": time_based_power(groups, levels),
        "folds": folds,
        "equal_loading_rmse_pp_on_common_support": {
            law: float(np.sqrt(np.mean([f["scores"][law]["rmse_percentage_points"]**2 for f in folds.values()])))
            for law in ("proportional", "power")},
        "scope": "Independent, shape-free time alignment using only other curves and targets. Scores compare the same supported holdout points, not extrapolated tails.",
    }


def vectors(groups, loadings, response, ratio, balanced):
    selected = [row for w in loadings for row in groups[w][1:]]
    weights = np.array([1/(len(groups[row["loading_pct"]])-1) if balanced else 1. for row in selected])
    weights *= len(weights)/weights.sum()
    return selected, (array(selected,"loading_pct"), array(selected,"time_min"), array(selected, ratio),
                      array(selected,response), weights)


def fit(model, law, data, reference, fixed_K=0., points=25):
    result = model.fit_loading(law, *data, reference_loading_pct=reference, fixed_K=fixed_K, profile_points=points)
    gradient = float(result.gradient_inf)
    record = {
        "law": law, "shape": float(result.shape), "k_ref_mol_kg_s": float(result.k_ref),
        "apparent_Kx": float(result.Kx), "weighted_sse": float(result.weighted_sse),
        "solver": str(result.solver), "conditional_gradient_inf": gradient if np.isfinite(gradient) else None,
    }
    if (not 0.0001 < record["k_ref_mol_kg_s"] < 10 or not 0.001 < record["apparent_Kx"] < 1000
            or record["solver"] != "least_squares" or gradient > 1e-8):
        raise RuntimeError(f"Selected loading fit did not resolve its nuisance parameters: {record}")
    return record


def batch_predictions(model, fit, data, reference):
    w,t,ratio,_,_ = data
    return np.array(model.predict_batch(fit["law"], fit["shape"], fit["k_ref_mol_kg_s"],
                                       fit["apparent_Kx"], w,t,ratio, reference_loading_pct=reference))


def qualification(folds, baseline, alternative, plan):
    b = [case[baseline]["rmse_percentage_points"] for case in folds.values()]
    a = [case[alternative]["rmse_percentage_points"] for case in folds.values()]
    rmse_b = float(np.sqrt(np.mean(np.square(b))))
    rmse_a = float(np.sqrt(np.mean(np.square(a))))
    improvement = 1-rmse_a/rmse_b
    worst = max(np.asarray(a)-b)
    limits = plan["cross_validation"]["qualification_rule"]
    passed = (improvement >= limits["relative_improvement_in_equal_loading_rmse_min"]
              and worst <= limits["maximum_worsening_of_any_fold_rmse_pp"])
    return {
        "equal_loading_rmse_pp": rmse_a, "relative_rmse_improvement": improvement,
        "largest_fold_rmse_worsening_pp": float(worst),
        "folds_improved": sum(x < y for x,y in zip(a,b,strict=True)),
        "passes_screen": bool(passed),
    }


def fit_batch_models(model, groups, plan):
    reference = plan["batch_source"]["reference_loading_pct"]
    all_loadings = list(groups)
    fits, predictions, fold_rows, full_records = [], [], [], {}
    variants, maximum_refinement = {}, 0.
    for name, (response, ratio, balanced) in VARIANTS.items():
        full_records[name] = {}
        folds, all_case_fits = {}, {}
        for held_out in (None,*all_loadings):
            train_loadings = [w for w in all_loadings if w != held_out]
            _, training = vectors(groups, train_loadings, response, ratio, balanced)
            eval_loadings = all_loadings if held_out is None else [held_out]
            eval_rows, validation = vectors(groups, eval_loadings, response, ratio, False)
            fold_name = "all_data" if held_out is None else str(held_out)
            fit_scores = {}
            all_case_fits[fold_name] = {}
            for law in MODELS:
                record = fit(model, law, training, reference)
                all_case_fits[fold_name][law] = record
                prediction = batch_predictions(model, record, validation, reference)
                if name == "printed" and held_out is None and law != "proportional":
                    refined = fit(model, law, training, reference, points=49)
                    change = np.max(abs(prediction-batch_predictions(model, refined, validation, reference)))
                    maximum_refinement = max(maximum_refinement, float(change))
                    if change > 1e-5:
                        raise RuntimeError("The loading-parameter profile did not converge.")
                fits.append({"variant": name, "held_out_loading_pct": held_out,
                             "training_post_start_points": len(training[0]), **record})
                errors = prediction - validation[3]
                fit_scores[law] = error_metrics(errors)
                for row,pred,observation in zip(eval_rows,prediction,validation[3],strict=True):
                    predictions.append({
                        "variant": name, "held_out_loading_pct": held_out,
                        "law": law, "loading_pct": row["loading_pct"], "time_min": row["time_min"],
                        "observed_conversion": observation, "predicted_conversion": pred,
                    })
                if held_out is None:
                    full_records[name][law] = record | {"in_sample_scores": fit_scores[law]}
                else:
                    fold_rows.append({"variant": name, "held_out_loading_pct": held_out,
                                      "law": law, **fit_scores[law]})
            if held_out is not None:
                folds[str(held_out)] = fit_scores
        primary_rmse = np.sqrt(np.mean([f["proportional"]["rmse_percentage_points"]**2 for f in folds.values()]))
        heldout_summary = {"proportional": {"equal_loading_rmse_pp": float(primary_rmse)}}
        for law in MODELS[1:]:
            test = qualification(folds, "proportional", law, plan)
            heldout_summary[law] = test
        variants[name] = {"full_fit": full_records[name], "fold_scores": folds,
                          "holdout_summary": heldout_summary, "fits_by_fold": all_case_fits}
    # Fixed published equilibrium is a sensitivity of the full fit only.
    # It is intentionally never supplied in the strict held-out folds.
    _, all_data = vectors(groups, all_loadings, "printed", "nominal_ratio", False)
    fixed_K = plan["reference_common_apparent_Kx_for_sensitivity"]
    fixed = {}
    for law in MODELS:
        record = fit(model, law, all_data, reference, fixed_K=fixed_K)
        fixed[law] = record | {"in_sample_scores":
            error_metrics(batch_predictions(model, record, all_data, reference)-all_data[3])}
    per_loading = {}
    for loading in all_loadings:
        _, selected = vectors(groups, [loading], "printed","nominal_ratio", False)
        record = fit(model, "proportional", selected, reference, fixed_K=fixed_K)
        per_loading[str(loading)] = record | {
            "observed_count": len(selected[3]),
            "in_sample_scores": error_metrics(batch_predictions(model, record, selected, reference)-selected[3]),
        }
    high, low = sorted(all_loadings)[-2:], sorted(all_loadings)[:-2]
    _, high_training = vectors(groups, high, "printed", "nominal_ratio", False)
    _, low_validation = vectors(groups, low, "printed", "nominal_ratio", False)
    high_pair = {}
    for law in MODELS:
        record = fit(model, law, high_training, reference)
        y = batch_predictions(model,record,low_validation,reference)
        errors = y-low_validation[3]
        high_pair[law] = record | {
            "high_pair_in_sample": error_metrics(batch_predictions(model,record,high_training,reference)-high_training[3]),
            "lower_loadings_held_out": {
                str(w): error_metrics(errors[low_validation[0]==w]) for w in low},
            "combined_lower_loading_score": error_metrics(errors),
        }
    return fits, predictions, fold_rows, {
        "variants": variants, "fixed_published_equilibrium_full_fit_sensitivity": fixed,
        "diagnostic_individual_loading_at_common_published_equilibrium": per_loading,
        "high_pair_fit_lower_loadings_held_out": {
            "training_loadings_pct": high, "held_out_loadings_pct": low,
            "cases": high_pair, "scope": "Closest measured concentrations alone; no reactor fit or qualification."
        },
        "maximum_full_primary_profile_refinement_absolute_conversion": maximum_refinement,
        "profile_points_primary_refined": [25,49],
    }


def transfer_factors(model, groups, source, plan, batch_fit):
    bed = source["bed"]
    if bed["catalyst_mass_kg"] != plan["transfer"]["bed_catalyst_mass_kg"]:
        raise ValueError("Keep the full 6 g printed inventory.")
    volume = pi * bed["diameter_m"]**2 / 4 * bed["length_m"]
    liquid = volume*bed["bed_porosity"]
    ref_load = plan["batch_source"]["reference_loading_pct"]
    ref = groups[ref_load][0]["nominal_concentration_g_L"]
    concentrations = {"void_liquid_basis": bed["catalyst_mass_kg"]/liquid,
                      "total_volume_proxy": bed["catalyst_mass_kg"]/volume}
    scales = {"mass_proportional": 1.}
    scenario_scope = {"mass_proportional": "reference"}
    sensitivity = {}
    for law, (basis, concentration) in product(MODELS[1:], concentrations.items()):
        case = f"{law}_{basis}"
        r = concentration/ref
        primary = batch_fit["variants"]["printed"]
        shape = primary["full_fit"][law]["shape"]
        scales[case] = float(model.rate_scale(law, shape, r))
        scenario_scope[case] = "all four 40 C batch loadings"
        from_folds = [float(model.rate_scale(law, fold[law]["shape"], r))
                      for key,fold in primary["fits_by_fold"].items() if key != "all_data"]
        from_variants = {name: float(model.rate_scale(law, variant["full_fit"][law]["shape"], r))
                         for name,variant in batch_fit["variants"].items()}
        high_pair_scale = float(model.rate_scale(
                law, batch_fit["high_pair_fit_lower_loadings_held_out"]["cases"][law]["shape"], r))
        diagnostic = f"high_pair_diagnostic_{case}"
        scales[diagnostic] = high_pair_scale
        scenario_scope[diagnostic] = "two highest loadings only; unqualified extrapolation"
        sensitivity[case] = {
            "primary_full_fit": scales[case], "range_from_training_folds": [min(from_folds),max(from_folds)],
            "full_fit_response_and_weight_sensitivities": from_variants,
            "full_fit_published_equilibrium_sensitivity": float(model.rate_scale(
                law, batch_fit["fixed_published_equilibrium_full_fit_sensitivity"][law]["shape"], r)),
            "high_pair_only_extrapolation_unqualified_diagnostic": high_pair_scale,
        }
    return scales, {
        "catalyst_mass_kg": bed["catalyst_mass_kg"], "cartridge_volume_ml": volume*1e6,
        "interbead_liquid_volume_ml": liquid*1e6, "batch_reference_concentration_g_L": ref,
        "largest_batch_concentration_g_L": max(groups[w][0]["nominal_concentration_g_L"] for w in groups),
        "bed_concentrations_g_L": concentrations, "factor_scenarios": scales,
        "factor_calibration_scope": scenario_scope,
        "factors_from_full_fits_and_batch_training_folds": sensitivity,
        "scope": "Extrapolation outside measured concentrations; all temperature-dependent kinetics independently anchored, relative loading trend held at 40 C. No alternative is accepted by a reactor score.",
    }


def anchor_batch_kinetics(model, source_rows, plan):
    T = {int(r["temperature_C"]):float(r["apparent_equilibrium_constant_Kx"])
         for r in read_rows(base.EQUILIBRIUM)}
    reference = plan["batch_source"]["reference_loading_pct"]
    anchors, records = {}, {}
    for temp,K in T.items():
        rows = [r for r in source_rows if float(r["temperature_C"]) == temp and
                float(r["catalyst_mass_pct_of_initial_acid"]) == reference and float(r["time_min"]) > 0]
        if len(rows) < 4:
            raise ValueError("The temperature-specific anchor requires 10% batch data.")
        fit = model.fit_at_temperature(temp, array(rows,"catalyst_mass_pct_of_initial_acid"),
            array(rows,"time_min"), array(rows,"acid_conversion"), K)
        anchors[temp] = fit.kinetics
        records[str(temp)] = base.fit_record(model,fit,temp) | {
            "batch_scores": error_metrics(np.array(fit.predicted)-array(rows,"acid_conversion"))}
    return anchors, records


def predict_bed(jl, source, film_props, factors, anchors, plan, pe):
    reactor, film_model, transfer = jl.LorenzoPackedBed, jl.BeadFilm, jl.LoadingFilmTransfer
    kinetics = reactor.PingPong(np.array([source["rate_parameters"][key] for key in published.PARAMETERS]))
    bed = source["bed"]
    mass, void, diameter = (bed[k] for k in ("catalyst_mass_kg","bed_porosity","bead_mean_diameter_m"))
    if mass != 0.006:
        raise ValueError("Preserve the full printed inventory.")
    area = pi*bed["diameter_m"]**2/4
    liquid = area*bed["length_m"]*void
    conc = mass/liquid
    a_l = 6*(1-void)/(diameter*void)
    film_pe = bed["length_m"]*0.20/(void**2*diameter)
    acid, alcohol = bed["feed_acid_mol_per_L"], bed["feed_alcohol_mol_per_L"]
    limits = plan["transfer"]
    by_T = {round(row["temperature_K"]-273.15):row for row in film_props["conditions"]}
    predictions = []
    batch_grid, film_grid, acid_balance, local_balance = 0.,0.,0.,0.
    for temp,flow in product(bed["reaction_temperature_C"],bed["feed_flow_ml_per_min"]):
        liq = by_T[temp]
        k = film_model.film_coefficients(flow/60e6/area, diameter, liq["density_kg_m3"],
            liq["viscosity_Pa_s"],liq["diffusivity_m2_s"]["central"],void,"wilson_geankoplis").coefficient_m_s
        eta = source["effectiveness_by_temperature"][str(temp)]["eta_internal"]
        inner = film_model.FilmKinetics(kinetics,60*float(k)*a_l*acid,eta)
        row = {"temperature_C": temp,"flow_ml_min": flow,"actual_catalyst_mass_kg": mass,
               "actual_catalyst_concentration_on_liquid_basis_g_L": conc}
        for case,scale in factors.items():
            params = np.array(jl.PentylPropanoate.coefficients(anchors[temp]))
            params[0] += np.log(scale)  # Only the per-mass rate constant changes.
            scaled = jl.PentylPropanoate.Kinetics(*params)
            args = (scaled,temp,mass,acid*flow/60000)
            coarse,fine = [jl.PentylPropanoate.bed_prediction(
                *args,peclet=pe,cells=n) for n in limits["bed_grid_batch_model"]]
            row[f"batch_{case}"] = float(fine.conversion)
            batch_grid = max(batch_grid,float(abs(coarse.conversion-fine.conversion)))
            acid_balance = max(acid_balance,float(coarse.balance_error),float(fine.balance_error))
            restricted = transfer.LoadingAdjustedFilm(inner,scale)
            coarse,fine = [reactor.packed_bed(restricted,temp,liquid*1e6/flow,conc,eta,film_pe,
                acid0=acid,alcohol0=alcohol,cells=n) for n in limits["bed_grid_film_model"]]
            row[f"film_{case}"] = float(fine.conversion)
            film_grid = max(film_grid,float(abs(coarse.conversion-fine.conversion)))
            acid_balance = max(acid_balance,float(coarse.balance_error),float(fine.balance_error))
            local = transfer.inlet_balance(restricted,temp,conc,acid0=acid,alcohol0=alcohol)
            local_balance = max(local_balance,float(local.film_balance_error_mol_L_min))
        predictions.append(row)
    if (max(batch_grid,film_grid) > limits["maximum_axial_grid_conversion_difference"]
            or acid_balance > limits["maximum_normalized_balance_error"]
            or local_balance > limits["maximum_film_balance_error_mol_L_min"]):
        raise RuntimeError("Loading-transfer axial grid or conservation/film checks failed.")
    previous = {published.condition(row):float(row["wilson_geankoplis_central"]) for row in read_rows(FILM_REFERENCE)}
    replay = max(abs(row["film_mass_proportional"]-previous[published.condition(row)]) for row in predictions)
    if replay > 1e-10:
        raise RuntimeError("The mass-proportional film predictions changed.")
    return predictions, {
        "max_batch_axial_grid_change": batch_grid, "max_film_axial_grid_change": film_grid,
        "max_normalized_acid_balance": acid_balance, "max_liquid_film_balance_mol_L_min": local_balance,
        "max_previous_central_film_replay_difference": replay,
        "measured_water_mixing_scenario_peclet": pe, "film_reconstruction_peclet": film_pe,
    }


def compare_bed(predictions, factors):
    rows = read_rows(published.HOLDOUT)
    measured = {published.condition(row):float(row["acid_conversion"]) for row in rows}
    if len(rows) != 20 or len(measured) != 20 or set(measured) != set(map(published.condition,predictions)):
        raise ValueError("Compare all 20 frozen reactor observations.")
    y = np.array([measured[published.condition(row)] for row in predictions])
    if not np.all(np.isfinite(y)&(y>=0)&(y<=1)):
        raise ValueError("Use finite measured fractional conversions.")
    results = [dict(row,observed_conversion=float(value)) for row,value in zip(predictions,y,strict=True)]
    cases = [f"{branch}_{case}" for branch,case in product(("batch","film"),factors)]
    metrics = {case:error_metrics(array(results,case)-y) for case in cases}
    return results, metrics


def run(output):
    plan = read_json(PLAN)
    source_rows = read_rows(base.BATCH)
    groups = batch_rows(source_rows,plan)
    matched,matched_summary = matched_conversion(groups,plan)
    aligned,aligned_summary = curve_alignment_check(groups,plan)
    jl = runtime()
    for name,path in (("PentylPropanoate",base.HERE/"model.jl"),("CatalystLoading",HERE/"model.jl"),
        ("LorenzoPackedBed",published.MODEL),("BeadFilm",film.HERE/"model.jl"),
        ("LoadingFilmTransfer",HERE/"film_transfer.jl")):
        if not jl.seval(f"isdefined(@__MODULE__, :{name})"):
            jl.include(str(path))
    fit_rows,prediction_rows,cv,batch_fit = fit_batch_models(jl.CatalystLoading,groups,plan)
    reactor,liquid = read_json(published.MODEL_INPUTS),read_json(film.PROPERTIES)
    factors,transfer = transfer_factors(jl.CatalystLoading,groups,reactor,plan,batch_fit)
    anchors,anchor_records = anchor_batch_kinetics(jl.PentylPropanoate,source_rows,plan)
    pe = read_json(TRACER)["conditional_danckwerts_transfer"]["nominal_rectangular_pulse"]["danckwerts_peclet"]
    bed_predictions,numerics = predict_bed(jl,reactor,liquid,factors,anchors,plan,pe)
    output.mkdir(parents=True,exist_ok=False)
    for name,rows in (
        ("matched_conversion.csv",matched),("shape_free_holdout.csv",aligned),
        ("batch_fits.csv",fit_rows),("batch_predictions.csv",prediction_rows),("batch_fold_scores.csv",cv),
        ("predictions_before_comparison.csv",bed_predictions)):
        write_rows(output/name,rows)
    # The reactor conversion markers are opened only after the independent
    # fits, loading factors and all full-inventory predictions are frozen.
    results,metrics = compare_bed(bed_predictions,factors)
    write_rows(output/"bed_comparison.csv",results)
    summary = {
        "batch_condition_summary": {str(w): {
            "assay": rows[0]["assay"], "rows_including_start": len(rows),
            "catalyst_g": rows[0]["mass_g"], "nominal_catalyst_concentration_g_L": rows[0]["nominal_concentration_g_L"],
            "printed_initial_alcohol_acid_ratio": rows[0]["printed_ratio"],
        } for w,rows in groups.items()},
        "liquid_only_aliquot_withdrawal_sensitivity": liquid_sample_screen(groups,plan),
        "matched_conversion": matched_summary,"shape_free_alignment": aligned_summary,
        "batch_fit": batch_fit,"bed_transfer": transfer,"ten_percent_batch_anchors": anchor_records,
        "bed_numerics": numerics,"bed_comparison": metrics,
        "sixty_C_three_ml_min": next(row for row in results if published.condition(row)==(60,3)),
        "julia": str(jl.seval("VERSION")),
        "source_sha256": {str(path.relative_to(ROOT)):sha(path) for path in (
            PLAN,HERE/"model.jl",HERE/"film_transfer.jl",HERE/"run.py",base.BATCH,base.EQUILIBRIUM,
            base.HERE/"model.jl",base.HERE/"run.py",published.MODEL,published.MODEL_INPUTS,
            published.HOLDOUT,film.PROPERTIES,film.HERE/"model.jl",FILM_REFERENCE,TRACER,
        )},
        "output_sha256": {path.name:sha(path) for path in sorted(output.glob("*.csv"))},
    }
    write_json(output/"summary.json",summary)
    return summary


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    parser.add_argument("--output",type=Path,default=OUTPUT/stamp)
    args = parser.parse_args()
    report = run(args.output.resolve())
    print("Output:",args.output)
    for law,score in report["batch_fit"]["variants"]["printed"]["holdout_summary"].items():
        print(law, "whole-loading held-out RMSE",f'{score["equal_loading_rmse_pp"]:.3f} pp')
    print("Rate factors:",report["bed_transfer"]["factor_scenarios"])
    for name,score in report["bed_comparison"].items():
        print(name,"held-out reactor RMSE",f'{score["rmse_percentage_points"]:.3f} pp')
