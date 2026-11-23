"""Compare matched-progress loading laws across independent batch studies."""

from argparse import ArgumentParser
from collections import Counter
from datetime import datetime, timezone
from itertools import product
from math import exp, log, pi
from pathlib import Path
import copy

import numpy as np

from examples.biocatalysis.pentyl_propanoate import run as base
from examples.biocatalysis.pentyl_propanoate.run import array, error_metrics, read_rows, write_rows
from examples.biocatalysis.pentyl_propanoate_loading import run as prior
from pore2process_demo.bridge import runtime
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT

HERE = Path(__file__).resolve().parent
PLAN = HERE / "protocol.json"
SOURCE = ROOT / "literature/sources"
BATCH_1 = SOURCE / "lorenzo2026/batch_1to1.csv"
OUTPUT = ROOT / "results/biocatalysis/pentyl_propanoate_loading_evidence"
LAWS = prior.MODELS
VARIANTS = {
    "primary": ("printed", "nominal_ratio", True),
    "point_weights": ("printed", "nominal_ratio", False),
    "printed_initial_ratio": ("printed", "printed_ratio", True),
    "rederived_acid": ("from_acid", "nominal_ratio", True),
    "rederived_ester": ("from_ester", "nominal_ratio", True),
    "published_K_training_anchors": ("printed", "nominal_ratio", True),
}


def pentyl_groups(plan):
    spec = plan["same_chemistry"]
    groups, basis = {}, {}
    for ratio, path in ((2, base.BATCH), (1, BATCH_1)):
        subplan = copy.deepcopy(read_json(prior.PLAN))
        acid_ml = spec["nominal_acid_charge_ml_by_alcohol_acid_ratio"][str(ratio)]
        alcohol_ml = spec["nominal_alcohol_charge_ml_by_alcohol_acid_ratio"][str(ratio)]
        density = spec["pure_liquid_densities_at_25_C_g_ml"]
        acid_g, alcohol_g = acid_ml*density["acid"], alcohol_ml*density["alcohol"]
        subplan["batch_source"].update({
            "assays": spec["assays_by_ratio"][str(ratio)],
            "nominal_acid_volume_ml": acid_ml, "nominal_alcohol_acid_molar_ratio": ratio,
        })
        selected = prior.batch_rows(read_rows(path), subplan)
        basis[ratio] = {"acid_g": acid_g, "alcohol_g": alcohol_g,
                        "nominal_liquid_density_g_L": (acid_g+alcohol_g)*1000/spec["nominal_liquid_charge_ml"]}
        for w,rows in selected.items():
            groups[(ratio,w)] = [row | {"feed": ratio, "loading_pct_total_substrates":
                w*acid_g/(acid_g+alcohol_g)} for row in rows]
    ref = groups[(2,10)][0]["nominal_concentration_g_L"]
    return groups,basis,ref


def crossing_bounds(time, conversion, margin, target):
    """Explicit measurement-reading scenarios, not statistical confidence limits."""
    t, y, m = map(np.asarray, (time,conversion,margin))
    if len(t) != len(y) or len(t) != len(m) or np.any(m < 0) or m[0] != 0:
        raise ValueError("Expected equal-length curve arrays and nonnegative margins, fixed origin.")
    value, bracket = prior.crossing_time(t,y,target)
    bounds = []
    for scenario in (np.minimum(1,y+m), np.maximum(0,y-m)):
        try:
            bounds.append(prior.crossing_time(t,scenario,target)[0])
        except ValueError:
            bounds.append(None)
    return value,bracket,bounds


def progress_records(series, curves, levels, anchor, reference, reported_deviation=None):
    """Curve-level data: reported dosage and a separately audited comparison basis."""
    records = []
    for target in levels:
        matched = []
        for loading,curve in curves.items():
            t,y,margin = (np.array(curve[k], dtype=float) for k in ("time_min","conversion","margin"))
            crossing,bracket,(early,late) = crossing_bounds(t,y,margin,target)
            if series.startswith("lorenzo") is False and bracket[0] == 0:
                raise ValueError(f"{series} target {target} is before the first observed positive time.")
            broad = [None,None]
            if reported_deviation is not None:
                wide = margin.copy()
                wide[1:] += reported_deviation/100
                broad = crossing_bounds(t,y,wide,target)[2]
            comparison = curve["comparison_pct"]
            matched.append({
                "series": series, "loading_pct_reported": loading,
                "loading_pct_total_substrates": curve["total_pct"],
                "relative_loading_for_model": comparison/reference,
                "target_conversion": target, "time_min": crossing,
                "bracket_low_min": bracket[0], "bracket_high_min": bracket[1],
                "time_early_reading_scenario_min": early, "time_late_reading_scenario_min": late,
                "time_early_with_reported_deviation_min": broad[0],
                "time_late_with_reported_deviation_min": broad[1],
                "anchor_loading_pct_reported": anchor,
            })
        a = next(r for r in matched if r["loading_pct_reported"] == anchor)
        for r in matched:
            fraction = anchor/r["loading_pct_reported"]
            r["relative_activity_per_gram"] = fraction*a["time_min"]/r["time_min"]
            for name,early_key,late_key in (
                ("reading", "time_early_reading_scenario_min", "time_late_reading_scenario_min"),
                ("with_reported_deviation", "time_early_with_reported_deviation_min", "time_late_with_reported_deviation_min"),
            ):
                endpoints = (a[early_key],r[late_key],a[late_key],r[early_key])
                lo = (fraction*endpoints[0]/endpoints[1]
                      if endpoints[0] is not None and endpoints[1] is not None else None)
                hi = (fraction*endpoints[2]/endpoints[3]
                      if endpoints[2] is not None and endpoints[3] is not None else None)
                if r is a:
                    lo,hi = 1.,1.
                r[f"relative_activity_{name}_low"] = lo
                r[f"relative_activity_{name}_high"] = hi
            records.append(r)
    return records


def all_progress(plan, groups, paroul_figure=4):
    external,internal = {},{}
    common = plan["external"]["common_reference_pct_of_total_substrates"]
    for series,spec in plan["external"]["series"].items():
        rows = read_rows(SOURCE/series/"loading_progress.csv")
        curves = {}
        for label in spec["loadings"]:
            column = f"conversion_pct_at_{label}"
            if series == "paroul2010" and label == "10" and paroul_figure != 4:
                column += f"_fig{paroul_figure}"
            selected = [r for r in rows if r[column] != ""]
            w = float(label)
            curves[w] = {"time_min": array(selected,"time_min"), "conversion": array(selected,column)/100,
                "margin": array(selected,"marker_reading_margin_pp")/100,
                "comparison_pct": w, "total_pct": w if spec["harmonised"] else None}
        reference = common if spec["harmonised"] else float(spec["anchor"])
        external[series] = progress_records(series,curves,spec["levels"],float(spec["anchor"]),
            reference,spec["reported_overall_conversion_deviation_pp"])
    for ratio in (2,1):
        series = f"lorenzo2026_{ratio}to1"
        curves = {}
        for (r,w), rows in groups.items():
            if r != ratio:
                continue
            margin = np.full(len(rows),plan["same_chemistry"]["reported_conversion_absolute_uncertainty"])
            margin[0] = 0
            total = rows[0]["loading_pct_total_substrates"]
            curves[w] = {"time_min": array(rows,"time_min"), "conversion": array(rows,"printed"),
                "margin": margin, "comparison_pct": total, "total_pct": total}
        internal[series] = progress_records(series,curves,plan["same_chemistry"]["matched_conversion_levels"],
            10.,common)
    return external,internal


def direct_exposure_checks(plan):
    comparisons = []
    for spec in plan["external"]["direct_exposure_checks"]:
        series = spec["series"]
        rows = read_rows(SOURCE/series/"loading_progress.csv")
        a,h = (next(r for r in rows if float(r["time_min"])==spec[key])
               for key in ("anchor_time_min","higher_time_min"))
        xa,xh = (float(row[f"conversion_pct_at_{spec[key]}"])
                 for row,key in ((a,"anchor_loading"),(h,"higher_loading")))
        exposure_a,exposure_h = (float(spec[w])*spec[t] for w,t in (
            ("anchor_loading","anchor_time_min"),("higher_loading","higher_time_min")))
        if exposure_h+1e-12 < exposure_a:
            raise ValueError("The higher-loading curve must have at least the anchor's mass exposure.")
        reported = plan["external"]["series"][series]["reported_overall_conversion_deviation_pp"]
        margin = float(a["marker_reading_margin_pp"])+float(h["marker_reading_margin_pp"])
        comparisons.append(spec | {"anchor_exposure_pct_min":exposure_a,
            "higher_exposure_pct_min":exposure_h,"anchor_conversion_pct":xa,"higher_conversion_pct":xh,
            "conversion_gap_pp":xa-xh,"minimum_gap_marker_reading_scenario_pp":xa-xh-margin,
            "minimum_gap_with_reported_deviation_scenario_pp":(
                xa-xh-margin-2*reported if reported is not None else None)})
    return comparisons


def reading_progress_scenario(external, toward_proportional):
    """Directional bounds on reading, with the reference curve unchanged."""
    result = {}
    for name,rows in external.items():
        shifted = []
        for r in rows:
            difference = r["loading_pct_reported"]-r["anchor_loading_pct_reported"]
            if difference == 0:
                shifted.append(r.copy())
                continue
            earlier = (difference > 0) == toward_proportional
            key = "time_early_reading_scenario_min" if earlier else "time_late_reading_scenario_min"
            t = r[key]
            if t is None:
                raise ValueError("Directional reading scenario lacks a covered crossing.")
            shifted.append(r | {"time_min":t})
        result[name] = shifted
    return result


def log_time_scores(rows):
    error = np.array([log(r["predicted_time_min"]/r["time_min"]) for r in rows])
    return {"count": len(rows), "log_time_rmse": float(np.sqrt(np.mean(error**2))),
            "typical_multiplicative_time_error": exp(float(np.sqrt(np.mean(error**2)))),
            "maximum_multiplicative_time_error": exp(float(np.max(abs(error))))}


def crossing_fit(model, law, rows):
    ids = {key:i+1 for i,key in enumerate(sorted({(r["series"],r["target_conversion"]) for r in rows}))}
    counts = Counter(r["series"] for r in rows)
    fit = model.fit_crossing_law(law, array(rows,"relative_loading_for_model"), array(rows,"time_min"),
        np.array([1/counts[r["series"]] for r in rows]),
        np.array([ids[(r["series"],r["target_conversion"])] for r in rows]))
    return {"law": law, "shape": float(fit.shape), "weighted_log_time_sse": float(fit.weighted_sse)}


def crossing_predict(model, fit, anchors, validation):
    # An intercept uses only the supplied training curves (or the single
    # anchor of an entirely unseen chemistry); no response from validation.
    rate = lambda r: float(model.log_total_rate(fit["law"],fit["shape"],r["relative_loading_for_model"]))
    transformed = {}
    for r in anchors:
        key = (r["series"],r["target_conversion"])
        transformed.setdefault(key,[]).append(log(r["time_min"])+rate(r))
    return [r | {"predicted_time_min": exp(float(np.mean(
        transformed[(r["series"],r["target_conversion"])]))-rate(r))}
        for r in validation]


def compare_progress_laws(model, plan, external, internal):
    fits, folds, predictions = {},[],[]
    for series,rows in external.items():
        per_law = {}
        for law in LAWS:
            fit = crossing_fit(model,law,rows)
            scores = {}
            for w in sorted({r["loading_pct_reported"] for r in rows}):
                training = [r for r in rows if r["loading_pct_reported"] != w]
                test = [r for r in rows if r["loading_pct_reported"] == w]
                fold = crossing_fit(model,law,training)
                values = crossing_predict(model,fold,training,test)
                score = log_time_scores(values)
                scores[str(w)] = score
                folds.append({"scope": "within_study", "held_out": f"{series}:{w}",
                              **fold, **score})
                predictions.extend(r | {"scope": "within_study", "law": law, "shape": fold["shape"]}
                                   for r in values)
            per_law[law] = fit | {"whole_loading_holdouts":scores,
                "equal_loading_log_time_rmse": float(np.sqrt(np.mean(
                    [s["log_time_rmse"]**2 for s in scores.values()])))}
        fits[series] = per_law
    pool = {name:rows for name,rows in external.items()
            if plan["external"]["series"][name]["harmonised"]}
    records = [r for rows in pool.values() for r in rows]
    common = {}
    for law in LAWS:
        full = crossing_fit(model,law,records)
        paper_scores,target_scores = {},{}
        for series,rows in pool.items():
            train = [r for r in records if r["series"] != series]
            fitted = crossing_fit(model,law,train)
            anchor = [r for r in rows if r["loading_pct_reported"] == r["anchor_loading_pct_reported"]]
            test = [r for r in rows if r["loading_pct_reported"] != r["anchor_loading_pct_reported"]]
            values = crossing_predict(model,fitted,anchor,test)
            score = log_time_scores(values)
            paper_scores[series] = score | {"training_shape": fitted["shape"]}
            folds.append({"scope": "held_out_paper", "held_out":series, **fitted, **score})
            predictions.extend(r | {"scope":"held_out_paper","law":law,"shape":fitted["shape"]}
                               for r in values)
        for series,rows in internal.items():
            anchor = [r for r in rows if r["loading_pct_reported"] == r["anchor_loading_pct_reported"]]
            test = [r for r in rows if r["loading_pct_reported"] != r["anchor_loading_pct_reported"]]
            values = crossing_predict(model,full,anchor,test)
            score = log_time_scores(values)
            target_scores[series] = score
            folds.append({"scope":"external_to_pentyl","held_out":series,**full,**score})
            predictions.extend(r | {"scope":"external_to_pentyl","law":law,"shape":full["shape"]}
                               for r in values)
        common[law] = full | {"held_out_external_papers":paper_scores,
            "external_to_pentyl":target_scores,
            "equal_paper_log_time_rmse": float(np.sqrt(np.mean(
                [s["log_time_rmse"]**2 for s in paper_scores.values()])))}
    null = common["proportional"]
    for law in LAWS[1:]:
        entry = common[law]
        improvement = 1-entry["equal_paper_log_time_rmse"]/null["equal_paper_log_time_rmse"]
        worst = max(score["log_time_rmse"]-null["held_out_external_papers"][series]["log_time_rmse"]
                    for series,score in entry["held_out_external_papers"].items())
        target = all(score["log_time_rmse"] < null["external_to_pentyl"][series]["log_time_rmse"]
                    for series,score in entry["external_to_pentyl"].items())
        entry["screen"] = {"external_paper_relative_improvement":improvement,
            "maximum_external_paper_worsening_log_time":worst, "improves_both_pentyl_ratios":target,
            "passes": bool(improvement >= .1 and worst <= .05 and target)}
    return {"separate_studies":fits,"common_external":common},folds,predictions


def select_rows(groups, keys):
    return [r for key in keys for r in groups[key][1:]]


def shared_fit(model,law,rows,ref,variant,*,fixed_shape=None,points=25):
    response, ratio,balanced = VARIANTS[variant]
    counts = Counter((r["feed"],r["loading_pct"]) for r in rows)
    weights = (np.array([len(rows)/len(counts)/counts[r["feed"],r["loading_pct"]] for r in rows])
               if balanced else np.ones(len(rows)))
    fixed_K = []
    if variant == "published_K_training_anchors":
        for feed in sorted({r["feed"] for r in rows}):
            if (feed,10) not in counts:
                raise ValueError("The published Kx is from the 10% curve; that curve must be in training.")
            fixed_K.append({1:11.15,2:4.53}[feed])
    result = model.fit_shared(law, array(rows,"nominal_concentration_g_L"),
        array(rows,"exposure_kg_s_mol_acid"), array(rows,ratio), array(rows,response),
        weights,np.array([r["feed"] for r in rows]),reference_concentration=ref,
        fixed_shape=float("nan") if fixed_shape is None else fixed_shape,profile_points=points,
        fixed_equilibria=np.array(fixed_K,dtype=float))
    params = {}
    for group,fit in zip(result.series,result.fits,strict=True):
        gradient = float(fit.gradient_inf)
        p = {"k_ref_mol_kg_s":float(fit.k_ref),"apparent_Kx":float(fit.Kx),
             "solver":str(fit.solver),"gradient_inf":gradient if np.isfinite(gradient) else None,
             "equilibrium_limit_probe_conversion_change":float(fit.limit_probe_change)}
        ordinary = (p["solver"]=="least_squares" and .001 < p["apparent_Kx"] < 1000
                    and np.isfinite(gradient) and gradient <= 1e-8)
        limit = (p["solver"]=="upper_equilibrium_limit" and p["apparent_Kx"] > 1e8
                 and p["equilibrium_limit_probe_conversion_change"] <= 1e-5)
        if not (.0001 < p["k_ref_mol_kg_s"] < 10 and (ordinary or limit)):
            raise RuntimeError(f"Unresolved selected nuisance fit: {law}, {variant}, {p}")
        params[str(int(group))] = p
    return {"law":law,"shape":float(result.shape),"weighted_sse":float(result.weighted_sse),
            "parameters_by_feed_ratio":params}


def shared_predict(model, fit, rows, ref, variant):
    response,ratio,_ = VARIANTS[variant]
    records = []
    for group in sorted({r["feed"] for r in rows}):
        selected = [r for r in rows if r["feed"] == group]
        p = fit["parameters_by_feed_ratio"][str(group)]
        x = np.asarray(model.predict_series(fit["law"],fit["shape"],p["k_ref_mol_kg_s"],
            p["apparent_Kx"],array(selected,"nominal_concentration_g_L"),
            array(selected,"exposure_kg_s_mol_acid"),array(selected,ratio),
            reference_concentration=ref))
        records.extend({"feed":group,"loading_pct":r["loading_pct"],"assay":r["assay"],
            "time_min":r["time_min"],"observed_conversion":r[response],"predicted_conversion":float(y)}
            for r,y in zip(selected,x,strict=True))
    return records


def shared_scores(records):
    return error_metrics(array(records,"predicted_conversion")-array(records,"observed_conversion"))


def combined_batch_fits(model,plan,groups,ref):
    keys = sorted(groups)
    summary,fit_rows,pred_rows = {},[],[]
    max_refinement = 0.
    def case(name,variant,training,testing,law,fixed_shape=None,refine=False):
        nonlocal max_refinement
        train,test = select_rows(groups,training),select_rows(groups,testing)
        fit = shared_fit(model,law,train,ref,variant,fixed_shape=fixed_shape)
        pred = shared_predict(model,fit,test,ref,variant)
        if refine and law != "proportional":
            fine = shared_fit(model,law,train,ref,variant,fixed_shape=fixed_shape,points=49)
            comparison = shared_predict(model,fine,test,ref,variant)
            change = float(np.max(abs(array(pred,"predicted_conversion")-array(comparison,"predicted_conversion"))))
            max_refinement = max(max_refinement,change)
            if change > plan["models"]["maximum_primary_refinement_absolute_conversion"]:
                raise RuntimeError("The shared concentration profile failed refinement.")
        for g,p in fit["parameters_by_feed_ratio"].items():
            fit_rows.append({"case":name,"variant":variant,"law":law,"shape":fit["shape"],
                "training_assays":"/".join(str(groups[k][0]["assay"]) for k in training),
                "feed_ratio":g,"weighted_sse":fit["weighted_sse"],**p})
        pred_rows.extend(r | {"case":name,"variant":variant,"law":law} for r in pred)
        return fit | {"scores":shared_scores(pred),"scores_by_assay": {
            str(a):shared_scores([r for r in pred if r["assay"]==a]) for a in sorted({r["assay"] for r in pred})}}
    for family,selected in (("2to1",[k for k in keys if k[0]==2]),
                            ("1to1",[k for k in keys if k[0]==1]),("joint",keys)):
        full,holdout = {},{}
        for law in LAWS:
            full[law] = case(f"{family}_all","primary",selected,selected,law,refine=family=="joint")
        for key in selected:
            assay = str(groups[key][0]["assay"])
            train = [k for k in selected if k!=key]
            holdout[assay] = {law:case(f"{family}_holdout_{assay}","primary",train,[key],law,
                                       refine=family=="joint" and key==(1,20)) for law in LAWS}
        score = {law:float(np.sqrt(np.mean([v[law]["scores"]["rmse_percentage_points"]**2
                                         for v in holdout.values()]))) for law in LAWS}
        summary[family] = {"full":full,"holdouts_by_assay":holdout,"equal_assay_rmse_pp":score}
    summary["ratio_transfer"] = {}
    for r,other in ((1,2),(2,1)):
        label = f"{other}to1_shape_to_{r}to1"
        summary["ratio_transfer"][label] = {
            law:case(label,"primary",[(r,10)],[k for k in keys if k[0]==r and k[1]!=10],law,
                     fixed_shape=summary[f"{other}to1"]["full"][law]["shape"]) for law in LAWS}
    high = (1,20)
    summary["sensitivity"] = {}
    for variant in list(VARIANTS)[1:]:
        summary["sensitivity"][variant] = {
            "full": {law:case("joint_all",variant,keys,keys,law) for law in LAWS},
            "highest_1to1_held_out":{law:case("joint_holdout_8",variant,[k for k in keys if k!=high],
                                                [high],law) for law in LAWS},
        }
    joint = summary["joint"]
    summary["qualification"] = {}
    for law in LAWS[1:]:
        improvement = 1-joint["equal_assay_rmse_pp"][law]/joint["equal_assay_rmse_pp"]["proportional"]
        worsening = max(v[law]["scores"]["rmse_percentage_points"]-v["proportional"]["scores"]["rmse_percentage_points"]
                        for v in joint["holdouts_by_assay"].values())
        high_tests = joint["holdouts_by_assay"]["8"]
        high_improvement = high_tests[law]["scores"]["rmse_percentage_points"] < high_tests["proportional"]["scores"]["rmse_percentage_points"]
        transfers = all(v[law]["scores"]["rmse_percentage_points"] < v["proportional"]["scores"]["rmse_percentage_points"]
                        for v in summary["ratio_transfer"].values())
        sensitivity = all(v["highest_1to1_held_out"][law]["scores"]["rmse_percentage_points"] <
                          v["highest_1to1_held_out"]["proportional"]["scores"]["rmse_percentage_points"]
                          for v in summary["sensitivity"].values())
        identified = all(all(p["solver"]=="least_squares" for p in v[law]["parameters_by_feed_ratio"].values())
                         for v in joint["holdouts_by_assay"].values())
        summary["qualification"][law] = {"joint_relative_rmse_improvement":improvement,
            "largest_joint_fold_worsening_pp":worsening,"improves_highest_1to1":high_improvement,
            "improves_both_ratio_transfers":transfers,"highest_improvement_in_all_sensitivities":sensitivity,
            "alternative_equilibrium_identified_in_every_joint_fold":identified,
            "passes":bool(improvement>=.1 and worsening<=.2 and high_improvement and transfers and sensitivity and identified)}
    summary["maximum_refinement_conversion_change"] = max_refinement
    return summary,fit_rows,pred_rows


def transfer_scenarios(jl, groups, basis, ref, batch, external):
    source = read_json(prior.published.MODEL_INPUTS)
    bed = source["bed"]
    volume = pi*bed["diameter_m"]**2/4*bed["length_m"]
    mass = bed["catalyst_mass_kg"]
    if mass != .006:
        raise ValueError("The actual catalyst inventory remains 6 g.")
    concentrations = {"void":mass/(volume*bed["bed_porosity"]),"total_volume_proxy":mass/volume}
    scale = lambda law,shape,relative: float(jl.CatalystLoading.rate_scale(law,shape,relative))
    factors = {"mass_proportional":1.}
    for law in LAWS[1:]:
        shape = batch["joint"]["full"][law]["shape"]
        factors[f"joint_pentyl_{law}_void"] = scale(law,shape,concentrations["void"]/ref)
    factors["joint_pentyl_crowding_total_volume_proxy"] = scale("crowding",
        batch["joint"]["full"]["crowding"]["shape"],concentrations["total_volume_proxy"]/ref)
    # Independent external law is on g/100g liquid reactants, not g/L.
    # The nominal 2:1 starting charge density is only a comparison proxy
    # for the unmeasured high-solids contact environment.
    w0 = groups[(2,10)][0]["loading_pct_total_substrates"]
    nominal_density = basis[2]["nominal_liquid_density_g_L"]
    w_bed = 100*concentrations["void"]/nominal_density
    for law in LAWS[1:]:
        shape = external["common_external"][law]["shape"]
        factors[f"external_only_{law}_mass_ratio_proxy"] = (
            scale(law,shape,w_bed/10)/scale(law,shape,w0/10))
    return factors, {"catalyst_mass_kg":mass,"concentrations_g_L":concentrations,
        "nominal_batch_reference_g_L":ref,"nominal_2to1_charge_density_g_L":nominal_density,
        "reference_pct_of_nominal_substrates":w0,"bed_pct_of_nominal_liquid_mass_proxy":w_bed,
        "factor_scenarios":factors,"scope":"Diagnostic extrapolations only; no accepted new packed-bed law."}


def run(output):
    plan = read_json(PLAN)
    groups,basis,ref = pentyl_groups(plan)
    jl = runtime()
    for name,path in (("PentylPropanoate",base.HERE/"model.jl"),("CatalystLoading",prior.HERE/"model.jl"),
            ("SharedLoadingEvidence",HERE/"shared_loading.jl"),
            ("LorenzoPackedBed",prior.published.MODEL),("BeadFilm",prior.film.HERE/"model.jl"),
            ("LoadingFilmTransfer",prior.HERE/"film_transfer.jl")):
        if not jl.seval(f"isdefined(@__MODULE__, :{name})"):
            jl.include(str(path))
    model = jl.SharedLoadingEvidence
    external,internal = all_progress(plan,groups)
    external_summary,folds,progress_pred = compare_progress_laws(model,plan,external,internal)
    alternatives = {}
    for fig in (3,5):
        alternate,_ = all_progress(plan,groups,paroul_figure=fig)
        alt,_,_ = compare_progress_laws(model,plan,alternate,internal)
        alternatives[str(fig)] = {"paroul_high_loading_matches":[r for r in alternate["paroul2010"]
            if r["loading_pct_reported"]==10],"common_external":alt["common_external"]}
    reading_sensitivities = {}
    for name,direction in (("more_high_loading_effective",True),("less_high_loading_effective",False)):
        alt,_,_ = compare_progress_laws(model,plan,reading_progress_scenario(external,direction),internal)
        reading_sensitivities[name] = alt["common_external"]
    batch,fit_rows,batch_pred = combined_batch_fits(model,plan,groups,ref)
    factors,transfer = transfer_scenarios(jl,groups,basis,ref,batch,external_summary)
    old_plan = read_json(prior.PLAN)
    anchors,anchor_records = prior.anchor_batch_kinetics(jl.PentylPropanoate,read_rows(base.BATCH),old_plan)
    pe = read_json(prior.TRACER)["conditional_danckwerts_transfer"]["nominal_rectangular_pulse"]["danckwerts_peclet"]
    bed,numerics = prior.predict_bed(jl,read_json(prior.published.MODEL_INPUTS),
        read_json(prior.film.PROPERTIES),factors,anchors,old_plan,pe)
    output.mkdir(parents=True,exist_ok=False)
    for name,rows in (("matched_progress.csv",[r for values in [*external.values(),*internal.values()] for r in values]),
        ("external_fold_scores.csv",folds),("external_progress_predictions.csv",progress_pred),
        ("batch_fits.csv",fit_rows),("batch_predictions.csv",batch_pred),("predictions_before_comparison.csv",bed)):
        write_rows(output/name,rows)
    summary = {
        "batch_groups":{f"{ratio}to1_{w}":{k:rows[0][k] for k in (
            "assay","mass_g","nominal_concentration_g_L","loading_pct_total_substrates","printed_ratio")}
            | {"rows_including_start":len(rows)} for (ratio,w),rows in groups.items()},
        "charge_basis":{str(k):v for k,v in basis.items()},
        "direct_observed_exposure_checks":direct_exposure_checks(plan),
        "external_studies":external_summary,"paroul_same_condition_sensitivity":alternatives,
        "external_marker_reading_sensitivity":reading_sensitivities,
        "equimolar_printing_consistency": {
            "max_printed_vs_acid_rederived_conversion": max(abs(r["printed"]-r["from_acid"])
                for k,rows in groups.items() if k[0]==1 for r in rows),
            "max_printed_vs_ester_rederived_conversion": max(abs(r["printed"]-r["from_ester"])
                for k,rows in groups.items() if k[0]==1 for r in rows),
            "max_mole_fraction_sum_error":max(abs(sum(float(r[k]) for k in
                ("x_acid","x_alcohol","x_ester","x_water"))-1) for r in read_rows(BATCH_1))
        },
        "pentyl_batch":batch,"transfer":transfer,"ten_percent_2to1_batch_anchors":anchor_records,
        "bed_numerics":numerics,"julia":str(jl.seval("VERSION")),
    }
    # Freeze all independent shapes and actual reactor predictions before
    # the measured reactor conversions are read.
    write_json(output/"frozen_calibration.json",summary)
    data,metrics = prior.compare_bed(bed,factors)
    write_rows(output/"bed_comparison.csv",data)
    summary["bed_comparison"] = metrics
    summary["sixty_C_three_ml_min"] = next(r for r in data if prior.published.condition(r)==(60,3))
    sources = [PLAN,HERE/"shared_loading.jl",HERE/"run.py",BATCH_1,base.BATCH,base.EQUILIBRIUM,
        base.HERE/"model.jl",base.HERE/"run.py",prior.PLAN,prior.HERE/"run.py",prior.HERE/"model.jl",
        prior.HERE/"film_transfer.jl",prior.published.MODEL,prior.published.MODEL_INPUTS,prior.published.HOLDOUT,
        prior.film.PROPERTIES,prior.film.HERE/"model.jl",prior.FILM_REFERENCE,prior.TRACER]
    for series in plan["external"]["series"]:
        sources.extend([SOURCE/series/"loading_progress.csv",SOURCE/series/"metadata.json",SOURCE/series/"README.md"])
    for series in ("richetti2010","paroul2010"):
        sources.append(SOURCE/series/"initial_rates.csv")
    sources.extend([SOURCE/"meneses2019/metadata.json",SOURCE/"meneses2019/README.md"])
    summary["source_sha256"] = {str(path.relative_to(ROOT)):sha(path) for path in sources}
    summary["output_sha256"] = {path.name:sha(path) for path in sorted(output.iterdir())}
    write_json(output/"summary.json",summary)
    return summary


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    parser.add_argument("--output",type=Path,default=OUTPUT/stamp)
    args = parser.parse_args()
    report = run(args.output.resolve())
    print("Output:",args.output)
    print("Joint batch held-out RMSE (pp):",report["pentyl_batch"]["joint"]["equal_assay_rmse_pp"])
    print("Same-chemistry screen:",report["pentyl_batch"]["qualification"])
    for law,score in report["external_studies"]["common_external"].items():
        print("External",law, score["shape"],"held-paper log-time RMSE",score["equal_paper_log_time_rmse"])
    for name,score in report["bed_comparison"].items():
        print(name, "held-out reactor RMSE",f'{score["rmse_percentage_points"]:.3f} pp')
