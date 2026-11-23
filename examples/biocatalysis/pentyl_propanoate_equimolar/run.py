"""Compare both feeds in the same six-gram bed with frozen batch-only models."""

from argparse import ArgumentParser
from datetime import datetime, timezone
from itertools import product
from math import pi
from pathlib import Path
from time import perf_counter
import platform

import numpy as np

from examples.biocatalysis.pentyl_propanoate.run import array, error_metrics, read_rows, write_rows
from examples.biocatalysis.pentyl_propanoate_composition import run as prior
from examples.biocatalysis.pentyl_propanoate_composition.check import check as check_prior
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT

HERE = Path(__file__).resolve().parent
PLAN = HERE / "protocol.json"
SOURCE = ROOT / "literature/sources/lorenzo2025"
EXTRA = SOURCE / "equimolar_extraction.json"
OBSERVATIONS = {1: SOURCE / "conversion_1to1.csv", 2: SOURCE / "conversion_2to1.csv"}
FROZEN = prior.HERE / "reference/summary.json"
PRIOR_PREDICTIONS = prior.HERE / "reference/predictions_before_comparison.csv"
PRIOR_CLOCK = prior.HERE / "reference/stationary_clock_screen.csv"
PUBLISHED = prior.loading.published
PROPERTIES = prior.loading.film.PROPERTIES
OUTPUT = ROOT / "results/biocatalysis/pentyl_propanoate_equimolar"
REFERENCE_CASE = "prior_central_film_reference"
ME_G_MOL = 74.0785 + 88.148 - 18.01528


def load_runtime():
    jl = prior.load_runtime()
    if not jl.seval("isdefined(@__MODULE__, :DualFeedBed)"):
        jl.include(str(HERE / "model.jl"))
    return jl


def case_specs(frozen):
    """Use only the saved training records, not the reactor comparison scores."""
    b = frozen["batch"]
    out = {REFERENCE_CASE: {"training": "published_Table_S3", "activity": 1.}}
    for label, train, laws in (("all_batch", "all_data", prior.LAWS),
                               ("active_batch", "hot_acid_challenge", prior.STATIC)):
        f = b[train]
        for law in laws:
            name = f"{label}_{law}" + ("_fully_active" if law in prior.CLOCK else "")
            out[name] = {"training": train, "law": law, "activity": 1.,
                         "training_assays": f["training_assays"],
                         "parameters": f["fit"][law]["parameters"],
                         "equilibrium_training": f["equilibrium_training"]}
    for law in ("residual_decay", "protected_residual_decay"):
        f = b["all_data"]
        out[f"all_batch_{law}_long_time_floor"] = {
            "training": "all_data", "law": law,
            "activity": f["fit"][law]["parameters"]["residual_active_fraction"],
            "training_assays": f["training_assays"],
            "parameters": f["fit"][law]["parameters"],
            "equilibrium_training": f["equilibrium_training"]}
    if len(out) != 22:
        raise ValueError("Expected all 22 frozen reactor scenarios.")
    return out


def conditions(jl, source, extra, liquid):
    bed, equimolar = source["bed"], extra["reaction_1to1"]
    if not (bed["catalyst_mass_kg"] == equimolar["catalyst_mass_kg"] == .006
            and bed["diameter_m"] == equimolar["column_bore_m"]
            and bed["length_m"] == equimolar["column_length_m"]):
        raise ValueError("Both feeds must use the same measured six-gram bed.")
    area = pi * bed["diameter_m"]**2 / 4
    void, diameter = bed["bed_porosity"], bed["bead_mean_diameter_m"]
    volume_L = area * bed["length_m"] * void * 1000
    catalyst_g_L = 1000 * bed["catalyst_mass_kg"] / volume_L
    a_l = 6 * (1-void) / (diameter * void)
    pe = bed["length_m"] * .20 / (void**2 * diameter)
    by_T = {round(r["temperature_K"]-273.15): r for r in liquid["conditions"]}
    rows = []
    for feed, T, Q in product((1, 2), bed["reaction_temperature_C"], bed["feed_flow_ml_per_min"]):
        p = by_T[T]
        if feed == 1:
            acid, alcohol = equimolar["inlet_acid_mol_L"], equimolar["inlet_alcohol_mol_L"]
            rho = acid * 74.0785 + alcohol * 88.148  # g/L = kg/m³, calculated nominal dry density
            eta = extra["effectiveness"]["by_temperature_C"][str(T)]["eta_internal"]
        else:
            acid, alcohol = bed["feed_acid_mol_per_L"], bed["feed_alcohol_mol_per_L"]
            rho = p["density_kg_m3"]
            eta = source["effectiveness_by_temperature"][str(T)]["eta_internal"]
        if alcohol != feed * acid:
            raise ValueError("The feed ratio does not match the batch model selection.")
        c = {"feed_ratio": feed, "temperature_C": T, "flow_ml_min": Q,
             "actual_mass_kg": bed["catalyst_mass_kg"], "catalyst_g_L": catalyst_g_L,
             "acid0_mol_L": acid, "alcohol0_mol_L": alcohol, "void_volume_L": volume_L,
             "residence_min": volume_L * 1000 / Q, "eta_internal": eta, "peclet": pe,
             "superficial_velocity_m_s": Q / 60e6 / area,
             "bead_diameter_m": diameter, "bed_void_fraction": void, "a_l_m_inv": a_l,
             "density_proxy_kg_m3": rho, "viscosity_proxy_Pa_s": p["viscosity_Pa_s"],
             "diffusivity_proxy_m2_s": p["diffusivity_m2_s"]["central"]}
        rows.append(transport(jl, c, "central"))
    if len(rows) != 40:
        raise ValueError("Expected twenty matched conditions for each feed.")
    return rows


def transport(jl, c, variant):
    multiplier = {"central": 1., "half_diffusivity": .5, "double_diffusivity": 2.,
                  "no_liquid_film": 1., "no_additional_internal_factor": 1.}[variant]
    D = c["diffusivity_proxy_m2_s"] * multiplier
    f = jl.BeadFilm.film_coefficients(c["superficial_velocity_m_s"], c["bead_diameter_m"],
                                    c["density_proxy_kg_m3"], c["viscosity_proxy_Pa_s"], D,
                                    c["bed_void_fraction"], "wilson_geankoplis")
    eta = 1. if variant == "no_additional_internal_factor" else c["eta_internal"]
    return c | {"transport_case": variant, "diffusivity_used_m2_s": D, "eta_internal": eta,
                "use_film": variant != "no_liquid_film", "reynolds": float(f.reynolds),
                "schmidt": float(f.schmidt), "sherwood": float(f.sherwood),
                "film_coefficient_m_s": float(f.coefficient_m_s),
                "capacity_mol_L_min": 60 * float(f.coefficient_m_s) * c["a_l_m_inv"] * c["acid0_mol_L"]}


def frozen_model(jl, spec, c, published_model, *, activity=None):
    if spec["training"] == "published_Table_S3":
        return published_model, None
    kinetics = jl.CompositionKinetics
    p = prior.parameters_from_record(kinetics, spec)
    K = prior.equilibrium_value(kinetics, spec["equilibrium_training"], c["feed_ratio"], c["temperature_C"])
    a = spec["activity"] if activity is None else activity
    return jl.DualFeedBed.FrozenKinetics(p, K, a, c["feed_ratio"]), K


def ester_g_h(conversion, c):
    # Inferred outlet production at constant Vdot, dry feed and one reaction.
    return conversion * c["acid0_mol_L"] * c["flow_ml_min"] / 1000 * 60 * ME_G_MOL


def solve(jl, kinetic, c, settings):
    args = (kinetic, c["temperature_C"], c["residence_min"], c["catalyst_g_L"],
            c["capacity_mol_L_min"], c["eta_internal"], c["peclet"])
    kwargs = {"acid0": c["acid0_mol_L"], "alcohol0": c["alcohol0_mol_L"],
              "use_film": c["use_film"]}
    coarse, fine = [jl.DualFeedBed.prediction(*args, **kwargs, cells=n) for n in settings["cells"]]
    grid = float(abs(coarse.conversion-fine.conversion))
    fields = ("balance_error", "max_film_balance_error_mol_L_min", "max_positive_rate_slope_over_film_capacity",
              "film_slope_sampling_refinement_change")
    balance, film, slope, slope_change = [max(float(getattr(coarse, k)), float(getattr(fine, k))) for k in fields]
    x, eq = float(fine.conversion), float(fine.equilibrium_conversion)
    nums = (grid, balance, film, slope, slope_change, x, eq)
    if not (np.all(np.isfinite(nums)) and 0 <= x <= eq + 1e-12 <= 1+1e-12
            and grid <= settings["max_absolute_grid_difference"]
            and balance <= settings["max_normalized_reaction_flux_balance_error"]
            and film <= settings["max_film_root_error_mol_L_min"]
            and slope < settings["max_sampled_positive_slope_over_film_capacity"]):
        raise RuntimeError(f"A numerical acceptance check failed: {c}; grid={grid}, balance={balance}, film={film}, slope={slope}, X={x}, Xeq={eq}")
    return {"predicted_conversion": x, "equilibrium_conversion": eq,
            "nominal_outlet_ester_g_h": ester_g_h(x,c), "grid_change_absolute_conversion": grid,
            "normalized_balance_error": balance, "film_balance_error_mol_L_min": film,
            "max_positive_rate_slope_over_film_capacity": slope,
            "film_slope_sampling_refinement_change": slope_change,
            "inlet_backdiffusion_fraction": float(fine.inlet_backdiffusion_fraction),
            "integrated_source_fraction": float(fine.integrated_source_fraction),
            "surface_progress_at_inlet": float(fine.surface_progress_at_inlet)}


def identity(c):
    return int(c["feed_ratio"]), int(c["temperature_C"]), float(c["flow_ml_min"])


def row_head(c):
    return {k: c[k] for k in ("feed_ratio", "temperature_C", "flow_ml_min", "actual_mass_kg",
            "catalyst_g_L", "acid0_mol_L", "alcohol0_mol_L", "residence_min", "eta_internal",
            "peclet", "transport_case", "diffusivity_used_m2_s", "capacity_mol_L_min", "use_film")}


def predict_campaign(jl, frozen, source, cases, plan):
    specs = case_specs(frozen)
    numerical = plan["numerical_acceptance"]
    p = jl.LorenzoPackedBed.PingPong(np.array([source["rate_parameters"][k] for k in PUBLISHED.PARAMETERS]))
    central, clock, sensitivity = [], [], []
    for name, spec in specs.items():
        print(f"Frozen dual-feed: {name}", flush=True)
        for c in cases:
            model, K = frozen_model(jl,spec,c,p)
            central.append({"case": name, **row_head(c), "apparent_Kx": K,
                            "activity": spec["activity"], **solve(jl,model,c,numerical)})
    for law in plan["age"]["laws"]:
        name = f"all_batch_{law}_fully_active"
        spec = specs[name]
        params = prior.parameters_from_record(jl.CompositionKinetics, spec)
        print(f"Stationary age endpoints: {law}", flush=True)
        for c in cases:
            feed,T,_ = identity(c)
            K = prior.equilibrium_value(jl.CompositionKinetics,spec["equilibrium_training"],feed,T)
            state = jl.CompositionKinetics.static_terms(params,T,float(feed),K,feed=feed)
            for age in plan["age"]["ages_min"]:
                envelope = jl.CompositionKinetics.activity_bounds(state,float(age))
                lo, hi = float(envelope.lower), float(envelope.upper)
                if not 0 <= lo <= hi <= 1:
                    raise RuntimeError("Nonphysical stationary activity interval.")
                make = lambda a: jl.DualFeedBed.FrozenKinetics(params,K,a,feed)
                lower = solve(jl,make(lo),c,numerical)
                upper = lower if lo == hi else solve(jl,make(hi),c,numerical)
                if lower["predicted_conversion"] > upper["predicted_conversion"] + 1e-8:
                    raise RuntimeError("Nonmonotone uniform-activity endpoints.")
                clock.append({"law": law, **row_head(c), "catalyst_age_min": age, "apparent_Kx": K,
                              "inlet_activity": lo, "equilibrium_activity": hi,
                              "inlet_hazard_min_inv": float(envelope.inlet_hazard_min_inv),
                              "equilibrium_hazard_min_inv": float(envelope.equilibrium_hazard_min_inv),
                              "recovery_fraction": float(envelope.recovery_fraction),
                              "conversion_uniform_inlet_activity": lower["predicted_conversion"],
                              "conversion_uniform_equilibrium_activity": upper["predicted_conversion"],
                              **{key: max(lower[key],upper[key]) for key in (
                                  "grid_change_absolute_conversion", "normalized_balance_error",
                                  "film_balance_error_mol_L_min", "max_positive_rate_slope_over_film_capacity",
                                  "film_slope_sampling_refinement_change")}})
    for name in plan["transport_sensitivities"]["cases"]:
        spec = specs[name]
        print(f"Transport-only sensitivity: {name}", flush=True)
        for c0, scenario in product(cases, plan["transport_sensitivities"]["scenarios"]):
            c = transport(jl,c0,scenario)
            model,K = frozen_model(jl,spec,c,p)
            sensitivity.append({"case": name, **row_head(c), "apparent_Kx": K,
                                "activity": spec["activity"], **solve(jl,model,c,numerical)})
    return central,clock,sensitivity


def compare_previous(central,clock,plan):
    """All existing 2:1 predictions must replay, including clock endpoints."""
    central_old = {(r["case"],PUBLISHED.condition(r)):float(r["predicted_conversion"])
                   for r in read_rows(PRIOR_PREDICTIONS)}
    central_new = {(r["case"],identity(r)[1:]):r["predicted_conversion"]
                   for r in central if r["feed_ratio"] == 2}
    key = lambda r: (r["law"],PUBLISHED.condition(r),float(r["catalyst_age_min"]))
    old_clock = {key(r):r for r in read_rows(PRIOR_CLOCK)}
    new_clock = {key(r):r for r in clock if r["feed_ratio"] == 2}
    if not (set(central_old) == set(central_new) and set(old_clock) == set(new_clock)
            and len(central_new) == 440 and len(new_clock) == 280):
        raise ValueError("The frozen 2:1 prediction ledger changed.")
    max_c = max(abs(central_new[k]-v) for k,v in central_old.items())
    max_a = max(abs(float(new_clock[k][field])-float(v[field])) for k,v in old_clock.items()
                for field in ("conversion_uniform_inlet_activity","conversion_uniform_equilibrium_activity"))
    if max(max_c,max_a) > plan["numerical_acceptance"]["prior_2to1_replay_tolerance_absolute_conversion"]:
        raise RuntimeError(f"Previous 2:1 predictions did not replay: {max_c}, {max_a}")
    return {"central_count":len(central_new),"clock_condition_count":len(new_clock),
            "max_central_absolute_difference":max_c,"max_clock_endpoint_absolute_difference":max_a}


def read_observations():
    selected = {}
    for feed,path in OBSERVATIONS.items():
        for row in read_rows(path):
            T,Q = PUBLISHED.condition(row)
            key = feed,T,Q
            conversion = float(row["acid_conversion"])
            if key in selected or not 0 <= conversion <= 1:
                raise ValueError("Duplicate or nonphysical reactor observation.")
            selected[key] = {"source_run_id":int(row["source_run_id"]),"observed_conversion":conversion}
        expected = set(range(45 if feed == 1 else 65, 65 if feed == 1 else 85))
        got = {v["source_run_id"] for k,v in selected.items() if k[0] == feed}
        if got != expected:
            raise ValueError(f"Incomplete held-out reactor cases for feed {feed}.")
    if set(selected) != set(product((1,2),(40,50,60,70,80),(.7,1.,2.,3.))):
        raise ValueError("The new comparison must contain the complete matched condition matrix.")
    return selected


def scored(data):
    return error_metrics(array(data,"predicted_conversion")-array(data,"observed_conversion"))


def score_groups(rows):
    answer = {"joint_40":scored(rows)}
    for feed in (1,2):
        scoped = [r for r in rows if r["feed_ratio"] == feed]
        answer[f"{feed}to1"] = scored(scoped)
        answer[f"{feed}to1_40to60C"] = scored([r for r in scoped if r["temperature_C"] <= 60])
        answer[f"{feed}to1_70to80C"] = scored([r for r in scoped if r["temperature_C"] >= 70])
    if answer["joint_40"]["count"] != 40 or answer["1to1"]["count"] != 20 or answer["2to1"]["count"] != 20:
        raise ValueError("Missing conditions in a model comparison.")
    return answer


def add_observation(row, observed):
    o = observed[identity(row)]
    return row | o | {"observed_nominal_outlet_ester_g_h":ester_g_h(o["observed_conversion"],row),
                      "conversion_error_pp":100*(row["predicted_conversion"]-o["observed_conversion"])}


def paired_comparison(central_data, observations):
    out, paired = {}, []
    for name in sorted({r["case"] for r in central_data}):
        rows = [r for r in central_data if r["case"] == name]
        by_key = {identity(r):r for r in rows}
        if set(by_key) != set(observations) or len(rows) != 40:
            raise ValueError(f"Incomplete predictions for {name}.")
        for T,Q in product((40,50,60,70,80),(.7,1.,2.,3.)):
            a,b = by_key[1,T,Q],by_key[2,T,Q]
            paired.append({"case":name,"temperature_C":T,"flow_ml_min":Q,
                           "observed_2to1_minus_1to1_pp":100*(b["observed_conversion"]-a["observed_conversion"]),
                           "predicted_2to1_minus_1to1_pp":100*(b["predicted_conversion"]-a["predicted_conversion"]),
                           "observed_2to1_minus_1to1_ester_g_h":b["observed_nominal_outlet_ester_g_h"]-a["observed_nominal_outlet_ester_g_h"],
                           "predicted_2to1_minus_1to1_ester_g_h":b["nominal_outlet_ester_g_h"]-a["nominal_outlet_ester_g_h"]})
        p = [r for r in paired if r["case"] == name]
        hot = {f"{feed}to1": {f"{Q:g}": {
            "observed_80_minus_60_pp":100*(by_key[feed,80,Q]["observed_conversion"]-by_key[feed,60,Q]["observed_conversion"]),
            "predicted_80_minus_60_pp":100*(by_key[feed,80,Q]["predicted_conversion"]-by_key[feed,60,Q]["predicted_conversion"])}
            for Q in (.7,1.,2.,3.)} for feed in (1,2)}
        out[name] = score_groups(rows) | {"paired_feed_difference":error_metrics(
            (array(p,"predicted_2to1_minus_1to1_pp")-array(p,"observed_2to1_minus_1to1_pp"))/100),
            "correct_feed_advantage_count":int(np.count_nonzero(array(p,"predicted_2to1_minus_1to1_pp")>0)),
            "hot_temperature_change":hot}
    return out, paired


def classify_clock(clock, observed):
    data = [r | observed[identity(r)] for r in clock]
    groups = {}
    for law,feed,age in product(prior.CLOCK,(1,2),(30,180)):
        sample = [r for r in data if r["law"] == law and r["feed_ratio"] == feed and r["catalyst_age_min"] == age]
        if len(sample) != 20:
            raise ValueError("Incomplete clock envelope.")
        def summarize(rows):
            lo,hi,y = [array(rows,key) for key in (
                "conversion_uniform_inlet_activity","conversion_uniform_equilibrium_activity","observed_conversion")]
            return {"count":len(rows),"covered_count":int(np.count_nonzero((lo <= y) & (y <= hi))),
                    "covered_with_0p5pp_reading_allowance":int(np.count_nonzero((lo-.005 <= y) & (y <= hi+.005))),
                    "observed_above_envelope_count":int(np.count_nonzero(y > hi)),
                    "observed_below_envelope_count":int(np.count_nonzero(y < lo)),
                    "mean_envelope_width_pp":float(100*np.mean(hi-lo)),
                    "max_distance_outside_envelope_pp":float(100*np.max(np.maximum.reduce((lo-y,y-hi,np.zeros_like(y)))))}
        groups[f"{law}_{feed}to1_{age}min"] = summarize(sample) | {
            "40to60C":summarize([r for r in sample if r["temperature_C"] <= 60]),
            "70to80C":summarize([r for r in sample if r["temperature_C"] >= 70])}
    return data,groups


def score_transport(rows, central, observed):
    lookup = {(r["case"],identity(r)):r["predicted_conversion"] for r in central}
    results = [add_observation(r,observed) | {
        "change_from_central_pp":100*(r["predicted_conversion"]-lookup[r["case"],identity(r)])} for r in rows]
    summary = {}
    for name,scenario in sorted({(r["case"],r["transport_case"]) for r in rows}):
        r = [r for r in results if (r["case"],r["transport_case"]) == (name,scenario)]
        summary[name+"__"+scenario] = score_groups(r) | {f"{feed}to1_change_pp":{
            "min":min(d["change_from_central_pp"] for d in r if d["feed_ratio"] == feed),
            "max":max(d["change_from_central_pp"] for d in r if d["feed_ratio"] == feed),
            "max_abs":max(abs(d["change_from_central_pp"]) for d in r if d["feed_ratio"] == feed)} for feed in (1,2)}
    return results,summary


def run(output):
    check_prior()
    plan,source,extra,liquid,frozen = [read_json(p) for p in (PLAN,PUBLISHED.MODEL_INPUTS,EXTRA,PROPERTIES,FROZEN)]
    jl = load_runtime()
    cases = conditions(jl,source,extra,liquid)
    output.mkdir(parents=True,exist_ok=False)
    central,clock,sensitivity = predict_campaign(jl,frozen,source,cases,plan)
    if (len(central),len(clock),len(sensitivity)) != (880,560,800):
        raise ValueError("Incomplete frozen campaign.")
    # Persist prediction ledgers before reading the held-out conversions.
    write_rows(output/"predictions_before_comparison.csv",central)
    write_rows(output/"clock_before_comparison.csv",clock)
    write_rows(output/"transport_before_comparison.csv",sensitivity)
    replay = compare_previous(central,clock,plan)
    observed = read_observations()
    data = [add_observation(r,observed) for r in central]
    comparison,paired = paired_comparison(data,observed)
    clock_data,clock_scores = classify_clock(clock,observed)
    transport_data,transport_scores = score_transport(sensitivity,central,observed)
    for name,rows in (("bed_comparison.csv",data),("paired_feed_comparison.csv",paired),
                      ("clock_comparison.csv",clock_data),("transport_comparison.csv",transport_data)):
        write_rows(output/name,rows)
    maxima = {key:max(r[key] for r in central+clock+sensitivity) for key in (
        "grid_change_absolute_conversion","normalized_balance_error","film_balance_error_mol_L_min",
        "max_positive_rate_slope_over_film_capacity","film_slope_sampling_refinement_change")}
    files = [PLAN,HERE/"run.py",HERE/"model.jl",EXTRA,*OBSERVATIONS.values(),FROZEN,
             PRIOR_PREDICTIONS,PRIOR_CLOCK,PUBLISHED.MODEL_INPUTS,PUBLISHED.MODEL,PROPERTIES,
             prior.HERE/"model.jl",prior.HERE/"run.py",prior.HERE/"bed_transfer.jl",
             prior.loading.HERE/"model.jl",prior.base.HERE/"model.jl",prior.loading.film.HERE/"model.jl"]
    report = {"model_boundary":"Zero new fits; previous batch parameters fixed; all 40 reactor measurements retained.",
              "python":platform.python_version(),"julia":str(jl.seval("VERSION")),"platform":platform.platform(),
              "catalyst_mass_kg":.006,"condition_count":len(cases),"central_count":len(central),
              "clock_count":len(clock),"transport_count":len(sensitivity),"reactor_conditions":cases,
              "prior_2to1_replay":replay,"numerical_maxima":maxima,"central_scores":comparison,
              "clock_scores":clock_scores,"transport_scores":transport_scores,
              "previous_batch_qualification":frozen["batch"]["qualification"],
              "source_sha256":{str(p.relative_to(ROOT)):sha(p) for p in files},
              "output_sha256":{p.name:sha(p) for p in sorted(output.iterdir())}}
    write_json(output/"summary.json",report)
    return report


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,default=OUTPUT/datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ"))
    args = parser.parse_args()
    t = perf_counter()
    record = run(args.output.resolve())
    for case,scores in record["central_scores"].items():
        print(case, ", ".join(f"{feed}: {scores[feed]['rmse_percentage_points']:.2f} pp"
                             for feed in ("1to1","2to1","joint_40")), flush=True)
    print(f"Completed {args.output} in {perf_counter()-t:.1f} s",flush=True)
