"""Fit reversible composition and catalyst-age hypotheses on whole batch curves."""

from argparse import ArgumentParser
from datetime import datetime, timezone
from itertools import product
from math import pi, sqrt
from pathlib import Path
from time import perf_counter

import numpy as np

from examples.biocatalysis.pentyl_propanoate import run as base
from examples.biocatalysis.pentyl_propanoate.run import array, error_metrics, read_rows, write_rows
from examples.biocatalysis.pentyl_propanoate_loading import run as loading
from examples.biocatalysis.pentyl_propanoate_loading_evidence import run as evidence
from pore2process_demo.bridge import runtime
from pore2process_demo.io import read_json,sha,write_json
from pore2process_demo.paths import ROOT

HERE=Path(__file__).resolve().parent
PLAN=HERE/"protocol.json"
SOURCE=ROOT/"literature/sources"
BATCH_FILES=(base.BATCH,evidence.BATCH_1,SOURCE/"lorenzo2026/batch_1to1_temperature.csv")
EQUILIBRIUM_FILES=(base.EQUILIBRIUM,SOURCE/"lorenzo2026/apparent_equilibrium_1to1.csv")
ACID_SOURCE=SOURCE/"vanhecke2024/acid_inhibition.json"
OUTPUT=ROOT/"results/biocatalysis/pentyl_propanoate_composition"
LAWS=("ratio_reference","reversible_shared","acid_analogue","acid","composition",
      "acid_thermal","acid_decay","acid_activation","acid_both","protected_decay","protected_both",
      "residual_decay","protected_residual_decay")
STATIC=LAWS[:6]
CLOCK=LAWS[6:]
VARIANTS={
    "primary":("printed","nominal_ratio",True),
    "point_weights":("printed","nominal_ratio",False),
    "printed_initial_ratio":("printed","printed_ratio",True),
    "from_acid":("from_acid","nominal_ratio",True),
    "from_ester":("from_ester","nominal_ratio",True),
}
PARAMETERS=("k2","k2_slope","k1","k1_slope","acid_b","alcohol_b","product_slope",
    "acid_heat_slope","hazard_ref_min_inv","hazard_temperature_slope","initial_deficit",
    "recovery_min","acid_exponent","hazard_exponent","protection_exponent","residual_active_fraction")


def read_batch(plan):
    groups={}
    expected=plan["source"]["assays"]
    nominal={i:(2 if i<=4 or 9<=i<=12 else 1) for i in expected}
    source=[r for path in BATCH_FILES for r in read_rows(path)]
    for assay in expected:
        rows=sorted([r for r in source if int(r["assay_id"])==assay],key=lambda r:float(r["time_min"]))
        if len(rows)<5 or float(rows[0]["time_min"])!=0:
            raise ValueError(f"Missing time-zero record or campaign: {assay}")
        feed=nominal[assay]
        w=(2.5,5,10,20)[(assay-1)%4] if assay<=8 else 10
        T=40 if assay<=8 else (50+10*((assay-9)%4))
        A0=float(rows[0]["x_acid"])
        R0=float(rows[0]["x_alcohol"])/A0
        t=array(rows,"time_min")
        if (not np.all(np.diff(t)>0) or not np.all(array(rows,"temperature_C")==T)
                or not np.all(array(rows,"catalyst_mass_pct_of_initial_acid")==w)):
            raise ValueError(f"Assay {assay} does not match Table 1.")
        data=[{
            "assay":assay,"feed":feed,"temperature_C":T,"loading_pct":w,
            "nominal_ratio":float(feed),"printed_ratio":R0,"time_min":float(row["time_min"]),
            "printed":float(row["acid_conversion"]),"from_acid":1-float(row["x_acid"])/A0,
            "from_ester":float(row["x_ester"])/A0,
            "mole_fraction_sum":sum(float(row[k]) for k in ("x_acid","x_alcohol","x_ester","x_water")),
        } for row in rows]
        for response in ("printed","from_acid","from_ester"):
            x=array(data,response)
            if x[0]!=0 or not np.all(np.isfinite(x)&(x>=0)&(x<=1)):
                raise ValueError(f"Nonphysical conversion response in assay {assay}: {response}")
        groups[assay]=data
    if set(groups)!=set(expected) or sum(map(len,groups.values()))!=225 or len(source)!=225:
        raise ValueError("Expected all 225 Table S2 enzyme records and no resin data.")
    return groups


def equilibrium_rows():
    return [{**{k:float(v) for k,v in r.items()},"feed":feed}
            for feed,path in ((2,EQUILIBRIUM_FILES[0]),(1,EQUILIBRIUM_FILES[1]))
            for r in read_rows(path)]


def equilibrium_training(model,rows,training):
    result={}
    for feed in (1,2):
        selected=[r for r in rows if r["feed"]==feed and int(r["assay_id"]) in training]
        if len(selected)<2:
            raise ValueError(f"Feed {feed} lacks independent training equilibrium temperatures.")
        fit=model.fit_equilibrium(array(selected,"temperature_C"),array(selected,"apparent_equilibrium_constant_Kx"))
        result[str(feed)]={
            "assays":[int(r["assay_id"]) for r in selected],
            "log_K_reference":float(fit.log_K_reference),
            "inverse_temperature_slope":float(fit.inverse_temperature_slope),
            "log_K_rmse":float(fit.log_K_rmse),
            "calibration_temperatures_C":[float(r["temperature_C"]) for r in selected],
        }
    return result


def equilibrium_value(model,fit,feed,T):
    r=fit[str(feed)]
    return float(model.equilibrium_at(r["log_K_reference"],r["inverse_temperature_slope"],T))


def julia_curves(model,groups,ids,eq,variant):
    response,ratio,balanced=VARIANTS[variant]
    selected=[groups[i][1:] for i in ids]
    heads=[r[0] for r in selected]
    flat=[r for curve in selected for r in curve]
    K=[equilibrium_value(model,eq,r["feed"],r["temperature_C"]) for r in heads]
    curves=model.make_curves(np.array(ids,dtype=int),
        np.array([r["feed"] for r in heads],dtype=int),array(heads,"temperature_C"),
        array(heads,ratio),array(heads,"loading_pct"),np.array(K),
        np.array([len(r) for r in selected],dtype=int),
        array(flat,"time_min"),array(flat,response),equal_curve_weights=balanced)
    return curves,flat


def fit_record(model,result):
    return {
        "model":str(result.name),
        "parameters":dict(zip(PARAMETERS,map(float,model.parameter_vector(result.parameters)),strict=True)),
        "optimizer_coordinates":dict(zip(map(str,result.names),map(float,result.coordinates),strict=True)),
        "bounds":dict(zip(map(str,result.names),map(str,result.active_bounds),strict=True)),
        "weighted_sse":float(result.weighted_sse),
        "converged":bool(result.converged),
        "projected_gradient_inf":float(result.projected_gradient_inf),
        "iterations":int(result.iterations),
        "numerical_rank":int(result.numerical_rank),
        "singular_values":list(map(float,result.singular_values)),
        "start_sse":list(map(float,result.start_sse)),
        "max_accepted_steps":int(result.max_accepted_steps),
        "max_rejected_steps":int(result.max_rejected_steps),
    }


def parameters_from_record(model,record):
    return model.Parameters(*[record["parameters"][p] for p in PARAMETERS])


def scores(rows):
    return error_metrics(array(rows,"predicted_conversion")-array(rows,"observed_conversion"))


def equal_curve_rmse(by_assay,ids):
    chosen=[by_assay[str(a)]["rmse_percentage_points"] for a in ids]
    return float(np.sqrt(np.mean(np.square(chosen))))


def fit_case(model,groups,eq_rows,name,variant,train,test,*,laws=LAWS,
             acid_exponent=3.5,hazard_exponent=3.5,analogue_KI=1.70,
             fixed_protection=None,refine=False):
    print(f"{variant}: {name} ({len(train)} training curves, {len(test)} validation curves, {len(laws)} models)",flush=True)
    eq=equilibrium_training(model,eq_rows,train)
    train_curves,_=julia_curves(model,groups,train,eq,variant)
    validation,rows=julia_curves(model,groups,test,eq,variant)
    response,ratio,_=VARIANTS[variant]
    fits=model.fit_models(train_curves,names=list(laws),acid_exponent=acid_exponent,
        hazard_exponent=hazard_exponent,analogue_KI=analogue_KI,
        fixed_protection=float("nan") if fixed_protection is None else fixed_protection)
    result,predictions={},[]
    for fit in fits:
        record=fit_record(model,fit)
        law=record["model"]
        y=model.predict_curves(fit.parameters,validation)
        refined_change=None
        if refine:
            tight=model.predict_curves(fit.parameters,validation,rtol=2e-9,atol=2e-11,max_step_min=10.)
            refined_change=float(np.max(abs(np.asarray(y.conversion)-np.asarray(tight.conversion))))
            if refined_change>1e-5:
                raise RuntimeError(f"Insufficient batch integration accuracy for {name}, {law}: {refined_change}")
        data=[{
            "case":name,"variant":variant,"law":law,"assay":r["assay"],"feed_ratio":r["feed"],
            "temperature_C":r["temperature_C"],"loading_pct":r["loading_pct"],"time_min":r["time_min"],
            "observed_conversion":r[response],"predicted_conversion":float(x),
            "relative_catalyst_activity":float(a),"integrated_damage_hazard":float(h),
        } for r,x,a,h in zip(rows,y.conversion,y.activity,y.hazards,strict=True)]
        if not all(0<=r["predicted_conversion"]<=1 and 0<=r["relative_catalyst_activity"]<=1 for r in data):
            raise RuntimeError("Nonphysical batch solution.")
        by_assay={str(i):scores([r for r in data if r["assay"]==i]) for i in test}
        result[law]=record | {
            "scores_by_assay":by_assay,"equal_assay_rmse_pp":equal_curve_rmse(by_assay,test),
            "point_weighted_scores":scores(data),"maximum_refinement_absolute_conversion":refined_change,
            "final_activity_by_assay":{str(i):next(r["relative_catalyst_activity"] for r in reversed(data) if r["assay"]==i)
                                       for i in test},
            "validation_max_accepted_steps":int(y.max_accepted_steps),
            "validation_max_rejected_steps":int(y.max_rejected_steps),
        }
        predictions.extend(data)
    return {"training_assays":list(train),"validation_assays":list(test),
        "equilibrium_training":eq,"fit":result},predictions


def fold_scores(folds,laws):
    merged={}
    for law in laws:
        by_assay={a:s for fold in folds.values() for a,s in fold["fit"][law]["scores_by_assay"].items()}
        ids=list(map(int,by_assay))
        target=[a for a in ids if a<=4 or 9<=a<=12]
        hot=[a for a in ids if a in (15,16)]
        active=[a for a in ids if a<15]
        merged[law]={
            "scores_by_assay":by_assay,"overall_equal_assay_rmse_pp":equal_curve_rmse(by_assay,ids),
            "active_equal_assay_rmse_pp":equal_curve_rmse(by_assay,active),
            "two_to_one_equal_assay_rmse_pp":equal_curve_rmse(by_assay,target),
            "hot_equimolar_equal_assay_rmse_pp":equal_curve_rmse(by_assay,hot) if hot else None,
            "all_fits_converged":all(fold["fit"][law]["converged"] for fold in folds.values()),
            "all_parameters_identified":all(
                fold["fit"][law]["numerical_rank"]==len(fold["fit"][law]["optimizer_coordinates"])
                and "upper" not in fold["fit"][law]["bounds"].values() for fold in folds.values()),
        }
    return merged


def screening(merged,sensitivities):
    null=merged["ratio_reference"]
    targets=[1,2,3,4,9,10,11,12]
    results={}
    for law,record in merged.items():
        if law=="ratio_reference":
            continue
        ratio=1-record["overall_equal_assay_rmse_pp"]/null["overall_equal_assay_rmse_pp"]
        target=1-record["two_to_one_equal_assay_rmse_pp"]/null["two_to_one_equal_assay_rmse_pp"]
        def error(a,x):
            return x["scores_by_assay"][str(a)]["rmse_percentage_points"]
        worse=max(error(a,record)-error(a,null) for a in targets)
        high=error(8,record)<error(8,null)
        hot=all(error(a,record)<error(a,null) for a in (15,16))
        robust=True
        for variant,cases in sensitivities.items():
            for case,fold in cases.items():
                if case=="all_data":
                    continue
                if law not in fold["fit"]:
                    continue
                a,b=fold["fit"][law],fold["fit"]["ratio_reference"]
                # A wholesale reversal beyond the reading-level tolerance
                # on either target critical condition is a failed robustness screen.
                if any(error(i,a)-error(i,b)>.2 for i in fold["validation_assays"]):
                    robust=False
        results[law]={"overall_relative_improvement":ratio,"two_to_one_relative_improvement":target,
            "worst_two_to_one_assay_worsening_pp":worse,"improved_highest_equimolar_loading":high,
            "improved_both_hot_equimolar_assays":hot,"critical_sensitivities_not_reversed":robust,
            "all_fits_converged":record["all_fits_converged"],
            "all_parameters_identified":record["all_parameters_identified"],
            "predictive_screen_passed":bool(ratio>=.1 and target>=.1 and worse<=.2 and high and hot
                                            and robust and record["all_fits_converged"]
                                            and record["all_parameters_identified"])}
    return results


def fit_campaign(model,plan,groups,eq_rows,source,*,quick=False):
    ids=sorted(groups)
    active=[i for i in ids if i<15]
    options={"acid_exponent":source["inhibition_exponent"],"hazard_exponent":3.5,
             "analogue_KI":source["recalculated_inhibition_concentration_mol_per_kg_liquid_substrates"]}
    predictions=[]
    def run_case(name,variant,train,test,**kwargs):
        settings=options|kwargs
        fold,rows=fit_case(model,groups,eq_rows,name,variant,train,test,**settings)
        predictions.extend(rows)
        return fold
    full=run_case("all_data","primary",ids,ids,refine=True)
    if quick:
        return {"all_data":full},predictions
    folds={}
    for name,test in plan["holdouts"]["paired_folds"].items():
        folds[name]=run_case(name,"primary",[a for a in ids if a not in test],test,
                            refine=name in ("loading_20_40C","temperature_80C"))
    challenge=run_case("hot_acid_challenge","primary",active,[15,16],refine=True)
    active_folds={}
    for name,held in plan["holdouts"]["paired_folds"].items():
        test=[i for i in held if i in active]
        active_folds[name]=run_case("active_"+name,"primary",[a for a in active if a not in test],test,laws=STATIC)
    sensitivities={}
    for variant in list(VARIANTS)[1:]:
        sensitivities[variant]={}
        for name in plan["holdouts"]["critical_sensitivities"]:
            test=ids if name=="all_data" else plan["holdouts"]["paired_folds"][name]
            train=ids if name=="all_data" else [i for i in ids if i not in test]
            sensitivities[variant][name]=run_case(name,variant,train,test)
    # Structural sensitivity for clock/inhibition shape exponents. The
    # other models are not changed. Use the same training boundaries.
    exponent_sensitivities={}
    for ex in (2.,5.):
        key=f"hazard_exponent_{ex:g}"
        exponent_sensitivities[key]={}
        for name in ("all_data","temperature_80C"):
            test=ids if name=="all_data" else plan["holdouts"]["paired_folds"][name]
            train=ids if name=="all_data" else [i for i in ids if i not in test]
            exponent_sensitivities[key][name]=run_case(f"{key}_{name}","primary",train,test,
                laws=("acid","acid_decay","acid_activation","acid_both"),hazard_exponent=ex)
    for ex in (2.,5.):
        key=f"acid_exponent_{ex:g}"
        exponent_sensitivities[key]={}
        for name in ("all_data","loading_20_40C"):
            test=ids if name=="all_data" else plan["holdouts"]["paired_folds"][name]
            train=ids if name=="all_data" else [i for i in ids if i not in test]
            exponent_sensitivities[key][name]=run_case(f"{key}_{name}","primary",train,test,
                laws=("acid","acid_decay","acid_both"),acid_exponent=ex)
    for ex in (2.,4.,8.):
        key=f"fixed_alcohol_protection_{ex:g}"
        exponent_sensitivities[key]={}
        for name in ("all_data","temperature_80C"):
            test=ids if name=="all_data" else plan["holdouts"]["paired_folds"][name]
            train=ids if name=="all_data" else [i for i in ids if i not in test]
            exponent_sensitivities[key][name]=run_case(f"{key}_{name}","primary",train,test,
                laws=("acid","acid_decay","protected_decay","protected_both","residual_decay","protected_residual_decay"),
                fixed_protection=ex)
    fixed_protection_folds={}
    for ex in (4.,8.):
        label=f"fixed_alcohol_protection_{ex:g}"
        cases={}
        for name,test in plan["holdouts"]["paired_folds"].items():
            if name=="temperature_80C":
                cases[name]=exponent_sensitivities[label][name]
            else:
                cases[name]=run_case(f"{label}_{name}","primary",[i for i in ids if i not in test],test,
                    laws=("acid","acid_decay","protected_decay","residual_decay","protected_residual_decay"),
                    fixed_protection=ex)
        fixed_protection_folds[label]={
            "folds":cases,"scores":fold_scores(cases,("protected_decay","protected_residual_decay"))}
    merged=fold_scores(folds,LAWS)
    active_merged=fold_scores(active_folds,STATIC)
    return {"all_data":full,"paired_holdouts":folds,"paired_scores":merged,
        "hot_acid_challenge":challenge,"active_domain_holdouts":active_folds,
        "active_domain_scores":active_merged,"sensitivities":sensitivities,
        "exponent_sensitivities":exponent_sensitivities,
        "fixed_protection_paired_sensitivity":fixed_protection_folds,
        "qualification":screening(merged,sensitivities)},predictions


def load_runtime():
    jl=runtime()
    for name,path in (("PentylPropanoate",base.HERE/"model.jl"),
            ("CatalystLoading",loading.HERE/"model.jl"),
            ("CompositionKinetics",HERE/"model.jl"),
            ("LorenzoPackedBed",loading.published.MODEL),
            ("BeadFilm",loading.film.HERE/"model.jl"),
            ("CompositionBedTransfer",HERE/"bed_transfer.jl")):
        if not jl.seval(f"isdefined(@__MODULE__, :{name})"):
            jl.include(str(path))
    return jl


def batch_summary():
    plan=read_json(PLAN)
    group=read_batch(plan)
    jl=load_runtime()
    results,pred=fit_campaign(jl.CompositionKinetics,plan,group,equilibrium_rows(),read_json(ACID_SOURCE),quick=True)
    for name,f in results["all_data"]["fit"].items():
        print(name,round(f["equal_assay_rmse_pp"],3),f["parameters"],f["converged"],f["projected_gradient_inf"],flush=True)
    return results,pred


def bed_conditions(jl,source,liquid):
    bed=source["bed"]
    if bed["catalyst_mass_kg"]!=.006:
        raise ValueError("Do not change the measured six gram catalyst inventory.")
    area=pi*bed["diameter_m"]**2/4
    void,diameter=bed["bed_porosity"],bed["bead_mean_diameter_m"]
    volume_L=area*bed["length_m"]*void*1000
    catalyst_g_L=1000*bed["catalyst_mass_kg"]/volume_L
    a_l=6*(1-void)/(diameter*void)
    pe=bed["length_m"]*.20/(void**2*diameter)
    acid,alcohol=bed["feed_acid_mol_per_L"],bed["feed_alcohol_mol_per_L"]
    by_T={round(r["temperature_K"]-273.15):r for r in liquid["conditions"]}
    cases=[]
    for T,Q in product(bed["reaction_temperature_C"],bed["feed_flow_ml_per_min"]):
        prop=by_T[T]
        coef=jl.BeadFilm.film_coefficients(Q/60e6/area,diameter,prop["density_kg_m3"],
            prop["viscosity_Pa_s"],prop["diffusivity_m2_s"]["central"],void,"wilson_geankoplis")
        eta=source["effectiveness_by_temperature"][str(T)]["eta_internal"]
        cases.append({"temperature_C":T,"flow_ml_min":Q,"catalyst_mass_kg":bed["catalyst_mass_kg"],
            "catalyst_g_L":catalyst_g_L,"acid0":acid,"alcohol0":alcohol,
            "residence_min":volume_L*1000/Q,"capacity_mol_L_min":60*float(coef.coefficient_m_s)*a_l*acid,
            "eta_internal":eta,"peclet":pe})
    return cases


def frozen_solve(jl,p,K,activity,condition,plan):
    model=jl.CompositionBedTransfer
    limits=plan["reactor"]
    c=condition
    args=(p,c["temperature_C"],K,activity,c["residence_min"],c["catalyst_g_L"],
          c["capacity_mol_L_min"],c["eta_internal"],c["peclet"])
    kwargs={"acid0":c["acid0"],"alcohol0":c["alcohol0"]}
    coarse,fine=[model.frozen_prediction(*args,**kwargs,cells=n) for n in limits["cells"]]
    grid=float(abs(coarse.conversion-fine.conversion))
    balance=max(float(coarse.balance_error),float(fine.balance_error))
    film=max(float(coarse.maximum_sampled_film_residual_mol_L_min),float(fine.maximum_sampled_film_residual_mol_L_min))
    slope=max(float(coarse.maximum_positive_rate_slope_over_film_capacity),
              float(fine.maximum_positive_rate_slope_over_film_capacity))
    if grid>limits["maximum_grid_conversion_difference"] or balance>limits["maximum_normalized_reaction_flux_balance_error"] or film>limits["maximum_film_balance_error_mol_L_min"]:
        raise RuntimeError(f"Failed frozen reactor refinement/balance: {c}, grid={grid}, balance={balance}, film={film}")
    if slope>=.99:
        raise RuntimeError(f"Film root monotonicity needs further investigation: slope/capacity={slope}")
    return {"predicted_conversion":float(fine.conversion),"grid_change_absolute_conversion":grid,
        "normalized_balance_error":balance,"film_balance_error_mol_L_min":film,
        "max_positive_rate_slope_over_film_capacity":slope,
        "film_slope_sampling_refinement_change":float(fine.film_slope_refinement_ratio_change),
        "surface_progress_at_inlet":float(fine.surface_progress_at_inlet)}


def replay_reference(jl,conditions,source,plan):
    reactor,film=jl.LorenzoPackedBed,jl.BeadFilm
    kinetic=reactor.PingPong(np.array([source["rate_parameters"][k] for k in loading.published.PARAMETERS]))
    rows=[]
    for c in conditions:
        f=film.FilmKinetics(kinetic,c["capacity_mol_L_min"],c["eta_internal"])
        a,b=[reactor.packed_bed(f,c["temperature_C"],c["residence_min"],c["catalyst_g_L"],
            c["eta_internal"],c["peclet"],acid0=c["acid0"],alcohol0=c["alcohol0"],cells=n)
            for n in plan["reactor"]["cells"]]
        local=film.inlet_balance(f,c["temperature_C"],c["catalyst_g_L"],acid0=c["acid0"],alcohol0=c["alcohol0"])
        rows.append({"temperature_C":c["temperature_C"],"flow_ml_min":c["flow_ml_min"],
            "case":"prior_central_film_reference","actual_mass_kg":c["catalyst_mass_kg"],"activity":1.,
            "apparent_Kx":None,"predicted_conversion":float(b.conversion),
            "grid_change_absolute_conversion":float(abs(a.conversion-b.conversion)),
            "normalized_balance_error":max(float(a.balance_error),float(b.balance_error)),
            "film_balance_error_mol_L_min":float(local.film_balance_error_mol_L_min),
            "max_positive_rate_slope_over_film_capacity":None,
            "film_slope_sampling_refinement_change":None,
            "surface_progress_at_inlet":float(1-local.surface_acid_relative_to_bulk)})
    # Verify previous predictions only; no reactor measurements inform a fit.
    previous={loading.published.condition(r):float(r["wilson_geankoplis_central"])
              for r in read_rows(loading.FILM_REFERENCE)}
    replay=max(abs(r["predicted_conversion"]-previous[loading.published.condition(r)]) for r in rows)
    if replay>1e-10:
        raise RuntimeError(f"Central-film replay changed {replay}")
    return rows,replay


def bed_campaign(jl,campaign,plan):
    print("Replaying and predicting fixed-inventory reactor scenarios",flush=True)
    source=read_json(loading.published.MODEL_INPUTS)
    conditions=bed_conditions(jl,source,read_json(loading.film.PROPERTIES))
    results,replay=replay_reference(jl,conditions,source,plan)
    sampled=[]
    for training,fit in (("all_batch",campaign["all_data"]),("active_batch",campaign["hot_acid_challenge"])):
        laws=LAWS if training=="all_batch" else STATIC
        for law in laws:
            case=f"{training}_{law}"+("_fully_active" if law in CLOCK else "")
            print(f"  {case}",flush=True)
            p=parameters_from_record(jl.CompositionKinetics,fit["fit"][law])
            for c in conditions:
                K=equilibrium_value(jl.CompositionKinetics,fit["equilibrium_training"],2,c["temperature_C"])
                solution=frozen_solve(jl,p,K,1.,c,plan)
                results.append({"temperature_C":c["temperature_C"],"flow_ml_min":c["flow_ml_min"],
                    "case":case,"actual_mass_kg":c["catalyst_mass_kg"],"activity":1.,"apparent_Kx":K,**solution})
    # Long-time nonzero activity is a model hypothesis, not a measured
    # fully wetted state at three hours. Include its limit explicitly.
    for law in ("residual_decay","protected_residual_decay"):
        fit=campaign["all_data"]
        record=fit["fit"][law]
        p=parameters_from_record(jl.CompositionKinetics,record)
        fraction=record["parameters"]["residual_active_fraction"]
        for c in conditions:
            K=equilibrium_value(jl.CompositionKinetics,fit["equilibrium_training"],2,c["temperature_C"])
            solution=frozen_solve(jl,p,K,fraction,c,plan)
            results.append({"temperature_C":c["temperature_C"],"flow_ml_min":c["flow_ml_min"],
                "case":f"all_batch_{law}_long_time_floor","actual_mass_kg":c["catalyst_mass_kg"],
                "activity":fraction,"apparent_Kx":K,**solution})
    # Stationary-clock envelopes deliberately do not use liquid residence
    # as enzyme age. No calibration of startup age to measured conversion.
    full=campaign["all_data"]
    for law in CLOCK:
        print(f"  Stationary catalyst age screen: {law}",flush=True)
        p=parameters_from_record(jl.CompositionKinetics,full["fit"][law])
        for c in conditions:
            T,Q=c["temperature_C"],c["flow_ml_min"]
            K=equilibrium_value(jl.CompositionKinetics,full["equilibrium_training"],2,T)
            s=jl.CompositionKinetics.static_terms(p,T,2.,K,feed=2)
            for age in (30.,180.):
                envelope=jl.CompositionKinetics.activity_bounds(s,age)
                lo,hi=float(envelope.lower),float(envelope.upper)
                lower=frozen_solve(jl,p,K,lo,c,plan)
                upper=lower if lo==hi else frozen_solve(jl,p,K,hi,c,plan)
                if lower["predicted_conversion"]>upper["predicted_conversion"]+1e-8:
                    raise RuntimeError("Frozen uniform activity did not have monotone reactor response.")
                sampled.append({"law":law,"temperature_C":T,"flow_ml_min":Q,"catalyst_age_min":age,
                    "inlet_composition_activity":lo,"equilibrium_composition_activity":hi,
                    "inlet_hazard_min_inv":float(envelope.inlet_hazard_min_inv),
                    "equilibrium_hazard_min_inv":float(envelope.equilibrium_hazard_min_inv),
                    "recovery_fraction":float(envelope.recovery_fraction),
                    "conversion_uniform_inlet_activity":lower["predicted_conversion"],
                    "conversion_uniform_equilibrium_activity":upper["predicted_conversion"],
                    "maximum_grid_change":max(lower["grid_change_absolute_conversion"],upper["grid_change_absolute_conversion"]),
                    "maximum_balance_error":max(lower["normalized_balance_error"],upper["normalized_balance_error"]),
                    "maximum_film_balance":max(lower["film_balance_error_mol_L_min"],upper["film_balance_error_mol_L_min"])})
    bounds={}
    for law in CLOCK:
        records={}
        for c in conditions:
            T,Q=c["temperature_C"],c["flow_ml_min"]
            a,b=[next(r for r in sampled if r["law"]==law and r["temperature_C"]==T
                       and r["flow_ml_min"]==Q and r["catalyst_age_min"]==age) for age in (30.,180.)]
            records[f"{T}C_{Q:g}mL_min"]={
                "inlet_exposure_activity_retention_30_to_180":b["inlet_composition_activity"]/a["inlet_composition_activity"]
                      if a["inlet_composition_activity"]>0 else None,
                "uniform_inlet_exposure_conversion_drift_pp":100*(b["conversion_uniform_inlet_activity"]-a["conversion_uniform_inlet_activity"]),
                "uniform_equilibrium_exposure_conversion_drift_pp":100*(b["conversion_uniform_equilibrium_activity"]-a["conversion_uniform_equilibrium_activity"]),
            }
        bounds[law]=records
    diagnostics={"previous_central_film_replay_difference":replay,
        "max_bed_grid_conversion_change":max([r["grid_change_absolute_conversion"] for r in results]+[r["maximum_grid_change"] for r in sampled]),
        "max_bed_normalized_balance_error":max([r["normalized_balance_error"] for r in results]+[r["maximum_balance_error"] for r in sampled]),
        "max_film_balance_error_mol_L_min":max([r["film_balance_error_mol_L_min"] for r in results]+[r["maximum_film_balance"] for r in sampled]),
        "max_sampled_positive_rate_slope_over_film_capacity":max(r["max_positive_rate_slope_over_film_capacity"]
            for r in results if r["max_positive_rate_slope_over_film_capacity"] is not None),
        "clock_age_scope":"Uniform activity endpoints are steady, fully wetted exposure scenarios. Actual bead history and transient concentration field are not predicted; do not select an age by fitting the bed.",
        "stationary_clock_drift_scenarios":bounds,
        "actual_mass_kg":conditions[0]["catalyst_mass_kg"],
        "catalyst_per_liquid_volume_g_L":conditions[0]["catalyst_g_L"],
        "reactor_conditions":conditions}
    return results,sampled,diagnostics


def compare_reactor(rows):
    source=read_rows(loading.published.HOLDOUT)
    observed={loading.published.condition(r):float(r["acid_conversion"]) for r in source}
    if len(observed)!=20 or len(source)!=20:
        raise ValueError("Every one of the 20 reactor measurements remains in the comparison.")
    results=[r | {"observed_conversion":observed[loading.published.condition(r)]} for r in rows]
    summary={}
    for case in sorted({r["case"] for r in results}):
        data=[r for r in results if r["case"]==case]
        if set(map(loading.published.condition,data))!=set(observed) or len(data)!=20:
            raise ValueError(f"Incomplete reactor prediction: {case}")
        summary[case]=scores(data)
    return results,summary


def progress_audit(groups,predictions):
    pairs=[]
    for tlo,thi in ((19.9,10.),(60.,30.)):
        low=next(r for r in groups[7] if r["time_min"]==tlo)
        high=next(r for r in groups[8] if r["time_min"]==thi)
        pairs.append({"nominal_1to1_40C_low_loading_time_min":tlo,"high_loading_time_min":thi,
            "low_loading_pct_min":10*tlo,"high_loading_pct_min":20*thi,
            **{f"{key}_conversion_high_minus_low_pp":100*(high[key]-low[key])
               for key in ("printed","from_acid","from_ester")}})
    tail=[]
    for a in (15,16):
        start=180.1
        for law in CLOCK:
            values=[r for r in predictions if r["case"]=="all_data" and r["variant"]=="primary"
                    and r["law"]==law and r["assay"]==a and r["time_min"]>=start]
            first,last=values[0],values[-1]
            tail.append({"assay":a,"law":law,"start_min":first["time_min"],"end_min":last["time_min"],
                "observed_conversion_increase_pp":100*(last["observed_conversion"]-first["observed_conversion"]),
                "predicted_conversion_increase_pp":100*(last["predicted_conversion"]-first["predicted_conversion"]),
                "model_final_activity":last["relative_catalyst_activity"]})
    return {"near_equal_mass_time_progress_pairs":pairs,"hot_equimolar_late_tail":tail}


def run(output):
    plan=read_json(PLAN)
    groups=read_batch(plan)
    jl=load_runtime()
    campaign,predictions=fit_campaign(jl.CompositionKinetics,plan,groups,equilibrium_rows(),read_json(ACID_SOURCE))
    # Save batch outputs before bed solves as a reproducible diagnostic
    # record even if a later nonlinear reactor branch needs investigation.
    output.mkdir(parents=True,exist_ok=False)
    write_rows(output/"batch_predictions.csv",predictions)
    write_json(output/"batch_fit.json",campaign)
    bed,clock,numerics=bed_campaign(jl,campaign,plan)
    write_rows(output/"predictions_before_comparison.csv",bed)
    write_rows(output/"stationary_clock_screen.csv",clock)
    summary={"batch":campaign,"progress_audit":progress_audit(groups,predictions),"printing_consistency":{str(a):{
        "records_including_start":len(r),"max_abs_mole_fraction_sum_error":max(abs(x["mole_fraction_sum"]-1) for x in r),
        "max_printed_vs_acid":max(abs(x["printed"]-x["from_acid"]) for x in r),
        "max_printed_vs_ester":max(abs(x["printed"]-x["from_ester"]) for x in r)} for a,r in groups.items()},
        "bed_numerics":numerics,"julia":str(jl.seval("VERSION"))}
    write_json(output/"frozen_calibration.json",summary)
    # The measured reactor file is consumed only here for scoring.
    comparison,metrics=compare_reactor(bed)
    write_rows(output/"bed_comparison.csv",comparison)
    summary["bed_comparison"]=metrics
    summary["sixty_C_three_ml_min"]=[r for r in comparison if loading.published.condition(r)==(60,3)]
    trace=ROOT/"examples/biocatalysis/pentyl_propanoate_conditioning/reference/summary.json"
    summary["independent_bed_time_trace"]=read_json(trace)["bed_time_trace"]
    files=[PLAN,HERE/"model.jl",HERE/"bed_transfer.jl",HERE/"run.py",*BATCH_FILES,*EQUILIBRIUM_FILES,
        ACID_SOURCE,SOURCE/"lorenzo2026/temperature_series.md",SOURCE/"vanhecke2024/README.md",
        SOURCE/"ascorbyl2019/README.md",
        base.HERE/"model.jl",base.HERE/"run.py",loading.HERE/"model.jl",loading.published.MODEL,
        loading.published.MODEL_INPUTS,loading.published.HOLDOUT,loading.film.PROPERTIES,
        loading.film.HERE/"model.jl",loading.FILM_REFERENCE,trace]
    summary["source_sha256"]={str(path.relative_to(ROOT)):sha(path) for path in files}
    summary["output_sha256"]={p.name:sha(p) for p in sorted(output.iterdir())}
    write_json(output/"summary.json",summary)
    return summary


if __name__=="__main__":
    parser=ArgumentParser(description=__doc__)
    stamp=datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    parser.add_argument("--output",type=Path,default=OUTPUT/stamp)
    args=parser.parse_args()
    started=perf_counter()
    report=run(args.output.resolve())
    print("Output:",args.output,flush=True)
    for law,fit in report["batch"]["paired_scores"].items():
        print(law,"paired held-out RMSE",round(fit["overall_equal_assay_rmse_pp"],3),
              "2:1",round(fit["two_to_one_equal_assay_rmse_pp"],3),"pp",flush=True)
    for case,metrics in report["bed_comparison"].items():
        print(case,"reactor RMSE",round(metrics["rmse_percentage_points"],3),"pp",flush=True)
    print("Elapsed seconds:",round(perf_counter()-started,1),flush=True)
