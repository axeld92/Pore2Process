"""Check immutable prediction partitions, then score both withheld reactor feeds."""

from argparse import ArgumentParser
from collections import Counter
from itertools import product
from pathlib import Path

import numpy as np

from examples.biocatalysis.pentyl_propanoate.run import read_rows, write_rows
from examples.biocatalysis.pentyl_propanoate_history import run
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT

STRINGS = {"task_id","scenario","law"}
INTS = {"task_index","feed_ratio","temperature_C","fine_cells","coarse_cells","numerical_solve_count",
        "accepted_steps","rejected_steps","max_newton_iterations"}
FLOWS = (.7,1.,2.,3.)
TEMPS = (40,50,60,70,80)


def typed(row):
    return {key:value if key in STRINGS else int(value) if key in INTS else float(value)
            for key,value in row.items()}


def identity(row):
    return row["scenario"],row["law"],*run.eq.identity(row)


def predicted_ledgers(part_paths,plan):
    rows=[]; nums=[]; profiles={}; manifests=[]
    wanted=run.sources()
    for path in part_paths:
        m=read_json(path/"manifest.json")
        if m["sources"]!=wanted or m["failed_cases"] or m["completed_count"]!=m["expected_count"]:
            raise RuntimeError(f"Incomplete or inconsistent prediction partition: {path}")
        for name,digest in m["ledger_sha256"].items():
            if sha(path/name)!=digest:
                raise RuntimeError(f"Prediction ledger changed: {path/name}")
        for name,digest in m["profile_sha256"].items():
            target=path/name
            if target.stem in profiles or sha(target)!=digest:
                raise RuntimeError(f"Profile is duplicated or changed: {target}")
            profiles[target.stem]=target
        new=[typed(r) for r in read_rows(path/"numerics.csv")]
        new_rows=[typed(r) for r in read_rows(path/"predictions.csv")]
        if not all(r["task_index"]%m["parts"]==m["part"] for r in new):
            raise RuntimeError("Wrong deterministic partition index.")
        nums.extend(new); rows.extend(new_rows); manifests.append(m)
    expected={(scenario,law,feed,T,Q) for scenario,laws in plan["scenarios"].items()
              for law,feed,T,Q in product(laws,(1,2),TEMPS,FLOWS)}
    if (len(manifests)!=manifests[0]["parts"] or {m["parts"] for m in manifests}!={len(manifests)} or
            {m["part"] for m in manifests}!=set(range(len(manifests))) or
            {identity(r) for r in nums}!=expected or len(nums)!=560 or len(expected)!=560 or
            len(profiles)!=560 or set(profiles)!={r["task_id"] for r in nums}):
        raise RuntimeError("The complete partition/condition matrix was not recovered.")
    rows.sort(key=lambda r:(r["task_index"],r["time_min"])); nums.sort(key=lambda r:r["task_index"])
    groups={r["task_id"]:[] for r in nums}
    for r in rows:
        groups[r["task_id"]].append(r)
    if len(rows)!=560*len(plan["times_min"]) or any(
            [r["time_min"] for r in data]!=plan["times_min"] for data in groups.values()):
        raise RuntimeError("Missing time-on-stream snapshots.")
    if set(r["task_index"] for r in nums)!=set(range(560)) or any(r["actual_mass_kg"]!=.006 for r in rows+nums) or \
            len({(m["python"],m["julia"],m["platform"]) for m in manifests})!=1:
        raise RuntimeError("Invalid task ordering or catalyst mass.")
    return groups,rows,nums,profiles,manifests


def trace_reference():
    spec=read_json(run.CONDITIONING)["bed_time_trace"]
    dig=spec["digitization"]; y0,y100=dig["y_ticks_0_100_pct_px"]
    times=60*np.asarray(dig["times_h"]); observed=(y0-np.asarray(dig["marker_y_px"]))/(y0-y100)
    use=times>=30
    span=100*float(np.ptp(observed[use])); half=100*dig["reading_half_range_px"]/(y0-y100)
    return {"source_doi":spec["doi"],"source_location":spec["location"],"temperature_C":spec["temperature_C"],
            "flow_ml_min":spec["flow_ml_min"],"times_min":times.tolist(),"conversion":observed.tolist(),
            "span_30_to_180_pp":span,"change_30_to_180_pp":100*float(observed[-1]-observed[use][0]),
            "reading_half_range_pp":half,"reading_only_expanded_span_pp":span+2*half}


def score(groups,nums,plan):
    observed=run.eq.read_observations()  # first use of held-out values in this campaign
    prior={ (r["case"],*run.eq.identity(r)):float(r["predicted_conversion"])
            for r in read_rows(run.PREVIOUS/"bed_comparison.csv") }
    clocks={ (r["law"],*run.eq.identity(r),float(r["catalyst_age_min"])):r
             for r in read_rows(run.PREVIOUS/"clock_comparison.csv") }
    trace=trace_reference(); comp=[]; envelope=[]
    for d in nums:
        data=groups[d["task_id"]]
        feed,T,Q=run.eq.identity(d)
        y=observed[feed,T,Q]
        windows={r["time_min"]:r for r in data}
        use=[windows[t] for t in plan["reported_steady_sample_times_min"]]
        pred=np.array([r["predicted_conversion"] for r in use])
        a,b=windows[30],windows[180]
        name=f"all_batch_{d['law']}_fully_active"
        steady=prior[name,feed,T,Q]
        if d["scenario"]=="fresh_bulk":
            for time in (30,180):
                old=clocks[d["law"],feed,T,Q,time]; current=windows[time]
                lo,hi=float(old["inlet_activity"]),float(old["equilibrium_activity"])
                violation=max(0.,lo-current["minimum_activity"],current["maximum_activity"]-hi)
                if violation>1e-10:
                    raise RuntimeError("Local activity escaped the earlier uniform activity limits.")
                xlo,xhi=float(old["conversion_uniform_inlet_activity"]),float(old["conversion_uniform_equilibrium_activity"])
                envelope.append({"task_id":d["task_id"],"law":d["law"],"feed_ratio":feed,"temperature_C":T,
                    "flow_ml_min":Q,"age_min":time,"activity_bound_violation":violation,
                    "uniform_steady_low_conversion":xlo,"uniform_steady_high_conversion":xhi,
                    "transient_conversion":current["predicted_conversion"],
                    "outside_uniform_steady_interval_pp":100*max(0.,xlo-current["predicted_conversion"],current["predicted_conversion"]-xhi)})
        comp.append({"task_index":d["task_index"],"task_id":d["task_id"],"scenario":d["scenario"],"law":d["law"],
            "feed_ratio":feed,"temperature_C":T,"flow_ml_min":Q,"source_run_id":y["source_run_id"],
            "actual_mass_kg":d["actual_mass_kg"],"acid0_mol_L":d["acid0_mol_L"],"residence_min":d["residence_min"],
            "observed_conversion":y["observed_conversion"],"previous_fully_active_steady_conversion":steady,
            "predicted_10min_conversion":windows[10]["predicted_conversion"],"predicted_30min_conversion":float(pred[0]),
            "predicted_180min_conversion":float(pred[-1]),"predicted_30to180_sample_mean":float(np.mean(pred)),
            "minimum_30to180_conversion":float(np.min(pred)),"maximum_30to180_conversion":float(np.max(pred)),
            "span_30to180_pp":float(100*np.ptp(pred)),"change_30to180_pp":float(100*(pred[-1]-pred[0])),
            "observed_in_sample_window":bool(np.min(pred)<=y["observed_conversion"]<=np.max(pred)),
            "observed_in_window_with_half_point_reading_allowance":bool(np.min(pred)-.005<=y["observed_conversion"]<=np.max(pred)+.005),
            "mean_activity_30min":a["mean_activity"],"mean_activity_180min":b["mean_activity"],
            "first_quarter_activity_30min":a["first_quarter_activity"],"last_quarter_activity_30min":a["last_quarter_activity"],
            "first_quarter_activity_180min":b["first_quarter_activity"],"last_quarter_activity_180min":b["last_quarter_activity"],
            "inlet_cell_hazard_180min":b["inlet_cell_hazard"],"outlet_cell_hazard_180min":b["outlet_cell_hazard"],
            "nominal_ester_g_h_sample_mean":eq_ester(pred.mean(),d),"observed_nominal_ester_g_h":eq_ester(y["observed_conversion"],d),
            "generated_180min_mol":b["generated_mol"],"effluent_180min_mol":b["effluent_product_mol"],
            "backdiffusion_180min_mol":b["backdiffused_product_mol"],"liquid_product_180min_mol":b["liquid_product_mol"],
            "initial_liquid_product_mol":d["initial_liquid_product_mol"]})
    info={}
    for scenario,laws in plan["scenarios"].items():
        for law in laws:
            rows=[r for r in comp if (r["scenario"],r["law"])==(scenario,law)]
            if len(rows)!=40:
                raise RuntimeError("A hypothesis must retain both complete 20-point holdouts.")
            metrics={field:run.eq.score_groups([r|{"predicted_conversion":r[field]} for r in rows]) for field in
                     ("previous_fully_active_steady_conversion","predicted_30min_conversion","predicted_180min_conversion","predicted_30to180_sample_mean")}
            by_feed={}
            for fd in (1,2):
                v=[r for r in rows if r["feed_ratio"]==fd]
                spans=np.array([r["span_30to180_pp"] for r in v]); drift=np.array([r["change_30to180_pp"] for r in v])
                item={"window_coverage_count":sum(r["observed_in_sample_window"] for r in v),
                      "window_with_reading_allowance_coverage_count":sum(r["observed_in_window_with_half_point_reading_allowance"] for r in v),
                      "span_min_pp":float(spans.min()),"span_median_pp":float(np.median(spans)),"span_max_pp":float(spans.max()),
                      "drift_min_pp":float(drift.min()),"drift_max_pp":float(drift.max()),
                      "mean_activity_180_min":min(r["mean_activity_180min"] for r in v),
                      "mean_activity_180_max":max(r["mean_activity_180min"] for r in v),
                      "rear_minus_front_activity_30_max":max(r["last_quarter_activity_30min"]-r["first_quarter_activity_30min"] for r in v),
                      "rear_minus_front_activity_180_max":max(r["last_quarter_activity_180min"]-r["first_quarter_activity_180min"] for r in v)}
                if fd==2:
                    item["unknown_condition_trace_span_count"]=int(np.count_nonzero(spans<=trace["span_30_to_180_pp"]))
                    item["unknown_condition_trace_with_reading_allowance_count"]=int(np.count_nonzero(spans<=trace["reading_only_expanded_span_pp"]))
                by_feed[f"{fd}to1"]=item
            info[f"{scenario}__{law}"]={"scores":metrics,"history":by_feed}
    # Exactly compare each prespecified sensitivity against its own fresh-bulk history.
    primary={(r["law"],*run.eq.identity(r)):r for r in comp if r["scenario"]=="fresh_bulk"}
    sensitivities={}
    for scenario,laws in plan["scenarios"].items():
        if scenario=="fresh_bulk":
            continue
        for law in laws:
            for feed in (1,2):
                rows=[r for r in comp if r["scenario"]==scenario and r["law"]==law and r["feed_ratio"]==feed]
                fields=("predicted_10min_conversion","predicted_30min_conversion","predicted_180min_conversion","predicted_30to180_sample_mean",
                        "mean_activity_30min","mean_activity_180min")
                changes={f"max_abs_{f}_change_pp":100*max(abs(r[f]-primary[law,*run.eq.identity(r)][f]) for r in rows) for f in fields}
                sensitivities[f"{scenario}__{law}__{feed}to1"]=changes
    return comp,info,sensitivities,envelope,trace


def eq_ester(conversion,row):
    return run.eq.ester_g_h(float(conversion),row)


def selected_profiles(nums,profiles):
    rows=[]
    for d in nums:
        if not (d["scenario"]=="fresh_bulk" and d["law"] in ("acid_decay","acid_both","protected_residual_decay")
                and d["temperature_C"] in (40,60,80) and d["flow_ml_min"] in (.7,3.)):
            continue
        with np.load(profiles[d["task_id"]],allow_pickle=False) as p:
            n=len(p["cell_center_fraction"]); t=p["times_min"]
            # Volume-weighted 32-bin reduction, not invented bead-scale data.
            reduced={k:p[k].reshape(32,n//32,len(t)).mean(axis=1) for k in run.PROFILE_NAMES}
            for age in (10,30,180):
                i=int(np.flatnonzero(t==age)[0])
                for j in range(32):
                    rows.append({"task_id":d["task_id"],"law":d["law"],"feed_ratio":d["feed_ratio"],
                        "temperature_C":d["temperature_C"],"flow_ml_min":d["flow_ml_min"],"time_min":age,
                        "axial_fraction":(j+.5)/32,**{k:float(v[j,i]) for k,v in reduced.items()}})
    return rows


def reaction_zone_metrics(nums,profiles):
    rows=[]
    for d in nums:
        if d["scenario"]!="fresh_bulk":
            continue
        with np.load(profiles[d["task_id"]],allow_pickle=False) as p:
            z=p["cell_center_fraction"]; t=p["times_min"]; n=len(z)
            for i,age in enumerate(t):
                if age<10:
                    continue
                # Nonnegative reaction rate except for equilibrium roundoff.
                rate=np.maximum(0.,p["local_rate_mol_L_min"][:,i]); total=float(rate.sum())
                generated=float(rate.mean()*d["void_volume_L"])
                fraction=generated/(d["acid0_mol_L"]*d["flow_ml_min"]/1000)
                # Reaction location is meaningless once net production is
                # negligible. Do not interpret a roundoff-level centroid.
                valid=fraction>=1e-6
                rows.append({"task_id":d["task_id"],"law":d["law"],"feed_ratio":d["feed_ratio"],
                    "temperature_C":d["temperature_C"],"flow_ml_min":d["flow_ml_min"],"time_min":float(age),
                    "instantaneous_bed_generation_mol_min":generated,"generation_fraction_of_current_acid_feed":fraction,
                    "reaction_weighted_axial_fraction":float(np.dot(z,rate)/total) if valid else None,
                    "reaction_fraction_last_quarter":float(rate[-n//4:].sum()/total) if valid else None})
    return rows


def aggregate(part_paths,output):
    plan=read_json(run.PLAN)
    groups,rows,nums,profiles,manifests=predicted_ledgers(part_paths,plan)
    output.mkdir(parents=True,exist_ok=False)
    # This joined prediction file is saved before reading experimental X.
    write_rows(output/"predictions_before_comparison.csv",rows)
    write_rows(output/"numerics.csv",nums)
    comp,info,sensitivities,envelope,trace=score(groups,nums,plan)
    write_rows(output/"bed_comparison.csv",comp)
    write_rows(output/"uniform_envelope_check.csv",envelope)
    write_rows(output/"selected_axial_profiles.csv",selected_profiles(nums,profiles))
    write_rows(output/"reaction_zone.csv",reaction_zone_metrics(nums,profiles))
    profile_index={key: {"path":str(p.relative_to(ROOT)),"sha256":sha(p)} for key,p in profiles.items()}
    write_json(output/"profile_index.json",profile_index)
    limits=plan["numerics"]
    for num in nums:
        for scope in ("grid","time"):
            if not run.passes({k:num[f"{scope}_{k}"] for k in run.LIMIT_KEYS},limits):
                raise RuntimeError("A reported refinement did not pass.")
    numerical={key:max(r[key] for r in nums) for key in
        (*[f"{prefix}_{k}" for prefix,k in product(("grid","time"),run.LIMIT_KEYS)],
         "max_stage_residual","max_normalized_balance_error","independent_normalized_balance_error",
         "max_absolute_balance_mol","max_film_error_mol_L_min","max_positive_film_slope_ratio",
         "sampled_global_film_slope_ratio","sampled_global_film_slope_refinement_change","max_normalized_exposure_bound_error")}
    numerical.update({"mesh_counts":dict(sorted(Counter(str(r["fine_cells"]) for r in nums).items())),
        "step_cap_counts":dict(sorted(Counter(str(r["fine_max_step_min"]) for r in nums).items())),
        "numerical_solve_count":sum(r["numerical_solve_count"] for r in nums),
        "accepted_steps_final_histories":sum(r["accepted_steps"] for r in nums),
        "rejected_steps_final_histories":sum(r["rejected_steps"] for r in nums)})
    # The fitted acid_activation law has zero damage and recovery deficit;
    # 180 min should approach the former steady result, at the refined mesh.
    static=[r for r in comp if r["scenario"]=="fresh_bulk" and r["law"]=="acid_activation"]
    static_error=max(abs(r["predicted_180min_conversion"]-r["previous_fully_active_steady_conversion"]) for r in static)
    if static_error>limits["max_sample_conversion_change_fraction"]:
        raise RuntimeError(f"Clock-free 180-min steady replay changed by {static_error}")
    numerical["all_40_clock_free_180min_vs_steady_max_difference"]=static_error
    record={"scope":"No reactor or time-trace fitting; predeclared frozen batch laws and sensitivity scenarios; 560 qualified histories.",
        "python":manifests[0]["python"],"julia":manifests[0]["julia"],"platform":manifests[0]["platform"],
        "condition_count":40,"primary_history_count":280,"history_count":len(nums),"snapshot_count":len(rows),"actual_mass_kg":.006,
        "score_scenarios":info,"initialization_and_surface_sensitivities":sensitivities,"unknown_condition_stability_trace":trace,
        "numerical":numerical,"uniform_activity_bound_check":{"cases":len(envelope),
            "maximum_activity_violation":max(r["activity_bound_violation"] for r in envelope),
            "maximum_transient_distance_outside_steady_conversion_envelope_pp":max(r["outside_uniform_steady_interval_pp"] for r in envelope)},
        "source_sha256":run.sources()|{str(p.relative_to(ROOT)):sha(p) for p in (Path(__file__),*run.eq.OBSERVATIONS.values())},
        "partition_manifest_sha256":{str(p.relative_to(ROOT)):sha(p/"manifest.json") for p in part_paths},
        "output_sha256":{p.name:sha(p) for p in sorted(output.iterdir())}}
    write_json(output/"summary.json",record)
    return record


if __name__=="__main__":
    p=ArgumentParser(description=__doc__)
    p.add_argument("--partitions",type=Path,nargs="+",required=True)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args(); s=aggregate([x.resolve() for x in args.partitions],args.output.resolve())
    for name,d in s["score_scenarios"].items():
        m=d["scores"]["predicted_30to180_sample_mean"]
        print(name,"RMSE 1:1 / 2:1",round(m["1to1"]["rmse_percentage_points"],2),round(m["2to1"]["rmse_percentage_points"],2),flush=True)
