"""Predict time-on-stream catalyst exposure: no reactor conversion input to the solver."""

from argparse import ArgumentParser
from datetime import datetime, timezone
from itertools import product
from pathlib import Path
from time import perf_counter
import platform

import numpy as np

from examples.biocatalysis.pentyl_propanoate.run import read_rows, write_rows
from examples.biocatalysis.pentyl_propanoate_equimolar import run as eq
from examples.biocatalysis.pentyl_propanoate_equimolar.check import check as check_previous
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT

HERE = Path(__file__).resolve().parent
PLAN = HERE / "protocol.json"
OUTPUT = ROOT / "results/biocatalysis/pentyl_propanoate_history"
PREVIOUS = eq.HERE / "reference"
CONDITIONING = eq.prior.base.HERE.parent / "pentyl_propanoate_conditioning/inputs.json"
ARRAYS = ("times_min", "cell_center_fraction", "bulk_progress", "surface_progress",
          "cumulative_hazard", "activity", "local_rate_mol_L_min", "local_hazard_min_inv",
          "outlet_conversion", "generated_mol", "effluent_product_mol", "backdiffused_product_mol",
          "liquid_product_mol")
SCALARS = ("initial_liquid_product_mol", "catalyst_mass_kg", "catalyst_mass_per_cell_kg",
           "catalyst_g_L", "residence_min", "equilibrium_conversion", "inlet_hazard_min_inv",
           "equilibrium_hazard_min_inv", "min_step_min", "max_step_min", "max_stage_residual",
           "max_normalized_balance_error", "max_absolute_balance_mol", "max_normalized_exposure_bound_error",
           "max_film_error_mol_L_min", "max_positive_film_slope_ratio")
INTS = ("accepted_steps", "rejected_steps", "max_newton_iterations")
PROFILE_NAMES = ("bulk_progress", "surface_progress", "cumulative_hazard", "activity",
                 "local_rate_mol_L_min", "local_hazard_min_inv")
LEDGERS = ("generated_mol", "effluent_product_mol", "backdiffused_product_mol", "liquid_product_mol")
LIMIT_KEYS = {
    "sample_conversion_change_fraction": "max_sample_conversion_change_fraction",
    "startup_conversion_change_fraction": "max_startup_conversion_change_fraction",
    "mean_activity_change_fraction": "max_mean_activity_change_fraction",
    "profile_activity_L1_change_fraction": "max_profile_activity_L1_change_fraction",
    "cumulative_equivalent_change_per_fed_mole": "max_cumulative_equivalent_change_per_fed_mole",
}


def load_runtime():
    jl = eq.load_runtime()
    if not jl.seval("isdefined(@__MODULE__, :SpatialCatalystHistory)"):
        jl.include(str(HERE / "model.jl"))
    return jl


def sources():
    files=(PLAN,HERE/"model.jl",HERE/"run.py",eq.FROZEN,eq.EXTRA,eq.PUBLISHED.MODEL_INPUTS,
           eq.PROPERTIES,eq.PUBLISHED.MODEL,eq.HERE/"run.py",eq.HERE/"model.jl",
           eq.prior.HERE/"run.py",eq.prior.HERE/"model.jl",eq.prior.base.HERE/"model.jl",
           eq.prior.loading.HERE/"model.jl",eq.prior.loading.film.HERE/"model.jl",
           PREVIOUS/"summary.json",PREVIOUS/"bed_comparison.csv",PREVIOUS/"clock_comparison.csv",
           CONDITIONING)
    return {str(p.relative_to(ROOT)):sha(p) for p in files}


def prepare(jl,plan):
    cases=eq.conditions(jl,read_json(eq.PUBLISHED.MODEL_INPUTS),read_json(eq.EXTRA),read_json(eq.PROPERTIES))
    frozen=read_json(eq.FROZEN)["batch"]["all_data"]
    if len(frozen["training_assays"])!=16 or plan["scenarios"]["fresh_bulk"]!=list(eq.prior.CLOCK):
        raise ValueError("The frozen batch calibration boundary changed.")
    laws={law:eq.prior.parameters_from_record(jl.CompositionKinetics,frozen["fit"][law])
          for law in plan["scenarios"]["fresh_bulk"]}
    tasks=[]
    for scenario, names in plan["scenarios"].items():
        for law,c in product(names,cases):
            feed,T,Q=eq.identity(c)
            K=eq.prior.equilibrium_value(jl.CompositionKinetics,frozen["equilibrium_training"],feed,T)
            task_id=f"{scenario}__{law}__r{feed}_T{T}_Q{Q:g}".replace(".","p")
            if c["actual_mass_kg"]!=plan["scientific_boundary"]["actual_mass_kg"]:
                raise ValueError("The weighed catalyst inventory changed.")
            tasks.append({"task_index":len(tasks),"task_id":task_id,"scenario":scenario,"law":law,
                          "apparent_Kx":K,**c})
    if len(tasks)!=560 or len({c["task_id"] for c in tasks})!=560:
        raise ValueError("Expected all 560 prespecified condition histories.")
    return tasks,laws


def conditions_flags(task):
    return {"surface_hazard":task["scenario"]=="fresh_surface", "recovered":task["scenario"]=="recovered_bulk",
            "equilibrium_prefill":task["scenario"]=="equilibrium_prefill_bulk"}


def one_history(jl,p,c,plan,n,dt):
    result=jl.SpatialCatalystHistory.simulate(p,c["temperature_C"],c["feed_ratio"],c["apparent_Kx"],
        np.asarray(plan["times_min"],dtype=float),catalyst_kg=c["actual_mass_kg"],liquid_L=c["void_volume_L"],
        flow_L_min=c["flow_ml_min"]/1000,acid0=c["acid0_mol_L"],alcohol0=c["alcohol0_mol_L"],
        capacity=c["capacity_mol_L_min"],eta_internal=c["eta_internal"],peclet=c["peclet"],
        cells=n,max_step_min=dt,**conditions_flags(c))
    data={name:np.array(getattr(result,name),dtype=np.float64,copy=True) for name in ARRAYS}
    data.update({name:float(getattr(result,name)) for name in SCALARS})
    data.update({name:int(getattr(result,name)) for name in INTS})
    data.update({"cells":n,"requested_max_step_min":dt})
    validate_result(data,c,plan)
    return data


def validate_result(r,c,plan):
    n=r["cells"]; times=np.asarray(plan["times_min"],dtype=float); nt=len(times)
    if not (np.array_equal(r["times_min"],times) and
            np.array_equal(r["cell_center_fraction"],(np.arange(n)+.5)/n) and
            all(r[name].shape==(n,nt) for name in PROFILE_NAMES) and
            all(r[name].shape==(nt,) for name in (*LEDGERS,"outlet_conversion")) and
            np.array_equal(r["outlet_conversion"],r["bulk_progress"][-1,:])):
        raise RuntimeError(f"Returned array contract failed: {c['task_id']}")
    if not all(np.all(np.isfinite(r[key])) for key in (*ARRAYS,*SCALARS)):
        raise RuntimeError("Returned nonfinite transient values.")
    x,y,a,H=r["bulk_progress"],r["surface_progress"],r["activity"],r["cumulative_hazard"]
    if not (np.all((x>=0)&(x<=r["equilibrium_conversion"])) and np.all(y>=x-1e-12) and
            np.all(y<=r["equilibrium_conversion"]+1e-12) and np.all((a>=0)&(a<=1)) and
            np.all(H>=0) and np.all(np.diff(H,axis=1)>=-1e-12)):
        raise RuntimeError("Liquid or stationary catalyst state violated physical bounds.")
    if not (r["catalyst_mass_kg"]==c["actual_mass_kg"]==.006 and
            abs(r["catalyst_mass_per_cell_kg"]*n-.006)<1e-14 and
            abs(r["catalyst_g_L"]-c["catalyst_g_L"])<1e-10 and
            abs(r["residence_min"]-c["residence_min"])<1e-10):
        raise RuntimeError("The six-gram bed inventory or liquid residence changed.")
    # Independently reconstruct the dimensional cumulative balance from
    # the saved physical quantities and cell liquid inventories.
    I=c["acid0_mol_L"]*c["void_volume_L"]*np.mean(x,axis=0)
    fed=c["acid0_mol_L"]*c["flow_ml_min"]/1000*np.maximum.reduce((times,np.full(nt,r["residence_min"]),np.ones(nt)))
    balance=np.abs(I-r["initial_liquid_product_mol"]+r["effluent_product_mol"]+
                   r["backdiffused_product_mol"]-r["generated_mol"])/fed
    r["independent_normalized_balance_error"]=float(np.max(balance))
    if not np.allclose(I,r["liquid_product_mol"],rtol=1e-12,atol=1e-14):
        raise RuntimeError("Liquid inventory unit reconstruction failed.")
    limits=plan["numerics"]
    checks={"max_stage_residual":"max_stage_residual_fraction",
            "max_normalized_balance_error":"max_normalized_cumulative_balance_error",
            "independent_normalized_balance_error":"max_normalized_cumulative_balance_error",
            "max_film_error_mol_L_min":"max_film_residual_mol_L_min",
            "max_positive_film_slope_ratio":"max_sampled_film_source_slope_over_capacity"}
    if not all(r[name]<=limits[bound] for name,bound in checks.items()):
        raise RuntimeError(f"Conservation/film numerical limit failed: {c['task_id']}")
    if r["max_normalized_exposure_bound_error"]>1e-10:
        raise RuntimeError("Stationary exposure history exceeded its physical range.")


def refinement(coarse,fine,c):
    t=fine["times_min"]
    dx=np.abs(coarse["outlet_conversion"]-fine["outlet_conversion"])
    a,b=coarse["activity"],fine["activity"]
    if len(a)!=len(b):
        if len(b)!=2*len(a):
            raise ValueError("Compare the same grid or a matched factor-two refinement.")
        b=b.reshape(len(a),2,len(t)).mean(axis=1)
    da=np.abs(a-b)
    denom=c["acid0_mol_L"]*c["flow_ml_min"]/1000*t[1:]
    return {"sample_conversion_change_fraction":float(np.max(dx[t>=10])),
            "startup_conversion_change_fraction":float(np.max(dx[t<10])),
            "mean_activity_change_fraction":float(np.max(np.abs(a.mean(axis=0)-b.mean(axis=0)))),
            "profile_activity_L1_change_fraction":float(np.max(da.mean(axis=0))),
            "cumulative_equivalent_change_per_fed_mole":max(float(np.max(np.abs(coarse[k][1:]-fine[k][1:])/denom))
                                                             for k in LEDGERS)}


def passes(changes,settings):
    return all(changes[k]<=settings[limit] for k,limit in LIMIT_KEYS.items())


class NumericalLimit(RuntimeError):
    pass


def solve_refined(jl,p,c,plan):
    limits=plan["numerics"]
    n=limits["initial_cells"][-1]; dt=limits["initial_max_steps_min"][-1]
    cache={}
    def get(n,dt):
        key=(n,dt)
        if key not in cache:
            cache[key]=one_history(jl,p,c,plan,n,dt)
        return cache[key]
    # Sample the local film root over the full physical progress domain
    # and all stage sizes up to the *coarser* time step. If the sufficient
    # monotonicity margin is too small, refine before solving the stage.
    while True:
        s=jl.SpatialCatalystHistory.film_slope_screen(p,c["temperature_C"],c["feed_ratio"],c["apparent_Kx"],
            catalyst_g_L=c["catalyst_g_L"],eta_internal=c["eta_internal"],capacity=c["capacity_mol_L_min"],
            max_step_min=2*dt,surface_hazard=conditions_flags(c)["surface_hazard"])
        screen=float(s.max_sampled_slope_ratio); screen_change=float(s.refinement_change)
        if screen<limits["max_sampled_film_source_slope_over_capacity"]:
            break
        dt/=2
        if dt<limits["minimum_max_step_min"]:
            raise NumericalLimit(f"Global sampled film-slope criterion failed: {c['task_id']}, {screen}")
    while True:
        fine=get(n,dt)
        time=refinement(get(n,2*dt),fine,c)
        if not passes(time,limits):
            dt/=2
            if dt<limits["minimum_max_step_min"]:
                raise NumericalLimit(f"Time refinement cap reached: {c['task_id']}: {time}")
            continue
        grid=refinement(get(n//2,dt),fine,c)
        if not passes(grid,limits):
            n*=2
            if n>limits["max_cells"]:
                raise NumericalLimit(f"Spatial refinement cap reached: {c['task_id']}: {grid}")
            continue
        break
    summary={"fine_cells":n,"fine_max_step_min":dt,"coarse_cells":n//2,"coarse_max_step_min":dt*2,
             "numerical_solve_count":len(cache),"sampled_global_film_slope_ratio":screen,
             "sampled_global_film_slope_refinement_change":screen_change,
             **{f"time_{k}":v for k,v in time.items()},**{f"grid_{k}":v for k,v in grid.items()},
             **{k:fine[k] for k in (*SCALARS,*INTS,"independent_normalized_balance_error")}}
    return fine,summary


def case_head(c):
    return {k:c[k] for k in ("task_index","task_id","scenario","law","feed_ratio","temperature_C","flow_ml_min",
                            "actual_mass_kg","void_volume_L","residence_min","acid0_mol_L","alcohol0_mol_L",
                            "catalyst_g_L","eta_internal","peclet","capacity_mol_L_min","apparent_Kx")}


def snapshots(r,c):
    t=r["times_min"]; n=r["cells"]; a=r["activity"]; H=r["cumulative_hazard"]
    return [case_head(c)|{"time_min":float(age),"predicted_conversion":float(r["outlet_conversion"][i]),
        "nominal_ester_outlet_g_h":eq.ester_g_h(float(r["outlet_conversion"][i]),c),
        "equilibrium_conversion":r["equilibrium_conversion"],"mean_liquid_progress":float(np.mean(r["bulk_progress"][:,i])),
        "mean_activity":float(np.mean(a[:,i])),"minimum_activity":float(np.min(a[:,i])),
        "maximum_activity":float(np.max(a[:,i])),"first_quarter_activity":float(a[:n//4,i].mean()),
        "last_quarter_activity":float(a[-n//4:,i].mean()),"mean_cumulative_hazard":float(H[:,i].mean()),
        "minimum_cumulative_hazard":float(H[:,i].min()),"maximum_cumulative_hazard":float(H[:,i].max()),
        "inlet_cell_hazard":float(H[0,i]),"outlet_cell_hazard":float(H[-1,i]),
        **{k:float(r[k][i]) for k in LEDGERS}}
        for i,age in enumerate(t)]


def prediction_partition(output,part,parts):
    check_previous()
    plan=read_json(PLAN); jl=load_runtime(); tasks,laws=prepare(jl,plan)
    if not 0<=part<parts:
        raise ValueError("Invalid deterministic partition.")
    chosen=[c for c in tasks if c["task_index"]%parts==part]
    if not chosen:
        raise ValueError("Empty partition.")
    output.mkdir(parents=True,exist_ok=False)
    (output/"profiles").mkdir(); (output/"cases").mkdir()
    input_hashes=sources(); rows=[]; numerical=[]; failure=[]; profile_hashes={}
    start=perf_counter()
    for index,c in enumerate(chosen,1):
        try:
            fine,diag=solve_refined(jl,laws[c["law"]],c,plan)
        except NumericalLimit as e:
            failure.append({"task_id":c["task_id"],"reason":str(e)})
            write_json(output/"failures.json",failure)
            print(f"Unqualified numerical case: {e}",flush=True)
            continue
        path=output/"profiles"/(c["task_id"]+".npz")
        np.savez_compressed(path,**{k:fine[k] for k in ARRAYS})
        profile_hashes[str(path.relative_to(output))]=sha(path)
        ss=snapshots(fine,c); dd=case_head(c)|diag
        # Keep each completed case immutable on disk in addition to the
        # final aggregate ledgers, so a numerical interruption is traceable.
        write_json(output/"cases"/(c["task_id"]+".json"),{"snapshots":ss,"numerical":dd,
            "profile_sha256":profile_hashes[str(path.relative_to(output))]})
        rows.extend(ss); numerical.append(dd)
        if index==1 or index%5==0 or index==len(chosen):
            print(f"Part {part+1}/{parts}: {index}/{len(chosen)} {c['task_id']} N{diag['fine_cells']} "
                  f"dt≤{diag['fine_max_step_min']:g} min ({perf_counter()-start:.1f} s)",flush=True)
    if not rows:
        raise RuntimeError("No qualified histories in this partition.")
    write_rows(output/"predictions.csv",rows); write_rows(output/"numerics.csv",numerical)
    manifest={"part":part,"parts":parts,"expected_count":len(chosen),"completed_count":len(numerical),
              "failed_cases":failure,"python":platform.python_version(),"julia":str(jl.seval("VERSION")),
              "platform":platform.platform(),"sources":input_hashes,"profile_sha256":profile_hashes,
              "ledger_sha256":{name:sha(output/name) for name in ("predictions.csv","numerics.csv")}}
    write_json(output/"manifest.json",manifest)
    return manifest


if __name__=="__main__":
    parser=ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,default=OUTPUT/datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ"))
    parser.add_argument("--part",type=int,default=0)
    parser.add_argument("--parts",type=int,default=1)
    args=parser.parse_args()
    manifest=prediction_partition(args.output.resolve(),args.part,args.parts)
    print(f"Saved {manifest['completed_count']} qualified histories at {args.output}",flush=True)
    if manifest["failed_cases"]:
        raise SystemExit(1)
