"""Freeze all short-reactor decisions before making a long-reactor prediction.

The measured blank already includes the valve pulse, tubing and biosensor.
Extra delay laws are diagnostic alternatives, calibrated on A only; they are
not independently measured hardware or pore-access parameters.
"""
import argparse
import math
from pathlib import Path

import numpy as np

from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT
from . import loop, native
from . import storage_audit, source_audit
from .handoff import matched_dt

DT = .05  # s; integration/binning cadence, not the marker integration step
BACKGROUNDS = (0, 2, 5)  # common approximate figure-reading sensitivity
MODELS = ("measured_blank", "extra_mixed_volume", "common_rare_delay",
          "solid_volume_rare_delay")


def load_trace(base, reactor, depth, dt, factor=1.):
    plan = read_json(base/"plan.json")
    n = plan["protocol"]["numerics"]
    tag = f"dt{dt:.8g}_D{factor:.4g}_n{n['particles']}_seed{n['particle_seed']}"
    directory = base/f"{reactor}_H{depth}"
    path = directory/(tag+".f64")
    record = read_json(directory/(tag+".json"))
    carrier = read_json(directory/"handoff.json")
    hydraulic = read_json(directory/"hydraulic.json")
    if (record["marker_sha256"] != sha(path) or
            record["handoff_sha256"] != sha(directory/"handoff.json") or
            carrier["field_sha256"] != hydraulic["field_sha256"] or
            (record["reactor"],record["depth_cells"],record["dt_s"],record["factor"]) !=
            (reactor,depth,dt,factor)):
        raise ValueError("Tracer record is not from the specified carrier and settings")
    times = np.fromfile(path,dtype="<f8")
    if len(times)!=n["particles"] or (np.isfinite(times).mean()!=record["recovery"]):
        raise ValueError("Exit-time mass ledger changed")
    return times, dict(tracer=record, handoff=carrier, hydraulic=hydraulic)


def numerical_ok(record, protocol):
    h, f, t = record["hydraulic"]["native"], record["handoff"], record["tracer"]
    n = protocol["numerics"]
    return bool(h["converged"] and h["collision"]=="TRT" and
        h["max_Mach"]<=n["max_Mach"] and h["max_density_departure"]<=n["max_density_departure"] and
        abs(h["flow_m3_s"]/record["hydraulic"]["target_flow_m3_s"]-1)<=n["max_native_flow_relative_error"] and
        f["carrier_qualified"] and t["carrier_qualified"] and
        t["recovery"]>=n["min_recovery"] and t["numerical_mass_recovery_passed"] and
        abs(t["mean_vs_volume_relative"])<=n["mean_refinement_relative"])


def refinement(coarse, fine, protocol):
    n=protocol["numerics"]
    a,b=coarse["tracer"],fine["tracer"]
    m=b["conditional_mean_s"]/a["conditional_mean_s"]-1
    v=b["conditional_variance_s2"]/a["conditional_variance_s2"]-1
    return dict(mean_relative=m,variance_relative=v,
                moments_pass=abs(m)<=n["mean_refinement_relative"] and abs(v)<=n["variance_refinement_relative"],
                qualified_carriers=numerical_ok(coarse,protocol) and numerical_ok(fine,protocol))


def fit_A_only(score, plan):
    """Translate A's *unexplained* first two moments into physical diagnostics."""
    m=score["incremental_observed"]["mean_s"]-score["incremental_predicted"]["mean_s"]
    v=score["incremental_observed"]["variance_s2"]-score["incremental_predicted"]["variance_s2"]
    q=plan["protocol"]["fluid"]["flow_ul_min"]/60
    # For a fixed-volume mixed connection, its variance must be m^2.
    # For exponential pauses, v = 2 * lambda * tau^2 and m = lambda * tau.
    rare=dict(available=False)
    if m>0 and v>0:
        tau=v/(2*m)
        rate=m/tau
        if rate<=1:  # a rare-delay diagnostic, not an unrestricted fit
            rare=dict(available=True,rate_A=rate,mean_event_s=tau,
                      delayed_fraction_A=-math.expm1(-rate),
                      moment_constraint="Two A incremental moments by construction; shape is not fitted.")
    return dict(unexplained_mean_A_s=m,unexplained_variance_A_s2=v,
                implied_positive_hold_up_A_ul=max(m,0.)*q,
                extra_mixed_volume_mean_s=max(m,0.),rare=rare,
                calibration_reactor="A",
                note="Extra delays on top of measured common loop. Physical origin unqualified.")


def solid_volume_ratio(plan, reactor):
    def volume(r):
        g=plan["reactors"][r]
        return float(np.prod(g["extent_m"])*(1-g["geometry"]["continuous_porosity"]))
    return volume(reactor)/volume("A")


def predict_curves(times, reactor, plan, parameters):
    """No measured A or C curve is accessed here; only the frozen blank."""
    scale=solid_volume_ratio(plan,reactor)
    r=parameters["rare"]
    kernels={"measured_blank":np.ones(1),
             "extra_mixed_volume":loop.exponential_connection(parameters["extra_mixed_volume_mean_s"],DT)}
    if r["available"]:
        kernels.update(common_rare_delay=loop.rare_delay(r["rate_A"],r["mean_event_s"],DT),
                       solid_volume_rare_delay=loop.rare_delay(scale*r["rate_A"],r["mean_event_s"],DT))
    bed=loop.tracer_kernel(times,DT)
    arrays={}
    for label,k in kernels.items():
        combined=np.convolve(bed,k)
        for bg in BACKGROUNDS:
            _,blank=loop.measured("blank",DT,bg)
            _,pred=loop.convolve(blank,combined,DT)
            arrays[f"{label}_bg{bg}"]=pred
    return arrays


def passes(score,protocol):
    compare=protocol["comparison"]
    keys={"whole_mean_relative":"whole_loop_mean_relative",
          "whole_variance_relative":"whole_loop_variance_relative",
          "whole_CDF_max":"whole_loop_cdf_max",
          "incremental_mean_relative":"incremental_mean_relative",
          "incremental_variance_relative":"incremental_variance_relative"}
    return {metric:score[metric] is not None and abs(score[metric])<=compare[limit]
            for metric,limit in keys.items()}


def evaluate_curves(curves, reactor, protocol):
    report={}
    for name,c in curves.items():
        t=np.arange(len(c))*DT
        score=loop.compare_prediction(t,c,reactor,DT,int(name.split("_bg")[1]))
        check=passes(score,protocol)
        report[name]=dict(**score,gates=check,all_gates=all(check.values()))
    return report


def table1(reactor, scores):
    row=next(r for r in native.source()["table1"]["rows"] if r[0]==reactor)
    return dict(reported_mean_s=60*row[6],reported_variance_s2=3600*row[7],reported_skewness=row[8],
                fig3b_incremental=scores["measured_blank_bg0"]["incremental_observed"],
                predicted_by_model={k:v["incremental_predicted"] for k,v in scores.items() if k.endswith("_bg0")},
                note="Printed Table 1 and moments of digitized Fig. 3b are not interchangeable.")


def freeze_A(base, depth, dt, diagnostic=False):
    """Both depth and dt are screened on A; no holdout measurement is read."""
    target=base/"frozen_A.json"
    if target.exists():
        raise ValueError("Short-reactor calibration already frozen")
    plan=read_json(base/"plan.json")
    p=plan["protocol"]
    if sha(native.PROTOCOL)!=plan["protocol_sha256"] or sha(native.SOURCE)!=plan["source_sha256"]:
        raise ValueError("Frozen inputs changed")
    expected_dt=matched_dt(p,max(p["geometry"]["depth_grids"]))
    if dt!=expected_dt:
        raise ValueError("Grid tests must use the matched physical tracer time step")
    times,selected=load_trace(base,"A",depth,dt)
    grids={}
    for h in p["geometry"]["depth_grids"]:
        try:
            _,g=load_trace(base,"A",h,dt)
        except FileNotFoundError:
            continue
        grids[str(h)]=g
    changes={f"H{a}_to_H{b}":refinement(grids[str(a)],grids[str(b)],p)
             for a,b in zip(p["geometry"]["depth_grids"][:-1],p["geometry"]["depth_grids"][1:])
             if str(a) in grids and str(b) in grids}
    # H12 is diagnostic if it fails face correction; H16 and H20 are the
    # intended independently qualified pair. H16 is sufficient if it agrees
    # with H20 at the prefixed 5% mean / 10% variance limits.
    final_grid=changes.get("H16_to_H20",{})
    grid_pass=bool(final_grid.get("moments_pass") and final_grid.get("qualified_carriers") and depth in (16,20))
    try:
        _,fine=load_trace(base,"A",depth,dt/2)
        time_screen=refinement(selected,fine,p)
    except FileNotFoundError:
        fine=None
        time_screen=None
    step_pass=bool(time_screen and time_screen["moments_pass"] and time_screen["qualified_carriers"])
    qualified=numerical_ok(selected,p) and grid_pass and step_pass
    if not qualified and not diagnostic:
        raise ValueError("A numerical screens are incomplete or failed; explicit diagnostic mode needed")
    _,_,baseline=loop.compare_times(times,"A",DT)
    params=fit_A_only(baseline,plan)
    apred=predict_curves(times,"A",plan,params)
    scores=evaluate_curves(apred,"A",p)
    sensitivity={}
    for factor in p["fluid"]["diffusivity_sensitivity_factors"]:
        try:
            t,s=load_trace(base,"A",depth,dt,factor)
        except FileNotFoundError:
            continue
        sensitivity[str(factor)]=dict(record=s,baseline=loop.compare_times(t,"A",DT)[2])
    np.savez_compressed(base/"A_predictions.npz",**apred)
    record=dict(reaction_enabled=False,calibration_reactor="A",holdout_reactor="C",
        protocol_sha256=sha(native.PROTOCOL),plan_sha256=sha(base/"plan.json"),
        figure_readings_sha256=sha(loop.MARKERS),comparison_code_sha256=sha(Path(__file__)),
        loop_code_sha256=sha(Path(loop.__file__)),
        selected_depth_cells=depth,selected_dt_s=dt,grid_records=grids,grid_changes=changes,
        grid_pass=grid_pass,fine_step_record=fine,step_change=time_screen,time_step_pass=step_pass,
        numerical_qualification=qualified,selected=selected,diffusivity_sensitivity=sensitivity,
        model_parameters=params,background_reading_pixels=list(BACKGROUNDS),comparison_dt_s=DT,
        dimensionless_source_audit=source_audit.audit(),
        independent_storage_parameters=storage_audit.parameters(),
        auxiliary_code_sha256={str(Path(module.__file__).relative_to(ROOT)):sha(Path(module.__file__))
                               for module in (source_audit,storage_audit)},
        solid_volume_ratio_C=solid_volume_ratio(plan,"C"),A_scores=scores,A_table1=table1("A",scores),
        A_predictions_sha256=sha(base/"A_predictions.npz"),
        independent_A_pass=qualified and scores["measured_blank_bg0"]["all_gates"])
    write_json(target,record)
    print("Frozen A:","numerical pass",qualified,"independent A pass",record["independent_A_pass"],
          "A-calibrated diagnostic parameters",params,flush=True)
    return record


def verify_freeze(base):
    frozen=read_json(base/"frozen_A.json")
    for key,path in (("plan_sha256",base/"plan.json"),("protocol_sha256",native.PROTOCOL),
                     ("figure_readings_sha256",loop.MARKERS),("comparison_code_sha256",Path(__file__)),
                     ("loop_code_sha256",Path(loop.__file__)),("A_predictions_sha256",base/"A_predictions.npz")):
        if sha(path)!=frozen[key]:
            raise ValueError(f"Frozen A input changed: {path}")
    for name,digest in frozen["auxiliary_code_sha256"].items():
        if sha(ROOT/name)!=digest:
            raise ValueError(f"Frozen auxiliary model changed: {name}")
    if (sha(source_audit.FIG3A)!=frozen["dimensionless_source_audit"]["fig3a_sha256"] or
        sha(storage_audit.EVIDENCE)!=frozen["independent_storage_parameters"]["evidence_sha256"]):
        raise ValueError("Frozen independent evidence changed")
    return frozen


def predict_C(base):
    """Write the complete long prediction to disk *before* scoring against C."""
    path=base/"C_prediction_before_comparison.npz"
    if path.exists():
        raise ValueError("Long prediction already exists")
    frozen=verify_freeze(base)
    plan=read_json(base/"plan.json")
    t,c=load_trace(base,"C",frozen["selected_depth_cells"],frozen["selected_dt_s"])
    arrays=predict_curves(t,"C",plan,frozen["model_parameters"])
    np.savez_compressed(path,**arrays)
    record=dict(reaction_enabled=False,holdout_reactor="C",calibration_reactor="A",
        frozen_A_sha256=sha(base/"frozen_A.json"),predictions_sha256=sha(path),
        carrier=c,numerical_ok=numerical_ok(c,plan["protocol"]),comparison_dt_s=DT,
        conditional_diagnostic=not frozen["independent_A_pass"],
        prediction_curves=list(arrays))
    write_json(base/"C_prediction_manifest.json",record)
    print("Saved C predictions before comparison; numerical pass",record["numerical_ok"],flush=True)
    return record


def score_C(base):
    frozen=verify_freeze(base)
    plan=read_json(base/"plan.json")
    record=read_json(base/"C_prediction_manifest.json")
    path=base/"C_prediction_before_comparison.npz"
    if sha(path)!=record["predictions_sha256"] or sha(base/"frozen_A.json")!=record["frozen_A_sha256"]:
        raise ValueError("Prediction changed after A was frozen")
    with np.load(path) as curves:
        scores=evaluate_curves(curves,"C",plan["protocol"])
    holdout=bool(record["numerical_ok"] and scores["measured_blank_bg0"]["all_gates"])
    summary=dict(reaction_enabled=False,A=frozen,C_manifest=record,C_scores=scores,
                 C_table1=table1("C",scores),independent_C_pass=holdout,
                 independent_storage_bounds={r:storage_audit.check(plan,r,s["measured_blank_bg0"],
                    frozen["independent_storage_parameters"]) for r,s in (("A",frozen["A_scores"]),("C",scores))},
                 reaction_gate_pass=frozen["independent_A_pass"] and holdout,
                 note="Calibration alternatives are diagnostics. A fit to two moments does not independently qualify the hardware or bed physics.")
    write_json(base/"summary.json",summary)
    print("Independent C pass",holdout,"reaction gate pass",summary["reaction_gate_pass"],flush=True)
    return summary


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("base",type=Path)
    parser.add_argument("action",choices=("freeze","predict","score"))
    parser.add_argument("--depth",type=int,choices=(16,20),default=16)
    parser.add_argument("--dt",type=float,default=.00125)
    parser.add_argument("--diagnostic",action="store_true")
    a=parser.parse_args()
    path=a.base.resolve()
    if a.action=="freeze":
        freeze_A(path,a.depth,a.dt,a.diagnostic)
    elif a.action=="predict":
        predict_C(path)
    else:
        score_C(path)


if __name__=="__main__":
    main()
