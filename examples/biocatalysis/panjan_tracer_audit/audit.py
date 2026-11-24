"""Audit visible pulse evidence, then test prior frozen tracer transfer laws.

All concentrations are mM, all time arrays s, liquid volumes µL. Wide
traces are partial integrals only. No reaction or new OpenLB simulation.
"""
import argparse
import csv
from datetime import datetime,timezone
import math
from pathlib import Path
import platform
import sys
import tomllib
from uuid import uuid4

import numpy as np

from pore2process_demo.io import read_json,sha,write_json
from pore2process_demo.paths import ROOT
from examples.biocatalysis.strnisa_tracer import loop

HERE=Path(__file__).resolve().parent
SOURCE=ROOT/"literature/sources/panjan2018/evidence.toml"
FIGURE=SOURCE.with_name("rtd_fig25.toml")
ARTICLE=ROOT/"literature/sources/strnisa2018/evidence.toml"
PRIOR=ROOT/"examples/biocatalysis/strnisa_tracer/reference/summary.json"
PROTOCOL=HERE/"protocol.toml"
REF=HERE/"reference"
FIGURES=HERE/"figures"
NARROW=("blank","R1","R2")
WIDE=("R3","R4")


def toml(path):
    with path.open("rb") as f:
        return tomllib.load(f)


def marker_values(name,record=None):
    r=toml(FIGURE) if record is None else record
    p=np.array(r[name]["points"],dtype=float)
    unit=5/(r["axis_conc_0_mM_px"]-r["axis_conc_5_mM_px"])
    return p[:,0].astype(int),np.maximum(0.,(r["axis_conc_0_mM_px"]-p[:,1])*unit),unit


def pulse(name,bg=0,close=True,dt=.05,record=None):
    """Concentration vs time for the *visible* window; no wide-tail closure."""
    r=toml(FIGURE) if record is None else record
    k,c,unit=marker_values(name,r)
    t=k*r["time_step_s"]
    c=np.maximum(0.,c-bg*unit)
    if name in NARROW:
        t=np.r_[0.,t]; c=np.r_[0.,c]
        if close:
            t=np.r_[t,t[-1]+toml(PROTOCOL)["readings"]["narrow_tail_closure_s"]]
            c=np.r_[c,0.]
    elif name not in WIDE:
        raise ValueError(name)
    grid=t[0]+np.arange(round((t[-1]-t[0])/dt)+1)*dt
    return grid,np.interp(grid,t,c)


def duplicate_audit():
    p=toml(FIGURE); a=toml(loop.MARKERS)
    unit=4/(a["axis_conc_0_mM_px"]-a["axis_conc_4_mM_px"])
    out={}
    for name,other,j in (("blank","blank",1),("R1","A",2),("R2","C",3)):
        keys,c,u=marker_values(name,p)
        first={int(k):float(v) for k,v in zip(keys,c)}
        second={int(row[0]):(a["axis_conc_0_mM_px"]-row[j])*unit for row in a["markers"]}
        matched=sorted(set(first)&set(second))
        core=[k for k in matched if k not in p[name]["ambiguous_indices"] and k not in a["overlap_indices"] and
              min(first[k],second[k])>toml(PROTOCOL)["readings"]["duplicate_core_min_mM"]]
        f=np.array([first[k] for k in core]); s=np.array([second[k] for k in core]); d=f-s
        allowance=p["distinct_reading_px"]*u+a["distinct_marker_reading_px"]*unit
        alld=np.array([first[k]-second[k] for k in matched])
        out[name]=dict(article_series=other,count_core=len(core),match_indices=core,
            difference_mM_core=d.tolist(),rmse_mM=float(np.sqrt(np.mean(d**2))),
            max_error_mM=float(np.max(abs(d))),correlation=float(np.corrcoef(f,s)[0,1]),
            all_matched_count=len(matched),all_matched_rmse_mM=float(np.sqrt(np.mean(alld**2))),
            all_matched_max_error_mM=float(np.max(abs(alld))),
            combined_distinct_pixel_allowance_mM=allowance,
            all_core_within_reading_allowance=bool(np.max(abs(d))<=allowance),
            interpretation="No time shift, amplitude scaling or normalization. Agreement suggests reuse, not proof of independent trials.")
    return out


def last_rising_crossing(t,c,threshold):
    """Last upward crossing before the principal peak, excluding early baseline noise."""
    t=np.asarray(t); c=np.asarray(c); stop=int(np.argmax(c))
    found=np.flatnonzero((c[:stop]<threshold)&(c[1:stop+1]>=threshold))
    if len(found)==0:
        return None
    i=int(found[-1])
    return float(t[i]+(t[i+1]-t[i])*(threshold-c[i])/(c[i+1]-c[i]))


def moment_scenarios():
    cfg=toml(PROTOCOL)["readings"]
    q=toml(SOURCE)["tracer"]["flow_ul_min"]
    out={}
    for bg in cfg["background_pixels"]:
        curves={k:pulse(k,bg,dt=cfg["quadrature_dt_s"]) for k in NARROW}
        raw={k:loop.moments(*v,q) for k,v in curves.items()}
        unclosed={k:loop.moments(*pulse(k,bg,close=False,dt=cfg["quadrature_dt_s"]),q) for k in NARROW}
        def increment(m):
            return {k:{key:m[k][key]-m["blank"][key] for key in ("mean_s","variance_s2")} for k in ("R1","R2")}
        truncations={}
        for limit in cfg["tail_cutoffs_s"]:
            trunc={}
            for name,(t,c) in curves.items():
                keep=t<=limit
                trunc[name]=loop.moments(t[keep],c[keep],q)
            truncations[str(limit)]=dict(moments=trunc,increment=increment(trunc),
                 area_fraction={k:trunc[k]["observed_glucose_umol"]/raw[k]["observed_glucose_umol"] for k in NARROW})
        wide={}
        for name in WIDE:
            t,c=pulse(name,bg,close=False,dt=cfg["quadrature_dt_s"])
            wide[name]=dict(window_s=[float(t[0]),float(t[-1])],first_mM=float(c[0]),last_mM=float(c[-1]),
                peak_mM=float(c.max()),partial_glucose_umol=float(np.trapezoid(c,t)*q/60000),
                full_moments_available=False)
        out[str(bg)]=dict(moments=raw,increment=increment(raw),unclosed_increment=increment(unclosed),
            truncated_windows=truncations,wide_partial=wide)
    return out


def onset_audit():
    cfg=toml(PROTOCOL)["readings"]
    result={}
    for bg in cfg["background_pixels"]:
        rows=[]
        curves={k:pulse(k,bg) for k in NARROW}
        for level in cfg["threshold_mM"]:
            crossing={k:last_rising_crossing(*v,level) for k,v in curves.items()}
            if crossing["blank"] is None:continue
            rows.append(dict(level_mM=level,crossing_s=crossing,delay_s={k:None if crossing[k] is None else
                crossing[k]-crossing["blank"] for k in ("R1","R2")}))
        result[str(bg)]=rows
    return result


def exponential_fall_bound(t,c,reading_mM,min_later=.2,min_span=15.):
    """Necessary decay bound for y=e^{-t/tau}/tau * u, with u(t)>=0.

    Use the earlier concentration minus and later concentration plus the
    reading allowance: a conservative test of a *specified* exponential.
    The bound is not a fitted time constant for the actual detector.
    """
    t=np.asarray(t,dtype=float);c=np.asarray(c,dtype=float)
    candidates=[]
    for i in range(int(np.argmax(c)),len(c)-1):
        for j in range(i+1,len(c)):
            earlier=c[i]-reading_mM;later=c[j]+reading_mM
            if t[j]-t[i]<min_span or c[j]<min_later or earlier<=later or later<=0:
                continue
            bound=(t[j]-t[i])/math.log(earlier/later)
            candidates.append((bound,i,j,earlier,later))
    if not candidates:return None
    bound,i,j,earlier,later=min(candidates)
    return dict(tau_upper_s=float(bound),T90_upper_s=float(bound*math.log(10)),
        earlier_time_s=float(t[i]),later_time_s=float(t[j]),earlier_reading_mM=float(c[i]),
        later_reading_mM=float(c[j]),earlier_lower_mM=float(earlier),later_upper_mM=float(later))


def detector_audit():
    r=toml(SOURCE)["detector"];cfg=toml(PROTOCOL)["readings"]
    f=r["figure4"]; readings=[]
    for concentration,strip,y in f["markers"]:
        if strip=="lower":
            current=(f["lower_axis_y_0_uA_px"]-y)*.3/(f["lower_axis_y_0_uA_px"]-f["lower_axis_y_0p3_uA_px"])
        else:
            current=.4+(f["upper_axis_y_0p4_uA_px"]-y)*.4/(f["upper_axis_y_0p4_uA_px"]-f["upper_axis_y_0p8_uA_px"])
        readings.append([concentration,current])
    a=np.array(readings); small=a[:,0]<=r["linear_max_mM_text"]
    slope,intercept=np.polyfit(a[small,0],a[small,1],1)
    tau=r["t90_step_s"]/math.log(10)
    bounds={}
    for name in NARROW:
        k,c,u=marker_values(name)
        witness=exponential_fall_bound(k*5,c,cfg["filter_reading_pixels"]*u,
                         cfg["filter_min_later_mM"],cfg["filter_min_interval_s"])
        witness["minimum_later_if_reported_first_order_mM"]=witness["earlier_lower_mM"]*math.exp(
             -(witness["later_time_s"]-witness["earlier_time_s"])/tau)
        witness["compatible_with_reported_T90_as_single_exponential"]=bool(witness["tau_upper_s"]>=tau)
        bounds[name]=witness
    return dict(figure4_readings_mM_uA=readings,illustrative_linear_fit_slope_uA_per_mM=float(slope),
        intercept_uA=float(intercept),fit_max_mM=r["linear_max_mM_text"],
        at_4_mM_measured_uA=float(a[-1,1]),at_4_mM_linear_extrapolated_uA=float(4*slope+intercept),
        pulse_peaks_mM={n:float(marker_values(n)[1].max()) for n in NARROW},
        reported_separate_T90_s=r["t90_step_s"],illustrative_first_order_tau_s=tau,
        first_order_decay_necessary_bound=bounds,
        note="No correction to the pulse concentration or secondary detector convolution. Different calibration, step transients or response shape may explain the incompatibility.")


def kernel_moments(k,dt):
    k=np.asarray(k,dtype=float); age=np.arange(len(k))*dt
    mass=float(k.sum()); mean=float((k@age)/mass)
    return dict(recovery=mass,mean_s=mean,variance_s2=float(k@((age-mean)**2)/mass))


def prior_kernel(derive=False):
    """Archive actual saved A marker ages once; permit a compact source-only replay."""
    old=read_json(PRIOR); p=old["A"]["selected"]["tracer"]
    dt=toml(PROTOCOL)["transfer"]["short_age_dt_s"]
    REF.mkdir(exist_ok=True)
    path=REF/"A_impulse.csv"; manifest=REF/"A_impulse_manifest.json"
    if derive or not (path.exists() and manifest.exists()):
        folder=ROOT/old["run_directory"]/f"A_H{p['depth_cells']}"
        tag=f"dt{p['dt_s']:.8g}_D{p['factor']:.4g}_n{p['particle_count']}_seed{p['seed']}"
        marker=folder/(tag+".f64")
        if sha(marker)!=p["marker_sha256"] or sha(folder/"handoff.json")!=p["handoff_sha256"]:
            raise ValueError("Original A carrier or marker ledger changed")
        ages=np.fromfile(marker,dtype="<f8")
        if len(ages)!=p["particle_count"] or not np.isfinite(ages).all():
            raise ValueError("The archived A run must have complete recovery")
        k=loop.tracer_kernel(ages,dt)
        with path.open("w",newline="") as f:
            w=csv.writer(f,lineterminator="\n"); w.writerow(("age_s","exit_mass_fraction"))
            w.writerows((f"{i*dt:.10g}",f"{weight:.17g}") for i,weight in enumerate(k))
        write_json(manifest,dict(prior_reference_sha256=sha(PRIOR),marker_sha256=p["marker_sha256"],
             handoff_sha256=p["handoff_sha256"],source_run=str(marker.relative_to(ROOT)),dt_s=dt,
             kernel_sha256=sha(path),moments=kernel_moments(k,dt),particle_count=len(ages),
             note="Binned actual recovered A_H16 stock TRT/Julia marker ages; no synthetic replacement."))
    m=read_json(manifest)
    if m["prior_reference_sha256"]!=sha(PRIOR) or m["marker_sha256"]!=p["marker_sha256"] or m["kernel_sha256"]!=sha(path):
        raise ValueError("The compact impulse does not match frozen A")
    age,k=np.loadtxt(path,delimiter=",",skiprows=1).T
    np.testing.assert_allclose(age,np.arange(len(age))*dt,atol=1e-9,rtol=0)
    np.testing.assert_allclose(k.sum(),1,atol=1e-12,rtol=0)
    np.testing.assert_allclose(kernel_moments(k,dt)["mean_s"],p["conditional_mean_s"],rtol=1e-10)
    np.testing.assert_allclose(kernel_moments(k,dt)["variance_s2"],p["conditional_variance_s2"],rtol=1e-3)
    return k,m,old


def scale_kernel(k,ratio,dt):
    """Rescale a weighted delay distribution, retaining mass and mean."""
    k=np.asarray(k,dtype=float)
    if ratio<=0 or not math.isfinite(ratio) or np.any(k<0):
        raise ValueError("Invalid age scaling")
    z=np.arange(len(k))*ratio
    lo=np.floor(z).astype(int); fraction=z-lo
    m=max(lo)+2
    out=np.bincount(lo,weights=k*(1-fraction),minlength=m)
    out+=np.bincount(lo+1,weights=k*fraction,minlength=m)
    np.testing.assert_allclose(out.sum(),k.sum(),rtol=1e-13)
    np.testing.assert_allclose(kernel_moments(out,dt)["mean_s"],ratio*kernel_moments(k,dt)["mean_s"],rtol=1e-13)
    return out


def geometry(name):
    source=toml(SOURCE)["tracer"]
    s=dict(zip(source["columns"],next(row for row in source["rows"] if row[0]==name)))
    s["box_volume_ul"]=s["length_mm"]*s["width_mm"]*s["depth_mm"]
    s["nominal_time_s"]=s["reported_volume_ul"]*60/source["flow_ul_min"]
    s["envelope_volume_ul"]=s["box_volume_ul"]-s["reported_volume_ul"]
    return s


def predict_before_score(k,old):
    """Only the thesis blank, reported volumes and frozen article-A model are used."""
    cfg=toml(PROTOCOL); dt=cfg["readings"]["quadrature_dt_s"]
    a_mean=kernel_moments(k,dt)["mean_s"]
    params=old["A"]["model_parameters"]
    a_envelope=old["independent_storage_bounds"]["A"]["bead_envelope_volume_ul"]
    predictions={};records={}
    for name in ("R1","R2"):
        g=geometry(name); nominal=g["nominal_time_s"]
        one=scale_kernel(k,nominal/a_mean,dt)
        n=1 if name=="R1" else cfg["transfer"]["R2_independent_sections"]
        single=scale_kernel(k,nominal/(n*a_mean),dt)
        base=np.ones(1)
        for _ in range(n):base=np.convolve(base,single)
        models={"volume_plug":loop.tracer_kernel(np.array([nominal]),dt),"prior_trt_independent":base}
        if name=="R2":models["prior_trt_full_memory"]=one
        models["extra_mixed_volume"]=np.convolve(base,loop.exponential_connection(params["extra_mixed_volume_mean_s"],dt))
        if params["rare"]["available"]:
            rare=params["rare"]
            models["common_rare_delay"]=np.convolve(base,loop.rare_delay(rare["rate_A"],rare["mean_event_s"],dt))
            models["bead_volume_rare_delay"]=np.convolve(base,loop.rare_delay(
                rare["rate_A"]*g["envelope_volume_ul"]/a_envelope,rare["mean_event_s"],dt))
        records[name]={"geometry":g,"kernel_moments":{m:kernel_moments(v,dt) for m,v in models.items()}}
        for bg in cfg["readings"]["background_pixels"]:
            _,blank=pulse("blank",bg,dt=dt)
            for label,model in models.items():
                _,c=loop.convolve(blank,model,dt)
                predictions[f"{name}__{label}__bg{bg}"]=c
    return predictions,records


def score_prediction(name,c,bg):
    cfg=toml(PROTOCOL); dt=cfg["readings"]["quadrature_dt_s"]
    to,co=pulse(name,bg,dt=dt);tb,cb=pulse("blank",bg,dt=dt)
    tp=np.arange(len(c))*dt
    o,b,p=(loop.moments(t,y) for t,y in ((to,co),(tb,cb),(tp,c)))
    inc={v:{key:m[key]-b[key] for key in ("mean_s","variance_s2")} for v,m in (("observed",o),("predicted",p))}
    t=np.arange(max(len(to),len(tp)))*dt
    f=np.interp(t,to,loop.cumulative(to,co),left=0,right=1)
    g=np.interp(t,tp,loop.cumulative(tp,c),left=0,right=1)
    score=dict(observed=o,predicted=p,blank=b,incremental_observed=inc["observed"],incremental_predicted=inc["predicted"],
        incremental_mean_relative=inc["predicted"]["mean_s"]/inc["observed"]["mean_s"]-1,
        incremental_variance_relative=inc["predicted"]["variance_s2"]/inc["observed"]["variance_s2"]-1,
        whole_mean_relative=p["mean_s"]/o["mean_s"]-1,whole_variance_relative=p["variance_s2"]/o["variance_s2"]-1,
        whole_CDF_max=float(np.max(abs(f-g))),dose_ratio=p["observed_glucose_umol"]/o["observed_glucose_umol"])
    score["gates"]={key:bool(abs(score[key])<=limit) for key,limit in cfg["gates"].items()}
    score["all_gates"]=all(score["gates"].values())
    return score


def volume_audit(measurements,model_records,old):
    """Conservation-based hold-up requirement, not a fitted volume or RTD."""
    q=toml(SOURCE)["tracer"]["flow_ul_min"]/60
    gate=toml(PROTOCOL)["gates"]["incremental_mean_relative"]
    storage=old["A"]["independent_storage_parameters"]["porosity_scenarios"]
    tau_d=old["independent_storage_bounds"]["A"]["diffusion_time_s"]
    result={}
    for bg,m in measurements.items():
        by={}
        for name in ("R1","R2"):
            g=geometry(name); inc=m["increment"][name]; nominal=g["nominal_time_s"]; volume=inc["mean_s"]*q
            # This is an illustrative comparison of the independently published
            # related-bead pore fractions; they are not Panjan adsorption data.
            baseline=model_records[name]["kernel_moments"]["prior_trt_independent"]
            cases={}
            for key,phi in storage.items():
                add=phi*g["envelope_volume_ul"]/q
                cases[key]=dict(added_time_s=add,mean_s=nominal+add,
                           variance_s2=baseline["variance_s2"]+2*add*tau_d,
                           mean_relative=(nominal+add)/inc["mean_s"]-1)
            by[name]=dict(nominal_mean_s=nominal,measured_increment=inc,
                apparent_liquid_volume_ul=volume,reported_liquid_volume_ul=g["reported_volume_ul"],
                apparent_extra_hold_up_ul=volume-g["reported_volume_ul"],
                implied_additional_accessible_bead_fraction=(volume-g["reported_volume_ul"])/g["envelope_volume_ul"],
                nominal_mean_relative=nominal/inc["mean_s"]-1,
                additional_delay_interval_s=[max(0.,(1-gate)*inc["mean_s"]-nominal),(1+gate)*inc["mean_s"]-nominal],
                independent_storage_cases=cases)
        a=by["R1"]["additional_delay_interval_s"];b=by["R2"]["additional_delay_interval_s"]
        ratio=geometry("R2")["envelope_volume_ul"]/geometry("R1")["envelope_volume_ul"]
        def intersection(f):
            lo=max(a[0],b[0]/f);hi=min(a[1],b[1]/f)
            return dict(feasible=lo<=hi,lower_equivalent_R1_s=lo,upper_equivalent_R1_s=hi)
        c=old["C_scores"][f"measured_blank_bg{bg}"]
        result[bg]=dict(narrow=by,common_positive_delay=intersection(1.),bead_volume_scaled_positive_delay=intersection(ratio),
             envelope_ratio_R2_R1=ratio,prior_article_C=dict(measured_increment=c["incremental_observed"],
                 base_increment=c["incremental_predicted"],apparent_extra_hold_up_ul=(c["incremental_observed"]["mean_s"]-
                    c["incremental_predicted"]["mean_s"])*q))
    article=toml(ARTICLE)["table2"]
    return dict(scenarios=result,article_B_geometry=dict(zip(article["columns"],next(row for row in article["rows"] if row[0]=="B"))),
        note="Reported liquid inventory is used as that inventory, not fitted to the pulse. R1/A likely re-use measurements, so do not treat them as independent calibrations.")


def make_plots(predictions,measurements,detector,old):
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt
    plt.rcParams.update({"font.size":10,"axes.spines.top":False,"axes.spines.right":False})
    FIGURES.mkdir(exist_ok=True)
    dt=toml(PROTOCOL)["readings"]["quadrature_dt_s"]
    fig,axes=plt.subplots(2,2,figsize=(12,8),layout="constrained")
    ax=axes[0,0]
    for name,other,color in (("blank","blank","#666666"),("R1","A","#1773a6")):
        t,c=pulse(name); ax.plot(t/60,c,color=color,lw=2,label=f"Thesis {name}")
        ta,ca=loop.measured(other); ax.plot(ta/60,ca,color=color,ls="--",lw=1.5,label=f"Article {other}")
    ax.set(title="Blank and R1/A: overlapping published traces",xlabel="Time (min)",ylabel="Glucose (mM)",xlim=(.6,5.5))
    ax.legend(frameon=False,fontsize=8);ax.grid(alpha=.18)
    ax=axes[0,1]
    for name,col,lab in (("blank",".62","Thesis blank"),("R2","#1773a6","Thesis R2 (52 mm)")):
        t,c=pulse(name); ax.plot(t/60,c,color=col,lw=1.8,label=lab)
    t,c=loop.measured("C");ax.plot(t/60,c,ls="--",color="#d97706",lw=1.8,label="Article C (81 mm)")
    ax.set(title="Different lengths, nearly the same mean arrival",xlabel="Time (min)",ylabel="Glucose (mM)",xlim=(.7,6))
    ax.legend(frameon=False,fontsize=8);ax.grid(alpha=.18)
    ax=axes[1,0]
    names=("R1 (25 mm)","R2 (52 mm)","C (81 mm)"); x=np.arange(3)
    inc=measurements["0"]["increment"]
    obs=[inc["R1"]["mean_s"],inc["R2"]["mean_s"],old["C_scores"]["measured_blank_bg0"]["incremental_observed"]["mean_s"]]
    ax.bar(x-.18,[geometry("R1")["nominal_time_s"],geometry("R2")["nominal_time_s"],56.5*.6],.34,color="#666666",label="Published liquid V/Q")
    ax.bar(x+.18,obs,.34,color="#1773a6",label="Integrated mean minus blank")
    ax.scatter(x[:2],np.array([.21,.35])*60,marker="D",color="#d97706",zorder=5,label="Thesis reported onset delay")
    ax.set_xticks(x,names);ax.set(ylabel="Time (s)",title="An onset delay is not an integral mean",ylim=(0,42))
    ax.grid(axis="y",alpha=.18);ax.legend(frameon=False,fontsize=8,loc="upper left")
    ax=axes[1,1]
    t,c=pulse("R2");ax.plot(t/60,c,color="black",lw=1.7,label="Measured thesis R2")
    for model,lab,col in (("prior_trt_independent","Two independent A age sections","#1773a6"),
        ("prior_trt_full_memory","Scaled A age memory","#679b92"),("bead_volume_rare_delay","Frozen A delay, bead-volume scale","#d97706")):
        v=predictions[f"R2__{model}__bg0"];ax.plot(np.arange(len(v))*dt/60,v,color=col,lw=1.7,label=lab)
    ax.set(xlabel="Time (min)",ylabel="Glucose (mM)",xlim=(.6,6),title="R2 reduced transfer: mean liquid volume fixed")
    ax.grid(alpha=.18);ax.legend(frameon=False,fontsize=8)
    fig.suptitle("Panjan (2018) and Strniša et al. (2018): dimensional aqueous tracer audit",fontsize=14)
    fig.savefig(FIGURES/"measurement_audit.png",dpi=180);plt.close(fig)

    fig,axes=plt.subplots(1,2,figsize=(10.5,4.5),layout="constrained")
    ax=axes[0]; a=np.array(detector["figure4_readings_mM_uA"])
    ax.plot(a[:,0],a[:,1],"o-",color="#1773a6",label="Fig. 4 marker readings")
    x=np.linspace(0,4.6,100);ax.plot(x,detector["intercept_uA"]+x*detector["illustrative_linear_fit_slope_uA_per_mM"],
           ls="--",color=".5",label="Low-concentration fit, extrapolated")
    ax.axvspan(4.,4.7,color="#d97706",alpha=.12,label="Fig. 25 narrow peak region")
    ax.set(xlabel="Glucose (mM)",ylabel="Current (µA)",title="Separate static detector calibration",xlim=(0,4.8),ylim=(0,1.35))
    ax.legend(frameon=False,fontsize=8);ax.grid(alpha=.18)
    ax=axes[1];w=detector["first_order_decay_necessary_bound"]["blank"]; tau=detector["illustrative_first_order_tau_s"]
    t,c=pulse("blank");ax.plot(t,c,color="#1773a6",lw=2,label="Measured blank")
    # Later signal may not be below this lower envelope for a positive input
    # passed through a single first-order detector of the stated step time.
    ts=np.arange(w["earlier_time_s"],190,.1)
    lower=w["earlier_lower_mM"]*np.exp(-(ts-w["earlier_time_s"])/tau)
    ax.plot(ts,lower,color="#d97706",ls="--",lw=2,label="Minimum if one 58 s T90 filter")
    ax.scatter([w["later_time_s"]],[w["later_upper_mM"]],color="#222",label="Later reading + 5 px",zorder=5)
    ax.vlines(w["later_time_s"],w["later_upper_mM"],w["minimum_later_if_reported_first_order_mM"],color="#b91c1c",lw=2)
    ax.set(xlabel="Time (s)",ylabel="Glucose (mM)",title="Blank falls faster than that single filter",xlim=(125,195),ylim=(0,3.35))
    ax.grid(alpha=.18);ax.legend(frameon=False,fontsize=8)
    fig.suptitle("Detector evidence: calibration and dynamic model are not transferable by assumption",fontsize=12)
    fig.savefig(FIGURES/"detector_audit.png",dpi=180);plt.close(fig)


def run(derive_kernel=False):
    cfg=toml(PROTOCOL)
    k,km,old=prior_kernel(derive_kernel)
    run_dir=ROOT/"results/biocatalysis/panjan_tracer_audit/runs"/(datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")+"-"+uuid4().hex[:8])
    run_dir.mkdir(parents=True)
    inputs=(SOURCE,FIGURE,ARTICLE,loop.MARKERS,PRIOR,PROTOCOL)
    hashes={str(f.relative_to(ROOT)):sha(f) for f in inputs}
    prediction_records=run_dir/"prediction_manifest.json"
    arrays,models=predict_before_score(k,old)
    full=run_dir/"predictions_before_comparison.npz"
    np.savez_compressed(full,**arrays)
    write_json(prediction_records,dict(source_sha256=hashes,prior_kernel=km,
        prediction_sha256=sha(full),reaction_enabled=False,models=models,
        note="Prediction uses only the measured thesis blank and independently reported volumes, plus the already frozen article-A model. No thesis reactor concentrations enter predictions."))
    # All scores below use the full saved curves, not selected matching portions.
    with np.load(full) as saved:
        scores={}
        for key,c in saved.items():
            name,model,background=key.split("__")
            scores[key]=score_prediction(name,c,int(background[2:]))
    m=moment_scenarios();det=detector_audit();dups=duplicate_audit();onset=onset_audit()
    vol=volume_audit(m,models,old)
    make_plots(arrays,m,det,old)
    sample=REF/"predicted_narrow_samples.csv"
    with sample.open("w",newline="") as f:
        w=csv.writer(f,lineterminator="\n");w.writerow(("reactor","model","time_s","predicted_glucose_mM"))
        dt=cfg["readings"]["quadrature_dt_s"]
        for key,c in arrays.items():
            name,model,bg=key.split("__")
            if bg!="bg0":continue
            for time in np.arange(0,421,5):
                i=round(time/dt);w.writerow((name,model,int(time),f"{c[i] if i<len(c) else 0.:.10g}"))
    output=(sample,REF/"A_impulse.csv",REF/"A_impulse_manifest.json",FIGURES/"measurement_audit.png",FIGURES/"detector_audit.png")
    code=tuple(f for f in HERE.glob("*.py"))+(Path(loop.__file__),)
    if hashes!={str(f.relative_to(ROOT)):sha(f) for f in inputs}:
        raise ValueError("Source/protocol changed during the audit")
    final=dict(reaction_enabled=False,new_native_run=False,run_directory=str(run_dir.relative_to(ROOT)),
        versions=dict(platform=platform.platform(),python=sys.version.split()[0],numpy=np.__version__),
        source_sha256=hashes,code_sha256={str(f.relative_to(ROOT)):sha(f) for f in code},
        output_sha256={str(f.relative_to(HERE)):sha(f) for f in output},
        saved_prediction_manifest_sha256=sha(prediction_records),saved_prediction_sha256=sha(full),
        archived_A_impulse=km,duplicate_audit=dups,measurements=m,onset_audit=onset,detector_audit=det,
        geometry_volume_audit=vol,prediction_models=models,prediction_scores=scores,
        frozen_delay_parameters=old["A"]["model_parameters"],
        prior_article_C=old["C_scores"]["measured_blank_bg0"],
        any_R2_reduced_case_pass=any(v["all_gates"] for key,v in scores.items() if key.startswith("R2__")),
        note="R1 and article A are likely repeated published traces based on raw matched markers. R2 is additional dimensional evidence, not a full independent native geometry qualification. All reaction remains disabled.")
    write_json(run_dir/"summary.json",final);write_json(REF/"summary.json",final)
    print("Reference:",REF/"summary.json",flush=True)
    print("R2 reduced-case pass:",final["any_R2_reduced_case_pass"],"reaction enabled: False",flush=True)
    return final


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--derive-kernel",action="store_true",help="Verify and rearchive actual local A marker ages")
    run(parser.parse_args().derive_kernel)
