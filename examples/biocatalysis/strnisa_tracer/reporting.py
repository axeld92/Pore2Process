"""Publish compact retained evidence and figures from the frozen tracer run."""
import argparse
import csv
from pathlib import Path
import platform
import subprocess
import sys

import numpy as np

from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT
from . import compare, loop, native, source_audit, storage_audit
from .handoff import julia

HERE=native.HERE
REF=HERE/"reference"
FIGURES=HERE/"figures"


def compact(r):
    h=r["hydraulic"]
    f=r["handoff"]
    t=r["tracer"]
    return dict(hydraulic=h,handoff={k:f[k] for k in ("checks","carrier_qualified","pore_volume_m3",
        "ideal_mean_s","removed_corner_only_cells","seam_flow_m3_s","cells","dx_m",
        "field_sha256","face_projection_s")},tracer={k:v for k,v in t.items() if k!="command"})


def mean_delay_bound(A_scores,C_scores,plan,scale):
    """Post-holdout feasibility of any nonnegative additive mean-delay law.

    The independent mean gate imposes an interval even without a choice of
    the delay's distribution. A common kernel and bead-volume scaling
    must have the very same A parameter in the intersection of intervals.
    No fitted replacement or curve modification is made.
    """
    limit=plan["protocol"]["comparison"]["incremental_mean_relative"]
    out={}
    for bg in compare.BACKGROUNDS:
        intervals={}
        for name,scores in (("A",A_scores),("C",C_scores)):
            b=scores[f"measured_blank_bg{bg}"]
            obs=b["incremental_observed"]["mean_s"]
            base=b["incremental_predicted"]["mean_s"]
            intervals[name]=[max(0.,(1-limit)*obs-base),(1+limit)*obs-base]
        def intersection(factor):
            a,c=intervals["A"],intervals["C"]
            lo,hi=max(a[0],c[0]/factor),min(a[1],c[1]/factor)
            return dict(feasible=lo<=hi,required_A_lower_s=lo,allowed_A_upper_s=hi)
        out[str(bg)]=dict(added_mean_intervals_s=intervals,common=intersection(1.),
                          bead_volume_scaled=intersection(scale))
    return out


def prediction_samples(base):
    # These samples are derived from the already saved full forward predictions,
    # not reconstructed from fitted or observed C curves. Report at 1 s cadence.
    out=REF/"prediction_samples.csv"
    with out.open("w",newline="") as f:
        writer=csv.writer(f,lineterminator="\n")
        writer.writerow(("reactor","model","background_pixels","time_s","glucose_mM"))
        for name,file in (("A","A_predictions.npz"),("C","C_prediction_before_comparison.npz")):
            with np.load(base/file) as curves:
                for label,c in curves.items():
                    if not label.endswith("_bg0"):
                        continue
                    model=label.split("_bg")[0]
                    stride=round(1/compare.DT)
                    # Entire plotted experiment and any substantive predicted late tail.
                    n=max(401*stride,int(np.flatnonzero(c>1e-5)[-1])+1)
                    for i in range(0,n,stride):
                        writer.writerow((name,model,0,f"{i*compare.DT:.8g}",f"{c[i] if i<len(c) else 0.:.10g}"))
    return out


def plot(base,summary):
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt
    from matplotlib.ticker import PercentFormatter

    colors={"measured_blank":"#1773a6","extra_mixed_volume":"#679b92",
            "common_rare_delay":"#d97706","solid_volume_rare_delay":"#8355a3"}
    names={"measured_blank":"TRT + measured blank",
           "extra_mixed_volume":"Extra mixed volume (A mean)",
           "common_rare_delay":"Common rare delays (A moments)",
           "solid_volume_rare_delay":"Bead-volume scaling (A moments)"}
    plt.rcParams.update({"font.size":10,"axes.spines.top":False,"axes.spines.right":False})
    fig,ax=plt.subplots(2,2,figsize=(12,8),layout="constrained",gridspec_kw={"height_ratios":[1.2,1]})
    for col,(reactor,file,scores) in enumerate((("A","A_predictions.npz",summary["A"]["A_scores"]),
                              ("C","C_prediction_before_comparison.npz",summary["C_scores"]))):
        a=ax[0,col]; b=ax[1,col]
        t,c=loop.measured("blank")
        a.plot(t/60,c,color=".62",ls="--",lw=1.3,label="Measured blank loop")
        t,c,_=loop.digitized(reactor)
        a.plot(t/60,c,"o",ms=3.7,color="#222222",label=f"Measured reactor {reactor}",zorder=5)
        with np.load(base/file) as curves:
            for model in compare.MODELS:
                key=model+"_bg0"
                if key not in curves or (reactor=="A" and model=="solid_volume_rare_delay"):
                    continue
                p=curves[key]
                a.plot(np.arange(len(p))*compare.DT/60,p,lw=1.5,ls="-" if model=="measured_blank" else ":",
                       color=colors[model],label=names[model])
        a.set(xlim=(.5,6.5),ylim=(-.05,5.0),xlabel="Time after injection (min)",ylabel="Glucose (mM)",
              title="A: short, 23.9 mm" if reactor=="A" else "C: long, 81.0 mm (held out)")
        a.legend(frameon=False,fontsize=8,loc="upper right")
        a.grid(alpha=.18)
        available=[m for m in compare.MODELS if f"{m}_bg0" in scores and (reactor!="A" or m!="solid_volume_rare_delay")]
        short_names={"measured_blank":"TRT + blank","extra_mixed_volume":"Mixed volume",
                     "common_rare_delay":"Rare delays","solid_volume_rare_delay":"Bead-volume law"}
        # Moment ratios make the two distinct units and gates explicit.
        pos=np.arange(len(available))
        for i,metric in enumerate(("incremental_mean_relative","incremental_variance_relative")):
            ratio=[1+scores[m+"_bg0"][metric] for m in available]
            b.bar(pos+(-.18 if i==0 else .18),ratio,width=.34,
                  color=("#2466a8","#cc7a28")[i],label=("Mean","Variance")[i])
        b.axhline(1.,color="#333",lw=1)
        b.axhspan(.90,1.10,color="#999",alpha=.13,label="±10% variance limit")
        b.axhspan(.95,1.05,color="#444",alpha=.12,label="±5% mean limit")
        b.set_xticks(pos,[short_names[m] for m in available],rotation=12,ha="right")
        b.yaxis.set_major_formatter(PercentFormatter(xmax=1))
        b.set(ylabel="Predicted / blank-subtracted observation",title="Increment over shared measured blank",ylim=(0,max(1.35,
            max(1+scores[m+"_bg0"][x] for m in available for x in ("incremental_mean_relative","incremental_variance_relative"))+.12)))
        b.grid(axis="y",alpha=.18)
        b.legend(frameon=False,fontsize=8,loc="upper left")
    fig.suptitle("Strniša et al. (2018): inert glucose pulse, 100 µL/min, 24 °C\n"
                 "All long-reactor predictions frozen from A; calibrated delays are diagnostics",fontsize=14)
    fig.savefig(FIGURES/"loop_comparison.png",dpi=180)
    plt.close(fig)

    # Independent check of the source's normalized panel.
    import tomllib
    a=tomllib.loads(source_audit.FIG3A.read_text())
    markers=np.array(a["markers"],dtype=float)
    x=(markers[:,1]-a["axis_time_0_px"])/(a["axis_time_1_px"]-a["axis_time_0_px"])
    y=(a["axis_conc_0_px"]-markers[:,2])/(a["axis_conc_0_px"]-a["axis_conc_1_px"])
    t,c,_=loop.digitized("A")
    mean=loop.moments(*loop.measured("A"))["mean_s"]
    fig,ax=plt.subplots(figsize=(7,4.3),layout="constrained")
    ax.plot(t/mean,c/c.max(),color="#1773a6",lw=2,label=f"Fig. 3b A / full-loop mean ({mean:.1f} s)")
    ax.scatter(x,y,marker="^",s=18,color="#222",label="Fig. 3a experimental triangles",zorder=5)
    ax.set(xlim=(0,2.55),ylim=(-.03,1.05),xlabel="Time divided by full-loop mean",
           ylabel="Concentration / peak",title="The source's normalized experimental profile retains the measurement loop")
    ax.legend(frameon=False,loc="upper right",fontsize=9); ax.grid(alpha=.18)
    fig.savefig(FIGURES/"dimensionless_audit.png",dpi=180)
    plt.close(fig)


def publish(base):
    result=read_json(base/"summary.json")
    a=compare.verify_freeze(base)
    if a!=result["A"] or sha(base/"C_prediction_before_comparison.npz")!=result["C_manifest"]["predictions_sha256"]:
        raise ValueError("Frozen predictions changed")
    plan=read_json(base/"plan.json")
    REF.mkdir(exist_ok=True); FIGURES.mkdir(exist_ok=True)
    curve_file=prediction_samples(base)
    plot(base,result)
    short=a["selected"]
    long=result["C_manifest"]["carrier"]
    build=read_json(base/f"A_H{a['selected_depth_cells']}"/"build.json")
    # Keep reference small: numerical checks, all comparison scores and hashes,
    # not raw fields, exit-time arrays or native executable.
    b={k:v for k,v in a.items() if k not in ("grid_records","selected","fine_step_record","diffusivity_sensitivity")}
    b.update(grid_records={h:compact(r) for h,r in a["grid_records"].items()},
             selected=compact(short),fine_step_record=compact(a["fine_step_record"]) if a["fine_step_record"] else None,
             diffusivity_sensitivity={k:{"record":compact(v["record"]),"baseline":v["baseline"]}
                                      for k,v in a["diffusivity_sensitivity"].items()})
    c={k:v for k,v in result["C_manifest"].items() if k!="carrier"}
    c["carrier"]=compact(long)
    source_files=(native.SOURCE,native.PROTOCOL,loop.MARKERS,source_audit.FIG3A,storage_audit.EVIDENCE)
    code_files=tuple(f for f in HERE.rglob("*") if f.is_file() and f.suffix in (".py",".jl",".cpp",".mk")) + (
        HERE/"openlb/Makefile", *(ROOT/name for name in (
        "examples/adsorption/transport_root_cause/hydraulics/stock_omp.py",
        "examples/adsorption/breakthrough_hypotheses/hydrodynamics/handoff.py",
        "examples/adsorption/disordered_packing/mesh.py",
        "examples/adsorption/disordered_refinement/balance.jl",
        "examples/adsorption/path_exchange/track.jl",
        "julia/adsorption/PorousFaceTransport.jl",
        )))
    out_files=(curve_file,FIGURES/"loop_comparison.png",FIGURES/"dimensionless_audit.png")
    rejected={}
    for folder in sorted(base.glob("*_attempt_*")):
        if not folder.is_dir():
            continue
        raw=folder/"raw.json"
        req=folder/"request.json"
        rejected[folder.name]=dict(raw=read_json(raw) if raw.is_file() else None,
            raw_sha256=sha(raw) if raw.is_file() else None,
            request=read_json(req) if req.is_file() else None,
            note="Preserved attempt, not used as an accepted carrier or independent tracer result.")
    record=dict(reaction_enabled=False,run_directory=str(base.relative_to(ROOT)),
        plan=plan,versions=dict(platform=platform.platform(),python=sys.version.split()[0],numpy=np.__version__,
            julia=subprocess.check_output([*julia(),"--version"],text=True).strip(),
            compiler=build["compiler"].splitlines()[0],openlb="1.9r0",binary_sha256=build["binary_sha256"],
            stock_vendor_integrity=build["vendor_integrity"]),
        compact_record_sha256={name:sha(base/name) for name in ("plan.json","frozen_A.json","A_predictions.npz",
                                        "C_prediction_before_comparison.npz","C_prediction_manifest.json","summary.json")},
        source_sha256={str(f.relative_to(ROOT)):sha(f) for f in source_files},
        code_sha256={str(f.relative_to(ROOT)):sha(f) for f in code_files},
        output_sha256={str(f.relative_to(HERE)):sha(f) for f in out_files},
        rejected_native_attempts=rejected,
        A=b,C_manifest=c,C_scores=result["C_scores"],C_table1=result["C_table1"],
        independent_storage_bounds=result["independent_storage_bounds"],
        positive_delay_mean_screen=mean_delay_bound(a["A_scores"],result["C_scores"],plan,a["solid_volume_ratio_C"]),
        independent_C_pass=result["independent_C_pass"],reaction_gate_pass=result["reaction_gate_pass"],
        note=result["note"])
    write_json(REF/"summary.json",record)
    print(REF/"summary.json",flush=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("base",type=Path)
    publish(parser.parse_args().base.resolve())
