"""Source-based volume ledger and Julia tube first-passage moment test.

Equivalent lengths are calculated FROM the measured mean residual only to
test plausibility of tubing dispersion, never used as fitted reactor inputs.
"""
import argparse
import csv
from datetime import datetime,timezone
import math
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import tomllib
from uuid import uuid4

import numpy as np

from pore2process_demo.io import read_json,sha,write_json
from pore2process_demo.paths import ROOT
from examples.adsorption.breakthrough_hypotheses.hydrodynamics.handoff import julia

HERE=Path(__file__).resolve().parent
REF=HERE/"reference"
FIGURES=HERE/"figures"
PROTOCOL=HERE/"protocol.toml"
BAJIC=ROOT/"literature/sources/bajic2017/evidence.toml"
STRNISA=ROOT/"literature/sources/strnisa2018/plumbing_geometry.toml"
PANJAN=ROOT/"literature/sources/panjan2018/plumbing_geometry.toml"
ARTICLE_VOLUMES=ROOT/"literature/sources/strnisa2018/evidence.toml"
STRNISA_PROTOCOL=ROOT/"examples/biocatalysis/strnisa_tracer/protocol.toml"
STRNISA_REF=ROOT/"examples/biocatalysis/strnisa_tracer/reference/summary.json"
PANJAN_REF=ROOT/"examples/biocatalysis/panjan_tracer_audit/reference/summary.json"


def toml(path):
    with path.open("rb") as f:
        return tomllib.load(f)


def volume_ledger():
    b=toml(BAJIC)["inventory"];s=toml(ARTICLE_VOLUMES)["table2"]
    bajic=[];strnisa=[]
    for row in b["rows"]:
        r=dict(zip(b["columns"],row));v=r["length_mm"]*r["width_mm"]*r["depth_mm"]
        r.update(rectangular_product_ul=v,reported_minus_product_ul=r["reactor_volume_ul"]-v)
        bajic.append(r)
    for row in s["rows"]:
        r=dict(zip(s["columns"],row));v=r["length_mm"]*r["width_mm"]*r["depth_mm"]
        # Keep only geometric and liquid volume columns here; not reaction data.
        g={key:r[key] for key in ("reactor","length_mm","width_mm","depth_mm","empty_volume_ul","void_volume_ul","porosity")}
        g.update(rectangular_product_ul=v,reported_minus_product_ul=r["empty_volume_ul"]-v,
                 shape="rectangle" if r["reactor"] in ("A","B","C") else "wide with triangular distributors")
        strnisa.append(g)
    return dict(bajic_rows=bajic,strnisa_rows=strnisa,
        bajic_all_within_table_rounding=all(abs(r["reported_minus_product_ul"])<.05+1e-9 for r in bajic),
        bajic_max_absolute_difference_ul=max(abs(r["reported_minus_product_ul"]) for r in bajic),
        note="All Bajić Table 1 reported reactor volumes equal nominal rectangular L*W*d to rounding, even the distinct wide hexagonal cases. Article Strniša estimates empty volume from dimensions; pre-chambers and port/tube hold-up have no separately measured volume. Numerical agreement is a convention check, not a fixture-only volume measurement.")


def fluid_and_ports():
    article=toml(STRNISA);precedent=toml(BAJIC)["assembly"];thesis=toml(PANJAN)
    d=article["tube_id_mm"]
    if d!=precedent["tube_id_mm"]:
        raise ValueError("The diameter transfer must be reviewed if source records differ")
    f=toml(STRNISA_PROTOCOL)["fluid"]
    q=f["flow_ul_min"]/60
    D=f["glucose_diffusivity_m2_s"]
    area=math.pi*d**2/4;rad=d*1e-3/2
    u=q*1e-9/(area*1e-6)
    ports={}
    for label,src in (("Bajic",precedent),("Strnisa",article),("Panjan",thesis)):
        total=src["upper_plate_thickness_mm"]*src["port_count"]
        ports[label]=dict(nominal_plate_mm=src["upper_plate_thickness_mm"],two_full_depth_length_mm=total,
             equivalent_volume_ul=total*area,mean_s=total*area/q,
             assumed_equal_to_reported_tube_id_mm=d,
             note="Only an illustrative fully wetted bore of this ID. Actual bore and insertion length are unreported.")
    return dict(flow_ul_s=q,diameter_mm=d,radius_m=rad,area_mm2=area,velocity_m_s=u,
        reynolds=f["density_kg_m3"]*u*d*1e-3/f["viscosity_Pa_s"],
        diameter_peclet=u*d*1e-3/D,glucose_diffusivity_m2_s=D,
        radial_diffusion_scale_s=rad**2/D,continuum_variance_slope_s=rad**2/(24*D),
        axial_diffusion_Taylor_slope_term_s=2*D/u**2,
        tube_volume_per_cm_ul=area*10,detector_nominal_volume_ul=thesis["detector_chamber_ul"],
        detector_nominal_hydraulic_time_s=thesis["detector_chamber_ul"]/q,ports=ports,
        note="Article/precedent tube ID is directly reported; for thesis R1/R2 it is a related-diameter scenario. Nominal water/PBS properties from prior protocol, no fitted D.")


def residual_row(case,reading,observed,base,blank,geometry_source,fluid):
    dm=observed["mean_s"]-base["mean_s"]
    dv=observed["variance_s2"]-base["variance_s2"]
    extra=dm*fluid["flow_ul_s"]
    total_length=fluid["ports"][geometry_source]["two_full_depth_length_mm"]
    return dict(case=case,reading=reading,geometry_source=geometry_source,observed=observed,bed_baseline=base,
        blank_complete_mean_s=blank["mean_s"],blank_complete_variance_s2=blank["variance_s2"],
        mean_residual_s=dm,variance_residual_s2=dv,equivalent_extra_liquid_ul=extra,
        equivalent_same_id_length_mm=extra/fluid["area_mm2"],
        bore_id_if_two_full_depth_paths_mm=math.sqrt(4*extra/(math.pi*total_length)) if extra>0 else None,
        single_independent_mixed_volume_variance_s2=dm**2 if dm>0 else None,
        note="Equivalent signed length from residual mean, not an observed line length or a model fit to the pulse.")


def residual_rows():
    s=read_json(STRNISA_REF);p=read_json(PANJAN_REF);cfg=toml(PROTOCOL)["comparison"]
    fluid=fluid_and_ports();out=[]
    for bg in cfg["background_pixels"]:
        label=f"bg{bg}"
        for case,sc in (("A",s["A"]["A_scores"]),("C",s["C_scores"])):
            r=sc[f"measured_blank_{label}"]
            out.append(residual_row(case,label,r["incremental_observed"],r["incremental_predicted"],r["blank"],"Strnisa",fluid))
        m=p["measurements"][str(bg)]
        for name,model,case in (("R1","prior_trt_independent","R1"),("R2","prior_trt_independent","R2"),
                               ("R2","prior_trt_full_memory","R2_memory")):
            out.append(residual_row(case,label,m["increment"][name],p["prediction_models"][name]["kernel_moments"][model],
                                     m["moments"]["blank"],"Panjan",fluid))
    # Deliberate partial-window sensitivity only; not a complete RTD observation.
    cutoff=str(cfg["thesis_short_window_s"]);m=p["measurements"]["0"]["truncated_windows"][cutoff]
    for name,model,case in (("R1","prior_trt_independent","R1"),("R2","prior_trt_independent","R2"),
                           ("R2","prior_trt_full_memory","R2_memory")):
        out.append(residual_row(case,f"window{cutoff}",m["increment"][name],
            p["prediction_models"][name]["kernel_moments"][model],m["moments"]["blank"],"Panjan",fluid))
    return out


def variance_bound(row,radius_m,D):
    """Continuum upper slope; also excludes the sign of a shortening."""
    dm=row["mean_residual_s"];dv=row["variance_residual_s2"]
    slope=radius_m**2/(24*D)
    endpoints=sorted((0.,dm*slope))
    return dict(continuum_slope_s=slope,possible_increment_interval_s2=endpoints,
        necessary_variance_screen_pass=bool(endpoints[0]<=dv<=endpoints[1]),
        diameter_lower_bound_mm=2e3*math.sqrt(24*D*dv/dm) if dm>0 and dv>0 else None,
        note="Bound only for an independent added or contiguous same-ID developed straight tube in the stated high-Pe model. A necessary ID bound is not a sufficient fitting diameter.")


def required_times(rows):
    req={0.}
    for r in rows:
        m=r["mean_residual_s"];b=r["blank_complete_mean_s"]
        req.update((abs(m),b,b+m))
    # Extra reference times for the convergence and visualized finite length law.
    req.update(map(float,np.arange(0,121,2)))
    if min(req)<0:raise ValueError("Unphysical negative baseline plus change")
    return sorted(req)


def julia_pipe_table(run_dir,rows,fluid):
    cfg=toml(PROTOCOL)["tube"];ds=[f*fluid["glucose_diffusivity_m2_s"] for f in cfg["diffusivity_factors"]]
    times=required_times(rows)
    def csvarg(items):return ",".join(format(x,".17g") for x in items)
    command=[*julia(1),str(HERE/"tube_moments.jl"),csvarg(cfg["radial_cells"]),
             format(fluid["radius_m"],".17g"),csvarg(ds),csvarg(times)]
    path=run_dir/"tube_moments.csv"; log=run_dir/"julia.log"
    with path.open("w",newline="") as out,log.open("w") as err:
        subprocess.run(command,cwd=ROOT,stdout=out,stderr=err,check=True,timeout=180)
    table={}
    with path.open(newline="") as f:
        for v in csv.DictReader(f):
            r={k:(int(x) if k=="radial_cells" else float(x)) for k,x in v.items()}
            table[(r["radial_cells"],r["diffusivity_m2_s"],r["nominal_mean_s"])]=r
    if len(table)!=len(cfg["radial_cells"])*len(ds)*len(times):raise ValueError("Missing Julia moments")
    for n in cfg["radial_cells"]:
        for d in ds:
            for t in times:
                if (n,d,t) not in table:raise ValueError("Julia returned different inputs")
    return table,dict(command=command,table_sha256=sha(path),cells=cfg["radial_cells"],
        diffusivity_m2_s=ds,nominal_mean_times_s=times,solver="Julia radial spectral finite-volume analytical-flow reference")


def evaluate_table(table,record,rows,fluid):
    cfg=toml(PROTOCOL)["tube"];ns=record["cells"];ds=record["diffusivity_m2_s"];times=record["nominal_mean_times_s"]
    fine=ns[-1];previous=ns[-2]
    def get(n,d,t):return table[(n,d,t)]
    mean_err=max(abs(v["mean_s"]/v["nominal_mean_s"]-1) for v in table.values() if v["nominal_mean_s"]>0)
    slope_err=max(abs(get(fine,d,0.)["asymptotic_slope_s"]/get(fine,d,0.)["continuum_slope_s"]-1) for d in ds)
    refinement=max(abs(get(fine,d,t)["variance_s2"]/get(previous,d,t)["variance_s2"]-1)
        for d in ds for t in times if t>0.)
    convergence=dict(max_relative_mean_error=mean_err,max_finest_variance_change_relative=refinement,
        max_finest_slope_error_relative=slope_err,finest_cells=fine,
        mean_pass=mean_err<cfg["mean_relative"],refinement_pass=refinement<cfg["finest_variance_change_relative"],
        closed_form_slope_pass=slope_err<cfg["asymptotic_slope_relative"])
    if not all(convergence[k] for k in ("mean_pass","refinement_pass","closed_form_slope_pass")):
        raise ValueError("Julia independent tube reference failed qualification")
    results=[]
    for r in rows:
        m=r["mean_residual_s"];b=r["blank_complete_mean_s"]
        cases={}
        for factor,d in zip(cfg["diffusivity_factors"],ds):
            bound=variance_bound(r,fluid["radius_m"],d)
            vstand=get(fine,d,abs(m))["variance_s2"] if m>=0 else None
            vblank=get(fine,d,b)["variance_s2"]
            vext=get(fine,d,b+m)["variance_s2"]-vblank
            lo,hi=bound["possible_increment_interval_s2"]
            if not lo-1e-7<=vext<=hi+1e-7 or (vstand is not None and vstand>hi+1e-7):
                raise ValueError("Finite tube violates its continuum sign/slope bound")
            cases[str(factor)]=dict(**bound,independent_added_segment_variance_s2=vstand,
                    contiguous_extension_from_entire_blank_equivalent_s2=vext,
                    blank_all_assigned_tube_variance_s2=vblank,
                    finite_extension_upper_bound_s2=hi)
            # Also allow a hypothetical full-depth bore diameter chosen from
            # the mean alone. Do not confuse it with the reported tube ID.
            bore=r["bore_id_if_two_full_depth_paths_mm"]
            needed=bound["diameter_lower_bound_mm"]
            if bore is not None:
                bore_bound=variance_bound(r,bore*1e-3/2,d)
                cases[str(factor)]["full_depth_mean_equivalent_bore_upper_s2"]=bore_bound["possible_increment_interval_s2"][1]
                cases[str(factor)]["full_depth_mean_equivalent_bore_necessary_pass"]=bore_bound["necessary_variance_screen_pass"]
            else:
                cases[str(factor)]["full_depth_mean_equivalent_bore_upper_s2"]=None
                cases[str(factor)]["full_depth_mean_equivalent_bore_necessary_pass"]=False
            cases[str(factor)]["cylindrical_length_at_necessary_diameter_mm"]=r["equivalent_extra_liquid_ul"]/(math.pi*needed**2/4) if needed else None
        results.append(dict(**r,diffusivity_cases=cases))
    return results,convergence


def make_plots(results,fluid,table,n):
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt
    plt.rcParams.update({"font.size":10,"axes.spines.top":False,"axes.spines.right":False})
    FIGURES.mkdir(exist_ok=True)
    main={r["case"]:r for r in results if r["reading"]=="bg0" and r["case"]!="R2_memory"}
    names=("A","R1","R2","C");x=np.arange(len(names)); labels=("A (article)","R1 (thesis)","R2 (thesis)","C (article)")
    fig,(a,b)=plt.subplots(1,2,figsize=(11.5,4.8),layout="constrained")
    vals=[main[v]["equivalent_extra_liquid_ul"] for v in names]
    bars=a.bar(x,vals,color="#1773a6",width=.56)
    for bar,v in zip(bars,vals):
        a.text(bar.get_x()+bar.get_width()/2,v+(.46 if v>=0 else -.46),f"{v:+.2f}",ha="center",va="bottom" if v>=0 else "top",fontsize=9)
    a.axhline(0,color=".3",lw=1)
    a.axhline(fluid["ports"]["Strnisa"]["equivalent_volume_ul"],color="#d97706",ls="--",label="Two full 15 mm × 0.50 mm-ID bores*")
    a.axhline(fluid["ports"]["Panjan"]["equivalent_volume_ul"],color="#679b92",ls=":",lw=2,label="Two full 20 mm × 0.50 mm-ID bores*")
    a.set_xticks(x,labels,rotation=12);a.set(ylabel="Required apparent liquid change (µL)",title="Tracer mean-time residual",ylim=(-4.5,21.5));a.grid(axis="y",alpha=.2)
    a.legend(frameon=False,fontsize=8,loc="upper left")
    need=[main[v]["variance_residual_s2"] for v in names]
    upper=[main[v]["diffusivity_cases"]["1.0"]["finite_extension_upper_bound_s2"] for v in names]
    upper_lowD=[main[v]["diffusivity_cases"]["0.8"]["finite_extension_upper_bound_s2"] for v in names]
    b.bar(x-.22,need,.30,label="Measured variance residual",color="#b91c1c")
    b.bar(x+.08,upper,.26,label="0.50 mm tube maximum",color="#1773a6")
    b.bar(x+.32,upper_lowD,.20,label="Maximum at 0.8 × D",color="#679b92")
    b.set_xticks(x,labels,rotation=12);b.set(ylabel="Variance increment (s²)",title="Same-ID tubing cannot supply required spread",ylim=(0,max(need)*1.2));b.grid(axis="y",alpha=.2)
    b.legend(frameon=False,fontsize=8)
    fig.suptitle("Plumbing hold-up screen (A and R1 likely repeated traces)\n*Bore volumes are illustrative, not measured port dimensions",fontsize=12)
    fig.savefig(FIGURES/"plumbing_screen.png",dpi=180);plt.close(fig)

    fig,ax=plt.subplots(figsize=(8,4.8),layout="constrained")
    D=fluid["glucose_diffusivity_m2_s"];ts=np.arange(0,121,2,dtype=float)
    var=[table[(n,D,float(t))]["variance_s2"] for t in ts]
    ax.plot(ts,var,lw=2,color="#1773a6",label="Finite circular tube (Julia)")
    ax.plot(ts,ts*fluid["continuum_variance_slope_s"],ls="--",color=".5",lw=2,label="Continuum upper slope")
    r=main["R2"];m=r["mean_residual_s"];b=r["blank_complete_mean_s"]
    p=table[(n,D,m)]["variance_s2"]
    ax.scatter(m,p,color="#d97706",zorder=5,label=f"R2 independent added tube ({m:.1f} s)")
    vb=table[(n,D,b)]["variance_s2"];vplus=table[(n,D,b+m)]["variance_s2"]
    ax.plot([b,b+m],[vb,vplus],"o-",color="#679b92",lw=2,label="Extension from whole-blank mean equivalent")
    ax.set(xlabel="Nominal tube transit time V/Q (s)",ylabel="Tube-only exit-time variance (s²)",
           title="0.50 mm-ID Poiseuille tube with independent radial diffusion",xlim=(0,122),ylim=(0,465))
    ax.legend(frameon=False,fontsize=9,loc="upper left");ax.grid(alpha=.2)
    fig.savefig(FIGURES/"tube_moments.png",dpi=180);plt.close(fig)


def run():
    inputs=(BAJIC,STRNISA,PANJAN,ARTICLE_VOLUMES,STRNISA_PROTOCOL,STRNISA_REF,PANJAN_REF,PROTOCOL)
    hashes={str(f.relative_to(ROOT)):sha(f) for f in inputs}
    run_dir=ROOT/"results/biocatalysis/bajic_plumbing_audit/runs"/(datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")+"-"+uuid4().hex[:8])
    run_dir.mkdir(parents=True)
    REF.mkdir(exist_ok=True)
    fluid=fluid_and_ports();rows=residual_rows();ledger=volume_ledger()
    write_json(run_dir/"inputs.json",dict(source_sha256=hashes,fluid=fluid,rows=rows,protocol=toml(PROTOCOL)))
    table,record=julia_pipe_table(run_dir,rows,fluid)
    results,numerical=evaluate_table(table,record,rows,fluid)
    make_plots(results,fluid,table,numerical["finest_cells"])
    shutil.copyfile(run_dir/"tube_moments.csv",REF/"tube_moments.csv")
    s=read_json(STRNISA_REF);p=read_json(PANJAN_REF)
    codes=tuple(f for f in HERE.iterdir() if f.is_file() and f.suffix in (".py",".jl"))
    outputs=(REF/"tube_moments.csv",FIGURES/"plumbing_screen.png",FIGURES/"tube_moments.png")
    if hashes!={str(f.relative_to(ROOT)):sha(f) for f in inputs}:raise ValueError("Input evidence changed during audit")
    out=dict(reaction_enabled=False,new_porous_OpenLB_run=False,actual_new_Julia_tube_run=True,
        run_directory=str(run_dir.relative_to(ROOT)),versions=dict(platform=platform.platform(),
            python=sys.version.split()[0],numpy=np.__version__,julia=subprocess.check_output([*julia(1),"--version"],text=True).strip()),
        source_sha256=hashes,code_sha256={str(f.relative_to(ROOT)):sha(f) for f in codes},
        output_sha256={str(f.relative_to(HERE)):sha(f) for f in outputs},
        input_manifest_sha256=sha(run_dir/"inputs.json"),tube_run=record,numerical=numerical,fluid=fluid,
        source_volume_ledger=ledger,residuals=results,
        prior_common_mean_feasibility_A_C=s["positive_delay_mean_screen"],
        prior_common_mean_feasibility_R1_R2={k:{name:v[name] for name in ("common_positive_delay","bead_volume_scaled_positive_delay")}
            for k,v in p["geometry_volume_audit"]["scenarios"].items()},
        any_same_id_tube_necessary_variance_pass=any(v["necessary_variance_screen_pass"] for r in results for v in r["diffusivity_cases"].values()),
        note="Source-derived diameter; no measured relative tube length or port bore. Conservative equivalent-length, high-Pe straight-tube moment diagnostics only. Actual fixture geometry and detector transfer remain unqualified; no reaction.")
    write_json(run_dir/"summary.json",out);write_json(REF/"summary.json",out)
    print("Reference:",REF/"summary.json",flush=True)
    print("Any same-ID tube necessary-variance pass:",out["any_same_id_tube_necessary_variance_pass"],flush=True)
    return out


if __name__=="__main__":
    argparse.ArgumentParser(description=__doc__).parse_args()
    run()
