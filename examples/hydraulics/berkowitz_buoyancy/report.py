"""Compare simulated finite pulses with the independently reported dispersivity."""
import argparse
import math
import platform
import re
from pathlib import Path

import numpy as np

from examples.hydraulics.berkowitz_homogeneous import run as carrier
from pore2process_demo.io import read_json, sha, write_json
from .run import OUT, PROTOCOL
from pore2process_demo.paths import ROOT


def moments(entry, arrival, length_m):
    residence = arrival-entry
    variance_in = entry.var(ddof=1)
    variance_out = arrival.var(ddof=1)
    effective_variance = variance_out-variance_in
    covariance = np.cov(entry,residence,ddof=1)[0,1]
    if not math.isclose(effective_variance, residence.var(ddof=1)+2*covariance, rel_tol=1e-12):
        raise ValueError("Finite-pulse moment accounting failed")
    return dict(mean_residence_min=float(residence.mean()/60),
                outlet_sd_min=float(variance_out**.5/60),
                equivalent_bed_sd_min=float(effective_variance**.5/60),
                covariance_entry_residence_s2=float(covariance),
                alpha_mm=float(500*length_m*effective_variance/residence.mean()**2),
                arrival_p10_p50_p90_min=(np.quantile(arrival,[.1,.5,.9])/60).tolist())


def record(folder):
    c = read_json(folder/"inputs.json")
    started=(folder/"runtime.toml").stat().st_mtime_ns
    if any((folder/name).stat().st_mtime_ns<started for name in ("history.csv","arrivals.f64","entries.f64","outlet_z.f64","run.log")):
        raise ValueError(f"Output predates its input snapshot in {folder}")
    entry = np.fromfile(folder/"entries.f64",dtype="<f8")
    arrivals = np.fromfile(folder/"arrivals.f64",dtype="<f8").reshape(-1,3)
    z = np.fromfile(folder/"outlet_z.f64",dtype="<f8")
    if not (len(entry)==len(arrivals)==len(z)==c["markers"] and
            np.isfinite(arrivals).all() and np.isfinite(z).all() and np.all(arrivals[:,0]>entry)
            and np.all(np.diff(arrivals,axis=1)>0)):
        raise ValueError(f"Incomplete or invalid pulse in {folder}")
    h = np.genfromtxt(folder/"history.csv",delimiter=",",names=True)
    levels = [.05,.1,c["exit_m"]]
    values = [moments(entry,time,length) for time,length in zip(arrivals.T,levels,strict=True)]
    final=arrivals[:,-1]
    groups=[moments(entry[k::16],final[k::16],levels[-1])["alpha_mm"] for k in range(16)]
    time=PROTOCOL["diagnostics"]["height_time_min"]*60
    def interior(column):
        return float(np.interp(time,h["time_s"],h[column]))
    raw_field=np.memmap(Path(c["saved"])/"field.f64",mode="r",dtype="<f8").reshape(-1,7)
    kept=np.fromfile(folder/"active.i32",dtype="<i4")-1
    lower_liquid_fraction=float(np.clip(c["nt"]/2-raw_field[kept,2],0,1).mean())
    liquid_litres=len(kept)*c["repeats"]*c["h"]**3*1000
    mass=c["salt_mass_mg"]*interior("inside")/len(entry)
    lower_mass=mass*interior("fraction_below_midplane")
    interval=carrier.EVIDENCE["pulse"]["collection_interval_min"][0]*60
    edges=np.arange(0, np.ceil(final.max()/interval)+1)*interval
    bins,_=np.histogram(final,bins=edges)
    norm=bins/len(final)*c["feed_s"]/interval
    recovery=float(norm.sum()*interval/c["feed_s"])
    if not math.isclose(recovery,1.0,abs_tol=1e-12):
        raise ValueError("Five-minute samples do not recover the injected mass")
    r=dict(case=folder.name,grid=round(.004/c["h"]),tau=c.get("tau_even",float(read_json(Path(c["saved"])/"raw.json")["tau_even"])),
           gravity_multiplier=c["gravity_multiplier"],gravity_accel_at_feed_m_s2=c["gravity_accel"],
           markers=c["markers"],seed=c["seed"],injection_min=c["feed_s"]/60,
           liquid_Q_mL_min=c["q"]*6e7,feed_NaCl_kg_m3=c["feed_kg_m3"],pulse_mass_mg=c["salt_mass_mg"],
           delta_rho_feed_kg_m3=c["delta_rho_feed_kg_m3"],coupling_s=c["coupling_s"],
           native_sweeps=c["sweeps"],tracer_dt_s=c["dt"],kernel_mm=c["kernel_m"]*1000,bounded_density=bool(c.get("bounded",False)),
           geometry_sha256=c["geometry_sha256"],water_sha256=c["water_sha256"],
           binary_sha256=c["binary_sha256"],native_source_sha256=c["native_source_sha256"],
           stepper_sha256=c["stepper_sha256"],arrivals_sha256=sha(folder/"arrivals.f64"),
           distances_mm=[50,100,200],moments=values,within_realization_alpha_se_mm=float(np.std(groups,ddof=1)/4),
           one_hour=dict(centroid_mm=interior("centroid_z_m")*1000,fraction_below=interior("fraction_below_midplane"),
                         fraction_in_tube=interior("inside")/len(entry),lower_liquid_fraction=lower_liquid_fraction,
                         lower_half_mg_L=lower_mass/(lower_liquid_fraction*liquid_litres),
                         upper_half_mg_L=(mass-lower_mass)/((1-lower_liquid_fraction)*liquid_litres)),
           outlet_fraction_below=float((z<0).mean()),outlet_mean_height_mm=float(z.mean()*1000),
           mass_recovery=recovery,recovered_mg=recovery*c["salt_mass_mg"],
           native_Q_max_relative_error=float(np.max(np.abs(h["native_Q_m3_s"]/c["q"]-1))),
           native_max_Mach=float(h["max_Mach"].max()),native_max_lattice_density_departure=float(h["max_abs_density_departure"].max()),
           max_abs_flow_constraint_correction=float(np.max(abs(h["pump_correction_fraction"]))),
           max_cell_flow_residual=float(h["face_div_relative"].max()),max_plane_Q_error=float(h["plane_Q_relative"].max()),
           max_density_mass_residual_over_dose=float(np.max(abs(h["salt_volume_error_m3"]))/c["source_m3"]),
           max_feed_fraction_used_for_forcing=float(h["feed_fraction_max"].max()),
           saved_face_correction_rms_relative=read_json(Path(c["saved"])/"projection.json")["checks"]["face_rms_correction_relative"],
           concentration_bins=dict(width_min=interval/60,C_over_feed=norm.tolist()))
    early,late=np.quantile(final,[.1,.9])
    r["outlet_height_by_arrival_decile_mm"]={
        "earliest":float(1000*z[final<=early].mean()),
        "latest":float(1000*z[final>=late].mean()),
    }
    run_log=(folder/"run.log").read_text()
    excursion=re.search(r"Furthest upstream voxel (-?\d+); guard (\d+)",run_log)
    if excursion is None:
        raise ValueError("Completion and upstream-guard verification missing")
    r["furthest_upstream_mm"]=int(excursion[1])*c["h"]*1000
    r["downstream_guard_mm"]=int(excursion[2])*c["h"]*1000
    warm=re.search(r"Warm water velocity change relative to saved field ([0-9.eE+-]+)",run_log)
    r["native_warm_velocity_rms_change"]=float(warm[1])
    if "reverse_entry_Q_fraction" in h.dtype.names:
        feed=(h["time_s"]<=c["feed_s"])
        r["max_reverse_entry_during_pulse"]=float(h["reverse_entry_Q_fraction"][feed].max())
    for column in ("redistributed_fraction_of_dose","reverse_entry_Q_fraction"):
        if column in h.dtype.names:
            r["max_"+column]=float(h[column].max())
    if c["gravity_multiplier"]==0:
        old=read_json(Path(c["saved"])/"tracer_m4096_dt1_D1.json")
        archived=old["records"][-1]["alpha_mm"]
        r["archived_passive_alpha_mm"]=archived
        r["passive_relative_difference"]=values[-1]["alpha_mm"]/archived-1
        r["passive_consistency_pass"]=abs(r["passive_relative_difference"])<=PROTOCOL["checks"]["passive_agreement_relative"]
    if (r["max_cell_flow_residual"]>PROTOCOL["checks"]["max_face_divergence_relative"] or
            r["max_plane_Q_error"]>PROTOCOL["checks"]["max_plane_flow_relative"]):
        raise ValueError("Local or plane throughput check failed")
    if r["bounded_density"] and r["max_feed_fraction_used_for_forcing"]>1+1e-12:
        raise ValueError("Reconstructed salt is more concentrated than the inlet")
    checkpoints=np.arange(5,151,5)*60.
    sample=checkpoints[(checkpoints<=h["time_s"][-1]) &
                       (np.interp(checkpoints,h["time_s"],h["inside"])/len(entry)>=.50)]
    r["interior_history"]=dict(time_min=(sample/60).tolist(),
        centroid_mm=(1000*np.interp(sample,h["time_s"],h["centroid_z_m"])).tolist(),
        fraction_below=np.interp(sample,h["time_s"],h["fraction_below_midplane"]).tolist())
    return r


def comparison(records, active, control, target):
    hot,cold=records[active],records[control]
    for field in ("geometry_sha256","water_sha256","markers","seed","liquid_Q_mL_min","pulse_mass_mg","coupling_s","native_sweeps","tracer_dt_s"):
        if hot[field]!=cold[field]:
            raise ValueError(f"Mismatched paired input {field}")
    on=hot["moments"][-1]["alpha_mm"]; off=cold["moments"][-1]["alpha_mm"]
    return dict(active=active,control=control,alpha_off_mm=off,alpha_on_mm=on,
                change_mm=on-off,gap_closed_fraction=1-abs(target-on)/abs(target-off),
                target_relative_error=(on-target)/target,
                passes_published_comparison=abs(on/target-1)<=carrier.PROTOCOL["checks"]["literature_relative"],
                centroid_shift_at_one_hour_mm=hot["one_hour"]["centroid_mm"]-cold["one_hour"]["centroid_mm"],
                fraction_below_change_at_one_hour=hot["one_hour"]["fraction_below"]-cold["one_hour"]["fraction_below"])


def relaxation(folder):
    rows=np.genfromtxt(folder/"relaxation.csv",names=True,delimiter=",")
    last=rows[-1]
    if last["sweeps"]<1536 or (last["weighted_rms_previous_change_m_s"]>=1e-7 and last["sweeps"]<6144):
        raise ValueError(f"Incomplete frozen-plume relaxation in {folder}")
    label=folder.name.removeprefix("audit_")
    result={field:float(last[field]) for field in rows.dtype.names}
    result["density_field_sha256"]=sha(folder.parent/("salt_"+label+".bin"))
    result["online_water_sha256"]=sha(folder.parent/("water_"+label+".bin"))
    return result


def figure(summary, pair, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ref=summary["paper_alpha_mm"]
    cold,hot=(summary["cases"][pair[k]] for k in ("control","active"))
    fig,axes=plt.subplots(1,3,figsize=(11.6,3.7),layout="constrained",gridspec_kw={"width_ratios":[1,1.18,.72]})
    for case,color,label in ((cold,"#586370","Buoyancy off"),(hot,"#007d77","Buoyancy on")):
        series=case["interior_history"]
        axes[0].plot(series["time_min"],series["centroid_mm"],lw=2.3,color=color,label=label)
        curve=case["concentration_bins"]; values=np.array(curve["C_over_feed"])
        axes[1].stairs(values*100,np.arange(len(values)+1)*curve["width_min"],color=color,lw=1.9)
    axes[0].axhline(0,color="#989d9f",lw=.8); axes[0].set(xlabel="Time from pulse start (min)",ylabel="Salt centre height (mm)",ylim=(-13.5,3),title="Height inside the tube")
    axes[0].legend(frameon=False,fontsize=8.4)
    axes[1].set(xlabel="Time from pulse start (min)",ylabel="Outlet salt (% of feed)",title="Five-minute outlet samples")
    a=[pair["alpha_off_mm"],pair["alpha_on_mm"],ref]
    bars=axes[2].bar(["Off","On","Paper"],a,color=["#586370","#007d77","#b57845"],width=.62)
    axes[2].bar_label(bars,fmt="%.1f",padding=3,fontsize=10)
    axes[2].set(title="Outlet spreading",ylabel="Dispersivity (mm)",ylim=(0,max(a)*1.22))
    for ax in axes:
        ax.spines[["right","top"]].set_visible(False)
        ax.grid(axis="y",alpha=.16); ax.set_axisbelow(True)
    fig.savefig(output,dpi=180)
    plt.close(fig)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output",type=Path,default=OUT/"summary.json")
    p.add_argument("--plot",type=Path,default=OUT/"pulse.png")
    args=p.parse_args()
    names=("N6_off","N6_down","N6_up","N8_off","N8_down","N8_low_off","N8_low_down",
           "N8_bounded","N8_kernel","N8_half_off","N8_half","N12_off","N12_bounded")
    target=carrier.EVIDENCE["single_packing"]["coarse_dispersivity_mm"]
    cases={n:record(OUT/n) for n in names}
    if len({v["geometry_sha256"] for v in cases.values()})!=1:
        raise ValueError("Packing changed between runs")
    pairs=[comparison(cases,n,off,target) for n,off in
           (("N6_down","N6_off"),("N6_up","N6_off"),("N8_down","N8_off"),("N8_low_down","N8_low_off"),
            ("N8_bounded","N8_low_off"),("N8_kernel","N8_low_off"),("N8_half","N8_half_off"),("N12_bounded","N12_off"))]
    audits={f"{n}/{stamp}":relaxation(OUT/n/("audit_"+stamp)) for n,stamp in
            (("N8_low_down","early"),("N12_bounded","early"),("N12_bounded","hour"))}
    hydro=read_json(OUT/"uniform_salt/control.json")
    versions=dict(Python=platform.python_version(),NumPy=np.__version__,
                  Julia=read_json(ROOT/".local/julia_runtime.json")["julia_version"],
                  OpenLB=read_json(Path(read_json(OUT/"N12_bounded/inputs.json")["saved"])/"raw.json")["solver"])
    def sensitivity(reference,check):
        previous=cases[reference]["moments"][-1]["alpha_mm"]
        current=cases[check]["moments"][-1]["alpha_mm"]
        return dict(reference=reference,check=check,change_mm=current-previous,
                    relative_change=(current-previous)/previous)
    sensitivities=dict(grid=sensitivity("N8_bounded","N12_bounded"),
                       coupling_interval=sensitivity("N8_bounded","N8_half"),
                       concentration_kernel=sensitivity("N8_bounded","N8_kernel"),
                       feed_bound=sensitivity("N8_low_down","N8_bounded"))
    result=dict(sensitivities=sensitivities,versions=versions,uniform_salt_control=hydro,paper_alpha_mm=target,source="Berkowitz et al. (2009), section 2, paragraph 14",
                cases=cases,pairs=pairs,frozen=audits)
    write_json(args.output,result)
    figure(result,pairs[-1],args.plot)


if __name__=="__main__":
    main()
