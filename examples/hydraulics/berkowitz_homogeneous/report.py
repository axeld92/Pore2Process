"""Compare completed single-packing tracer predictions with the published estimates."""
import argparse
from functools import lru_cache
import math
from pathlib import Path
import platform

import numpy as np

from .run import EVIDENCE, HERE, OUTPUT, PROTOCOL, centres_path
from examples.adsorption.breakthrough_hypotheses.hydrodynamics import native, packing
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT


@lru_cache
def sphere_check(path, period, radius):
    """Independent all-pairs minimum-image check in bead diameters."""
    points = np.loadtxt(path)
    gap = np.inf
    for i in range(len(points)-1):
        offset = points[i+1:]-points[i]
        offset[:,0] -= period*np.rint(offset[:,0]/period)
        gap = min(gap, float(np.min(np.einsum("ij,ij->i",offset,offset))))
    separation = math.sqrt(gap)
    clearance = float(radius - np.max(np.hypot(points[:,1],points[:,2])) - .5)
    if separation < 1-1e-7 or clearance < -1e-7:
        raise ValueError("Independent bead or wall overlap check failed")
    return dict(min_centre_distance_in_beads=separation, min_wall_clearance_in_beads=clearance)


def finite_pulse_bins(arrivals, edges, width):
    """Effluent C/Cfeed in collection bins for a rectangular inlet pulse."""
    return np.array([np.maximum(0, np.minimum(arrivals+width,hi)-np.maximum(arrivals,lo)).mean()
                     /(hi-lo) for lo,hi in zip(edges[:-1],edges[1:],strict=True)])


@lru_cache
def check_native_voxels(folder):
    p = read_json(folder/"carrier.json")
    dp = p["bead_diameter_m"]
    material = packing.voxel_material(np.loadtxt(centres_path(p))*dp,tuple(p["shape"]),p["dx_m"],
                                      radius=p["tube_radius_m"],diameter=dp)
    field = np.memmap(folder/"field.f64",mode="r",dtype="<f8").reshape(-1,7)
    expected = np.argwhere(material==1)
    if not np.array_equal(expected,field[:,:3]):
        raise ValueError(f"Independent physical-coordinate pore raster disagrees with {folder}")
    return len(expected)


def reference_cdf(t, mean, length_mm, alpha_mm):
    if t<=0:
        return 0.
    ratio = t/mean
    k = math.sqrt(length_mm/(4*alpha_mm*ratio))
    return .5*(math.erfc((1-ratio)*k)+math.exp(length_mm/alpha_mm)*math.erfc((1+ratio)*k))


def reference_bins(edges, width, mean, length_mm, alpha_mm):
    x,w = np.polynomial.legendre.leggauss(8)
    data = []
    for lo,hi in zip(edges[:-1],edges[1:],strict=True):
        t = (hi+lo)/2 + x*(hi-lo)/2
        data.append(sum(weight*(reference_cdf(z,mean,length_mm,alpha_mm)-
                                reference_cdf(z-width,mean,length_mm,alpha_mm))
                        for z,weight in zip(t,w,strict=True))/2)
    return np.array(data)


def case(beads, period, seed, grid, dt=1., diffusion=1., low_compressibility=False):
    suffix = "_tau060" if low_compressibility else ""
    folder = OUTPUT/f"{beads}_L{period:g}_s{seed}_N{grid}{suffix}"
    root = f"tracer_m{PROTOCOL['transport']['markers']}_dt{dt:g}_D{diffusion:g}"
    p, face, tracer = (read_json(folder/name) for name in ("carrier.json","projection.json",root+".json"))
    raw = read_json(folder/"raw.json")
    path = folder/(root+".f64")
    if not (sha(path)==tracer["arrivals_sha256"] and sha(folder/"field.f64")==p["native_field_sha256"]
            and sha(folder/"carrier.json")==face["carrier_sha256"]
            and sha(folder/"projection.json")==tracer["projection_sha256"]
            and sha(centres_path(p))==p["packing_sha256"]):
        raise ValueError(f"Changed inputs in {folder}")
    levels = tracer["records"]
    t = np.fromfile(path,dtype="<f8").reshape(tracer["markers"],len(levels))
    for column,r in enumerate(levels):
        x = t[:,column]
        if not (np.isfinite(x).all() and len(x)==r["arrivals"]==tracer["markers"]):
            raise ValueError("Incomplete tracer recovery")
        np.testing.assert_allclose([x.mean(),x.var(ddof=1),
            1000*r["length_m"]*x.var(ddof=1)/(2*x.mean()**2)],
            [r["mean_s"],r["variance_s2"],r["alpha_mm"]], rtol=1e-12)
    last = levels[-1]
    target = EVIDENCE["single_packing"][beads+"_dispersivity_mm"]
    age_error = last["mean_s"]/last["nominal_mean_s"]-1
    delta = last["alpha_mm"]/target-1
    audit = sphere_check(centres_path(p),period,p["tube_radius_m"]/p["bead_diameter_m"])
    water = PROTOCOL["water"]
    scale = p["native_length_per_physical_length"]
    velocity_scale = p["native_speed_per_physical_speed"]
    gradient = raw["gradient_Pa_m"][0]*water["density_kg_m3"]/1000*scale/velocity_scale**2
    return dict(case=folder.name,grid=grid,cells_per_bead=p["bead_diameter_m"]/p["dx_m"],
        tau_even=raw["tau_even"],seed=seed,period_beads=period,dt_factor=dt,
        D_factor=diffusion,D_m2_s=tracer["D_m2_s"],dt_s=tracer["dt_s"],diameter_m=p["bead_diameter_m"],
        independently_checked_native_pore_voxels=check_native_voxels(folder),
        geometric_porosity=p["geometric_porosity"],voxel_porosity=p["voxel_porosity"],
        connected_porosity=face["connected_porosity"],sphere_count=p["sphere_count"],**audit,
        flow_mL_min=p["achieved_flow_m3_s"]*60e6,pore_Pe=p["pore_Pe"]/diffusion,pore_Re=p["pore_Re"],
        face_checks=face["checks"],corner_only_fraction=face["corner_only_fraction"],
        excluded_axial_flux_fraction=face["excluded_native_axial_flow_fraction"],
        inlet_reverse_fraction=face["entry_reverse_flux_fraction"],native_max_Mach=raw["max_Mach"],
        pressure_gradient_Pa_m=gradient,max_native_density_departure=raw["max_density_departure"],
        steps=raw["steps"],markers=tracer["markers"],exit=levels,reported_alpha_mm=target,
        mean_time_relative_error=age_error,alpha_relative_error=delta,
        within_literature_screen=abs(delta)<=PROTOCOL["checks"]["literature_relative"],
        packing_sha256=p["packing_sha256"],native_field_sha256=p["native_field_sha256"],
        arrivals_sha256=sha(path)),t[:,-1]


def summarize():
    config = PROTOCOL["packing"]
    base_seed,other_seed = config["baseline_seed"],config["second_seed"]
    pulse = EVIDENCE["pulse"]
    width, interval = pulse["pulse_min"][0]*60,pulse["collection_interval_min"][0]*60
    storage = EVIDENCE["column"]
    length = EVIDENCE["single_packing"]["length_mm"]
    Q = pulse["flow_ml_min"][0]*1e-6/60
    age = math.pi*(storage["inner_diameter_mm"]*.0005)**2*length*.001*storage["porosity"]/Q
    summary = dict(schema_version=1,flow_mL_min=pulse["flow_ml_min"][0],inventory_mean_s=age,
        injected_salt_mg=pulse["volume_ml"]*.001*pulse["feed_NaCl_mg_L"],
        injection_variance_s2=width**2/12,ideal_collection_variance_s2=interval**2/12,
        pulse_and_collection_alpha_equivalent_mm=length*(width**2+interval**2)/(24*age**2),
        water=PROTOCOL["water"],checks=PROTOCOL["checks"],python_version=platform.python_version(),
        julia_version=read_json(ROOT/".local/julia_runtime.json")["julia_version"],openlb_version="1.9r0",
        native_solver="D3Q27, stock TRT, Kupershtokh forcing and bounce-back",
        source_doi=EVIDENCE["doi"],source_sha256={p.name:sha(p) for p in
          [HERE/"packing.jl",HERE/"run.py",HERE/"first_passage.jl",HERE/"protocol.toml",HERE/"report.py"]},beds={})
    build = read_json(native.BUILD/"build_record/build.json")
    summary["native_build"] = dict(source_sha256=sha(native.APP),binary_sha256=build["binary_sha256"],
        compiler=build["compiler"].splitlines()[0],vendor_integrity=build["vendor_integrity"])
    for path in (ROOT/"uv.lock", ROOT/"julia/Manifest.toml", ROOT/"literature/sources/berkowitz2009/evidence.toml",
                 native.HERE/"balance.jl", ROOT/"julia/adsorption/PorousFaceTransport.jl",
                 ROOT/"examples/adsorption/path_exchange/track.jl"):
        summary["source_sha256"][str(path.relative_to(ROOT))] = sha(path)
    density = PROTOCOL["density_screen"]
    concentration = pulse["feed_NaCl_mg_L"]*.001 / density["NaCl_molar_mass_kg_mol"]
    contrast = concentration*(density["NaCl_molar_mass_kg_mol"] -
        PROTOCOL["water"]["density_kg_m3"]*density["NaCl_limiting_apparent_molar_volume_m3_mol"])
    summary["undiluted_feed_density_contrast_kg_m3"] = contrast
    for beads,grids in (("fine",[6,8,10]),("coarse",[6,8,12,16,20])):
        period = config[beads+"_period_in_beads"]
        grid_cases = []
        for g in grids:
            record,arrivals = case(beads,period,base_seed,g)
            grid_cases.append(record)
        if len({r["packing_sha256"] for r in grid_cases}) != 1:
            raise ValueError("Mesh comparison changed the sphere coordinates")
        primary = grid_cases[-1]
        checks = {}
        for key, args in (("second_packing",(period,other_seed,8,1.,1.)),
                          ("double_period",(2*period,base_seed,8,1.,1.)),
                          ("half_timestep",(period,base_seed,8,.5,1.)),
                          ("higher_diffusivity",(period,base_seed,8,1.,1.05))):
            checks[key] = case(beads,*args)[0]
        for key in ("half_timestep","higher_diffusivity"):
            if checks[key]["native_field_sha256"] != grid_cases[1]["native_field_sha256"]:
                raise ValueError(f"{key} requires the identical native field")
        control,_ = case(beads,period,base_seed,8,low_compressibility=True)
        if control["packing_sha256"]!=grid_cases[1]["packing_sha256"] or control["tau_even"]!=.60:
            raise ValueError("Carrier time-scale control changed the physical packing")
        checks["lower_lattice_compressibility"] = control
        def alpha(r):
            return r["exit"][-1]["alpha_mm"]
        changes = {k:alpha(v)/alpha(grid_cases[1])-1 for k,v in checks.items()}
        changes["last_grid_pair"] = alpha(grid_cases[-2])/alpha(primary)-1
        changes["from_100_to_200_mm"] = primary["exit"][-2]["alpha_mm"]/alpha(primary)-1
        # Include the last numerical arrival and the broader comparison tail.
        horizon = max(arrivals.max()+width+interval,330*60+interval)
        while 1-reference_cdf(horizon-width,age,length,primary["reported_alpha_mm"]) > 1e-11:
            horizon *= 1.25
        edges = np.arange(0,horizon+interval,interval)
        concentration = finite_pulse_bins(arrivals,edges,width)
        recovery = float(np.dot(concentration,np.diff(edges))/width)
        reference = reference_bins(edges,width,age,length,primary["reported_alpha_mm"])
        reference_recovery = float(np.dot(reference,np.diff(edges))/width)
        if abs(reference_recovery-1)>1e-10:
            raise ValueError("Reference collection window truncates the published-coefficient response")
        if abs(recovery-1)>1e-12:
            raise ValueError("Finite-pulse collection did not preserve the dose")
        limits = PROTOCOL["checks"]
        screens = dict(
            throughput=abs(primary["flow_mL_min"]/pulse["flow_ml_min"][0]-1)<=limits["throughput_relative"],
            native_to_face_correction=primary["face_checks"]["face_rms_correction_relative"]<=limits["face_rms_correction_relative"],
            connected_pore_volume=primary["corner_only_fraction"]<=.001,
            mean_transit_time=abs(primary["mean_time_relative_error"])<=limits["mean_time_relative"],
            dose_recovery=abs(recovery-1)<=1e-12,
            published_alpha=primary["within_literature_screen"],
            **{key:abs(value)<=limits["relative_change_in_alpha"] for key,value in changes.items()})
        summary["beds"][beads] = dict(mesh=grid_cases,sensitivity=checks,
            relative_changes=changes,screens=screens,measured_length_mm=length,
            paper_implied_impulse_sd_min=age*math.sqrt(2*primary["reported_alpha_mm"]/length)/60,
            max_vertical_buoyancy_to_axial_driving_gradient=9.80665*contrast/primary["pressure_gradient_Pa_m"],
            predicted_impulse_sd_min=float(arrivals.std(ddof=1)/60),pulse_mass_recovery=recovery,
            five_minute_dose_at_achieved_flow_mg=primary["flow_mL_min"]*width/60*pulse["feed_NaCl_mg_L"]/1000,
            pulse_bin_mid_min=((edges[:-1]+edges[1:])/120).tolist(),
            pulse_C_over_feed=concentration.tolist(),
            ADE_from_published_alpha=reference.tolist(),ADE_integrated_fraction=reference_recovery)
    fine = summary["beds"]["fine"]["mesh"][-1]
    coarse = summary["beds"]["coarse"]["mesh"][-1]
    summary["coarse_to_fine_dispersivity_ratio"] = dict(
        predicted=coarse["exit"][-1]["alpha_mm"]/fine["exit"][-1]["alpha_mm"],
        published=coarse["reported_alpha_mm"]/fine["reported_alpha_mm"])
    return summary


def plot(summary, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1,2,figsize=(10,4.0),sharex=True,sharey=True,layout="constrained")
    for ax,(beads,bed) in zip(axes,summary["beds"].items(),strict=True):
        t = bed["pulse_bin_mid_min"]
        grid = bed["mesh"][-1]["grid"]
        ax.step(t,bed["pulse_C_over_feed"],where="mid",color="#176c79",lw=2,label=f"TRT flow + tracer, {grid} cells/bead")
        ax.plot(t,bed["ADE_from_published_alpha"],ls="--",color="#ac502b",lw=2,
                label="From published dispersivity")
        ax.set_title(f"{beads.capitalize()} bed  ·  {bed['mesh'][-1]['diameter_m']*1000:g} mm beads")
        ax.set(xlabel="Minutes after start of salt injection",xlim=(0,330))
        ax.grid(axis="y",alpha=.2)
        ax.spines[["top","right"]].set_visible(False)
        ax.legend(frameon=False,fontsize=8)
    axes[0].set_ylabel("Outlet / feed salt concentration")
    fig.savefig(path,dpi=165)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,default=OUTPUT/"comparison.json")
    parser.add_argument("--plot",type=Path,default=OUTPUT/"outlet_pulses.png")
    args = parser.parse_args()
    summary = summarize()
    write_json(args.output,summary)
    plot(summary,args.plot)
    for name,b in summary["beds"].items():
        r=b["mesh"][-1]
        print(f"{name}: {r['exit'][-1]['alpha_mm']:.3f} vs {r['reported_alpha_mm']} mm")


if __name__=="__main__":
    main()
