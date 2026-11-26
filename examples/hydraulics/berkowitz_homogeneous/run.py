"""Stock TRT flow, conservative face projection and salt first-passage transport."""
import argparse
import math
import os
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory
import tomllib

import numpy as np

from examples.adsorption.breakthrough_hypotheses.hydrodynamics import handoff, native
from examples.adsorption.disordered_packing.mesh import shared_faces
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT

HERE = Path(__file__).resolve().parent
OUTPUT = ROOT / "results/hydraulics/berkowitz_homogeneous"
PROTOCOL = tomllib.loads((HERE / "protocol.toml").read_text())
EVIDENCE = tomllib.loads((ROOT / "literature/sources/berkowitz2009/evidence.toml").read_text())


def run_logged(command, folder, logfile, threads):
    env = {**os.environ, "OMP_NUM_THREADS": str(threads), "OMP_DYNAMIC": "FALSE",
           "OPENBLAS_NUM_THREADS": "1"}
    with (folder / logfile).open("w") as log:
        subprocess.run(list(map(str, command)), cwd=folder, env=env, stdout=log,
                       stderr=subprocess.STDOUT, check=True)


def settings(beads, period, seed, grid):
    column, pulse, water = EVIDENCE["column"], EVIDENCE["pulse"], PROTOCOL["water"]
    dp = column[f"{beads}_bead_diameter_mm"] * .001
    R = column["inner_diameter_mm"] * .0005
    q = pulse["flow_ml_min"][0] * 1e-6 / 60
    eps = column["porosity"]
    speed = q / (math.pi*R**2*eps)
    nu = water["viscosity_Pa_s"]/water["density_kg_m3"]
    D = (water["NaCl_D25_m2_s"] * (273.15+pulse["temperature_C"]) /
         (273.15+water["temperature_reference_C"]) *
         water["water_viscosity_25_Pa_s"]/water["viscosity_Pa_s"])
    nx = round(period * grid)
    h = period * dp / nx
    nt = math.ceil(2*R/h) + 4
    return dict(beads=beads, seed=seed, requested_cells_per_bead=grid,
                period_beads=period, bead_diameter_m=dp, tube_radius_m=R,
                pore_fraction=eps, flow_m3_s=q, pore_speed_m_s=speed,
                pore_Re=speed*dp/nu, pore_Pe=speed*dp/D, nu_m2_s=nu, D_m2_s=D,
                shape=[nx,nt,nt], dx_m=h)


def centres_path(p):
    name = f"{p['beads']}_{p['period_beads']:g}_{p['seed']}.tsv"
    return OUTPUT / "packing" / name


def prepare(p, folder, threads):
    target = centres_path(p)
    target.parent.mkdir(parents=True, exist_ok=True)
    candidate = folder / "packing_candidate.tsv"
    # Recreate from source, including on restart, so a stale geometry cannot
    # be silently reused. Use the candidate only after packing has succeeded.
    run_logged([*handoff.julia(), HERE/"packing.jl", p["tube_radius_m"]/p["bead_diameter_m"],
                p["period_beads"], p["pore_fraction"], p["seed"], candidate],
               folder, "packing.log", threads)
    if target.exists() and sha(target) != sha(candidate):
        raise ValueError(f"Packing differs from saved {target}. Use a fresh output root.")
    candidate.replace(target)
    return target


def carrier(p, folder, threads):
    pack = prepare(p, folder, threads)
    buildfile = native.BUILD / "build_record/build.json"
    if not buildfile.is_file():
        native.build()
    build = read_json(buildfile)
    binary = Path(build["binary"])
    if sha(binary) != build["binary_sha256"] or sha(native.APP) != build["source_sha256"][native.APP.name]:
        raise ValueError("Native executable/source mismatch. Rebuild the external sphere case.")
    dp, R, q = p["bead_diameter_m"], p["tube_radius_m"], p["flow_m3_s"]
    nx, nt, _ = p["shape"]
    # Reuse the verified sphere-case executable by exact Navier–Stokes
    # similarity. Starred coordinates are in its existing external case units.
    s = native.p.DIAMETER / dp
    nu_star = native.p.VISCOSITY/native.p.DENSITY
    u_scale = (nu_star/p["nu_m2_s"])/s
    q_star = q * s**2 * u_scale
    points = np.loadtxt(pack)
    shifted = points*native.p.DIAMETER + [0,nt*p["dx_m"]*s/2,nt*p["dx_m"]*s/2]
    np.savetxt(folder/"centres_native_m.txt",shifted,fmt="%.17g")
    eps = p["pore_fraction"]
    superficial = q_star/(math.pi*(R*s)**2)
    gradient = 150*(1-eps)**2/eps**3 * native.p.VISCOSITY*superficial / native.p.DIAMETER**2
    tau = p.get("tau_even",.85)
    warm = []
    if tau != .85:
        baseline = folder.with_name(folder.name.removesuffix("_tau060"))
        old = read_json(baseline/"carrier.json")
        if old["packing_sha256"] != sha(pack) or old["shape"] != p["shape"] or sha(baseline/"field.f64") != old["native_field_sha256"]:
            raise ValueError("Low-compressibility control requires the same converged baseline field")
        raw = read_json(baseline/"raw.json")
        gradient = raw["gradient_Pa_m"][0]
        warm = [baseline/"field.f64",((tau-.5)/(raw["tau_even"]-.5))**2]
        p["initial_field_sha256"] = old["native_field_sha256"]
    command = [binary, nx, nt, p["period_beads"]*native.p.DIAMETER, R*s, tau, 3/16,
               q_star, gradient, folder/"centres_native_m.txt",*warm]
    # Stock OpenLB writes a construction frame after each sphere. All
    # scientific outputs come from the final geometry and converged field.
    with TemporaryDirectory(dir=folder,prefix=".openlb-") as scratch:
        work = Path(scratch)
        try:
            run_logged(command, work, folder/"native.log", threads)
        except subprocess.CalledProcessError:
            if (work/"history.csv").exists():
                (work/"history.csv").replace(folder/"failed_history.csv")
            raise
        for filename in ("raw.json","field.f64","history.csv"):
            (work/filename).replace(folder/filename)
    raw, field = native.field(folder)
    if (raw["collision"], raw["converged"], raw["shape"]) != ("TRT",True,p["shape"]):
        raise ValueError("Unexpected native run")
    if not np.isclose(raw["dx_m"],p["dx_m"]*s,rtol=1e-12,atol=0) or raw["max_Mach"] >= .05:
        raise ValueError("Native units or Mach check failed")
    mean_rate = float(field[:,3].sum()/u_scale * p["dx_m"]**2/nx)
    if abs(mean_rate/q-1) > PROTOCOL["checks"]["throughput_relative"]:
        raise ValueError("Native throughput outside tolerance")
    grain_porosity = 1-len(points)*dp**2/(6*R**2*p["period_beads"])
    voxel_porosity = len(field)*p["dx_m"]**2/(nx*math.pi*R**2)
    p.update(native_length_per_physical_length=s, native_speed_per_physical_speed=u_scale,
             native_Q_per_physical_Q=s*s*u_scale, achieved_flow_m3_s=mean_rate,
             native_max_Mach=raw["max_Mach"], converged_steps=raw["steps"], sphere_count=len(points),
             geometric_porosity=grain_porosity, voxel_porosity=voxel_porosity,
             native_field_sha256=sha(folder/"field.f64"), packing_sha256=sha(pack),
             native_source_sha256=sha(native.APP), binary_sha256=sha(binary),
             native_command=list(map(str,command)))
    write_json(folder/"carrier.json",p)


def balance(folder, threads):
    p = read_json(folder/"carrier.json")
    if sha(folder/"field.f64") != p["native_field_sha256"]:
        raise ValueError("Native field differs from recorded run")
    _, field = native.field(folder)
    shape, h = p["shape"], p["dx_m"]
    ijk = field[:,:3].astype(np.int32)
    vel = field[:,3:6] / p["native_speed_per_physical_speed"]
    mask = np.zeros(shape,dtype=np.uint8)
    mask[tuple(ijk.T)] = 1
    keep = handoff.connected_pore_component(ijk,mask,vel)
    excluded_flux = float(vel[~keep,0].sum()/vel[:,0].sum())
    ijk, vel = ijk[keep], vel[keep]
    faces = shared_faces(ijk,vel,shape,h)
    work = folder/"faces"
    work.mkdir(exist_ok=True)
    for name,a in faces.items():
        a.tofile(work/name)
    target = p["achieved_flow_m3_s"]
    run_logged([*handoff.julia(), native.HERE/"balance.jl",work,len(ijk),shape[0]*target],
               work,"balance.log",threads)
    q = np.fromfile(work/"balanced.f64",dtype="<f8")
    checks = handoff.conservation(ijk,dict(shape=shape,directional_flow_m3_s=[target,0.,0.]),faces,q)
    left, right = faces["left.i32"], faces["right.i32"]
    axis = faces["axis.u8"] - 1
    before = np.zeros((len(ijk),3),dtype="<i4")
    after = np.zeros_like(before)
    incoming = np.zeros(before.shape,dtype="<f8")
    outgoing = np.zeros_like(incoming)
    before[right-1,axis] = left
    after[left-1,axis] = right
    incoming[right-1,axis] = outgoing[left-1,axis] = q/h**3
    for name,matrix in (("before.i32",before),("after.i32",after),("low.f64",incoming),("high.f64",outgoing)):
        matrix.tofile(work/name)
    inlet = np.flatnonzero((ijk[:,0]==0) & (incoming[:,0]>0))
    (inlet+1).astype("<i4").tofile(work/"inlet.i32")
    reverse = -np.minimum(incoming[ijk[:,0]==0,0],0).sum()*h**3/target
    write_json(folder/"projection.json",dict(checks=checks,fluid_cells=len(ijk),
        corner_only_fraction=1-float(keep.mean()), excluded_native_axial_flow_fraction=excluded_flux,
        qualifies_connected_void=(1-float(keep.mean()) <= .001), entry_reverse_flux_fraction=float(reverse),
        connected_porosity=p["voxel_porosity"]*float(keep.mean()),carrier_sha256=sha(folder/"carrier.json")))


def transport(folder, threads, factor, markers, diffusion_factor):
    p = read_json(folder/"carrier.json")
    face = read_json(folder/"projection.json")
    conf = PROTOCOL["transport"]
    if face["carrier_sha256"] != sha(folder/"carrier.json"):
        raise ValueError("Projected flow belongs to another carrier")
    h, D = p["dx_m"], p["D_m2_s"]*diffusion_factor
    speed = p["achieved_flow_m3_s"]/(math.pi*p["tube_radius_m"]**2*face["connected_porosity"])
    dt = factor*min((conf["max_rms_diffusion_step_cells"]*h)**2/(2*D),
                    conf["max_mean_drift_step_cells"]*h/speed)
    levels = np.rint(np.array(conf["lengths_mm"])*.001/h).astype(int)
    stop = conf["end_mean_transit_times"]*levels[-1]*h/speed
    prefix = f"tracer_m{markers}_dt{factor:g}_D{diffusion_factor:g}"
    path = folder / (prefix+".f64")
    run_logged([*handoff.julia(threads), HERE/"first_passage.jl",folder/"faces",path,h,D,dt,stop,
                markers,conf["seed"],",".join(map(str,levels))],folder,prefix+".log",threads)
    times = np.fromfile(path,dtype="<f8").reshape(markers,len(levels))
    if not np.all(np.isfinite(times)) or np.any(np.diff(times,axis=1)<=0):
        raise ValueError(f"Unrecovered or out-of-order tracer in {path}. Extend the observation window.")
    records = []
    for t, level in zip(times.T,levels,strict=True):
        length = level*h
        mean, var = float(t.mean()),float(t.var(ddof=1))
        alpha = length*var/(2*mean**2)
        batches = t.reshape(16,-1)
        estimates = length*batches.var(axis=1,ddof=1)/(2*batches.mean(axis=1)**2)
        records.append(dict(length_m=length,mean_s=mean,variance_s2=var,
            nominal_mean_s=length/speed,alpha_mm=1000*alpha,alpha_mc_se_mm=1000*float(estimates.std(ddof=1))/4,
            arrivals=len(t),clock_quantiles_s=np.quantile(t,[.1,.5,.9]).tolist()))
    write_json(folder/(prefix+".json"),dict(records=records,markers=markers,
        seed=conf["seed"],dt_s=dt,dt_factor=factor,D_m2_s=D,diffusion_factor=diffusion_factor,
        stop_s=stop,arrivals_sha256=sha(path),projection_sha256=sha(folder/"projection.json")))
    print(folder.name,prefix,"alpha at 200 mm",records[-1]["alpha_mm"],flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage",choices=("all","flow","faces","tracer"))
    parser.add_argument("beads",choices=("fine","coarse"))
    parser.add_argument("--grid",type=int,default=8)
    parser.add_argument("--seed",type=int,default=PROTOCOL["packing"]["baseline_seed"])
    parser.add_argument("--period",type=float)
    parser.add_argument("--threads",type=int,default=6)
    parser.add_argument("--dt-factor",type=float,default=1)
    parser.add_argument("--diffusion-factor",type=float,default=1)
    parser.add_argument("--markers",type=int,default=PROTOCOL["transport"]["markers"])
    parser.add_argument("--low-compressibility",action="store_true",
                        help="Repeat the matched carrier at TRT relaxation time 0.60 using its converged baseline")
    args = parser.parse_args()
    if (args.grid < 6 or args.threads < 1 or args.markers < 32 or args.markers%16 or
        any(not math.isfinite(v) or v <= 0 for v in (args.dt_factor,args.diffusion_factor)) or
        (args.period is not None and (not math.isfinite(args.period) or args.period < 4))):
        parser.error("Use grid ≥6, positive threads and factors, and markers ≥32 divisible by 16.")
    period = args.period or PROTOCOL["packing"][f"{args.beads}_period_in_beads"]
    p = settings(args.beads,period,args.seed,args.grid)
    p["tau_even"] = .60 if args.low_compressibility else .85
    suffix = "_tau060" if args.low_compressibility else ""
    folder = OUTPUT/f"{args.beads}_L{period:g}_s{args.seed}_N{args.grid}{suffix}"
    folder.mkdir(parents=True,exist_ok=True)
    if args.stage in ("all","flow"):
        carrier(p,folder,args.threads)
    if args.stage in ("all","faces"):
        balance(folder,args.threads)
    if args.stage in ("all","tracer"):
        transport(folder,args.threads,args.dt_factor,args.markers,args.diffusion_factor)


if __name__ == "__main__":
    main()
