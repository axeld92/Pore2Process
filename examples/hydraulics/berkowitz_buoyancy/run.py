"""Launch the matched finite-dose runs. No fitted outlet inputs enter this module."""
import argparse
import fcntl
import json
import math
import os
from pathlib import Path
import subprocess

import numpy as np
import tomllib

from examples.adsorption.breakthrough_hypotheses.hydrodynamics import handoff
from examples.hydraulics.berkowitz_homogeneous import run as passive
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT
from . import build

HERE = Path(__file__).resolve().parent
OUT = ROOT/"results/hydraulics/berkowitz_buoyancy"
PROTOCOL = tomllib.loads((HERE/"protocol.toml").read_text())
DEFAULT = PROTOCOL["simulation"]


def inputs(options):
    saved = passive.OUTPUT/f"coarse_L7_s31_N{options.grid}"
    if options.low_mach:
        saved = saved.with_name(saved.name+"_tau060")
    if not (saved/"projection.json").is_file():
        raise ValueError(f"Prepare the homogeneous-bed flow and faces at {saved}. See the reproduction commands in README.md.")
    p = read_json(saved/"carrier.json")
    project = read_json(saved/"projection.json")
    raw = read_json(saved/"raw.json")
    if (sha(saved/"field.f64") != p["native_field_sha256"] or
            sha(saved/"carrier.json") != project["carrier_sha256"] or
            sha(passive.centres_path(p)) != p["packing_sha256"]):
        raise ValueError("Saved field, packing or projected face graph has changed")
    if not build.BINARY.is_file():
        build.build()
    record = read_json(build.OUT/"build.json")
    if record["binary_sha256"] != sha(build.BINARY) or record["source_sha256"] != sha(build.HERE/"openlb/salt3d.cpp"):
        raise ValueError("Rebuild the external native salt case after a change")
    field = np.memmap(saved/"field.f64", mode="r", dtype="<f8").reshape(-1,7)
    ijk = field[:,:3].astype("int32")
    mask = np.zeros(p["shape"],dtype="uint8")
    mask[tuple(ijk.T)] = 1
    keep = handoff.connected_pore_component(ijk,mask,field[:,3:6])
    if keep.sum() != project["fluid_cells"]:
        raise ValueError("Connected pore indexing differs from the saved projection")
    water = passive.PROTOCOL["water"]
    density = passive.PROTOCOL["density_screen"]
    pulse = passive.EVIDENCE["pulse"]
    salt = pulse["feed_NaCl_mg_L"]*0.001
    drho = salt*(1-water["density_kg_m3"]*density["NaCl_limiting_apparent_molar_volume_m3_mol"]/
                 density["NaCl_molar_mass_kg_mol"])
    g = -options.gravity*DEFAULT["gravity_m_s2"]*drho/water["density_kg_m3"]
    axial = raw["gradient_Pa_m"][0]/1000*p["native_length_per_physical_length"]/p["native_speed_per_physical_speed"]**2
    folder = OUT/options.name
    folder.mkdir(parents=True, exist_ok=True)
    (np.flatnonzero(keep)+1).astype("<i4").tofile(folder/"active.i32")
    native = [build.BINARY,*p["shape"][:2],options.repeats,p["dx_m"],p["nu_m2_s"],water["density_kg_m3"],
              p["native_speed_per_physical_speed"],axial,raw["tau_even"],saved/"field.f64"]
    (folder/"native_args.txt").write_text("\n".join(map(str,native))+"\n")
    q = p["flow_m3_s"]
    dose_m3 = pulse["volume_ml"]*1e-6
    settings = dict(saved=str(saved),nx=p["shape"][0],nt=p["shape"][1],repeats=options.repeats,
                    h=p["dx_m"],diffusivity=p["D_m2_s"],q=q,gravity_accel=g,period_m=0.028,
                    bounded=int(not options.unbounded),tau_even=raw["tau_even"],kernel_m=options.kernel,markers=options.markers,seed=options.seed,source_m3=dose_m3,
                    feed_s=dose_m3/q,stop_s=DEFAULT["end_s"],exit_m=DEFAULT["observation_m"],
                    dt=options.dt,coupling_s=options.interval,sweeps=options.sweeps,warm=DEFAULT["warm_sweeps"],
                    salt_mass_mg=salt*dose_m3*1e6,feed_kg_m3=salt,delta_rho_feed_kg_m3=drho,
                    native_speed_scale=p["native_speed_per_physical_speed"],native_rows=len(ijk),gravity_multiplier=options.gravity,geometry_sha256=p["packing_sha256"],
                    water_sha256=p["native_field_sha256"],binary_sha256=record["binary_sha256"],
                    native_source_sha256=record["source_sha256"],stepper_sha256=sha(HERE/"coupled.jl"),
                    faces_sha256={f.name:sha(f) for f in sorted((saved/"faces").iterdir()) if f.suffix in (".f64",".u8",".i32")},
                    bulk_Q_before_projection_m3_s=p["achieved_flow_m3_s"])
    write_json(folder/"inputs.json",settings)
    (folder/"runtime.toml").write_text("".join(f"{k} = {json.dumps(v)}\n" for k,v in settings.items() if isinstance(v,(float,int,str))))
    return folder


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("name")
    parser.add_argument("--unbounded",action="store_true",help="Reproduce the initial density-kernel sensitivity case")
    parser.add_argument("--low-mach",action=argparse.BooleanOptionalAction,default=DEFAULT["hydraulic_tau_even"]==.60,help="Use the archived tau=0.60 carrier for this geometry")
    parser.add_argument("--grid",type=int,default=DEFAULT["cells_per_bead"])
    parser.add_argument("--gravity",type=float,choices=(-1.,0.,1.),default=1.)
    for key,typ,val in (("repeats",int,DEFAULT["repeats"]),("markers",int,DEFAULT["markers"]),
                        ("seed",int,DEFAULT["seed"]),("kernel",float,DEFAULT["density_kernel_m"]),
                        ("dt",float,DEFAULT["tracer_step_s"]),("interval",float,DEFAULT["coupling_s"]),
                        ("sweeps",int,DEFAULT["native_sweeps"]),("threads",int,4)):
        parser.add_argument("--"+key,type=typ,default=val)
    a = parser.parse_args()
    positive=(a.grid,a.markers,a.sweeps,a.threads,a.kernel,a.dt,a.interval)
    if (Path(a.name).name!=a.name or a.name in ("",".","..") or
            not all(math.isfinite(x) and x>0 for x in positive) or a.repeats*0.028<=DEFAULT["observation_m"]):
        parser.error("Use a case name, positive control values, and a domain longer than the sampling distance")
    folder=OUT/a.name
    folder.mkdir(parents=True,exist_ok=True)
    with (folder/".run.lock").open("a") as lock:
        try:
            fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:
            parser.error(f"A simulation is already writing {folder}")
        inputs(a)
        env={**os.environ,"OMP_NUM_THREADS":str(a.threads),"OMP_DYNAMIC":"FALSE","OPENBLAS_NUM_THREADS":"1"}
        command=[*handoff.julia(a.threads),HERE/"coupled.jl"]
        with (folder/"run.log").open("w") as log:
            subprocess.run(command,cwd=folder,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)


if __name__=="__main__":
    main()
