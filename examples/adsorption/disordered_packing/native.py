"""Run the stock TRT fields in separate external build and result directories."""
import argparse
from datetime import datetime, timezone
import os
from pathlib import Path
import shutil
import subprocess
from time import perf_counter
from uuid import uuid4

import numpy as np

from examples.adsorption.porous_face_handoff.run import archived
from examples.adsorption.repeated_geometry.run import OUTPUT as ORDERED_OUTPUT
from examples.adsorption.transport_root_cause.hydraulics.stock_omp import build_case_omp
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT

from .packing import DIAMETER_M, PERIODS, SEED, SIDE_M, SWEEPS, describe, disordered, ordered

HERE = Path(__file__).resolve().parent
OUTPUT = ROOT/"results/adsorption/disordered_packing"
BUILD = ROOT/".local/adsorption/disordered_packing"
ORDERED = ORDERED_OUTPUT/"runs/20261112T153358Z-4f1e7011"
ORDERED_SHA = "0cb54c71bbc5247de8c4b67836b0eeb34dd811949b498d985d2d22ffd0e547be"


def build():
    return build_case_omp(HERE/"openlb", BUILD/"build", BUILD/"build_record")


def prepare():
    path = ORDERED/"summary.json"
    if sha(path) != ORDERED_SHA or not read_json(path)["analysis"]["longer_holdout_passed"]:
        raise ValueError("Pinned eight-cell ordered tracer evidence changed")
    plan = read_json(ORDERED/"frozen_prediction.json")
    if sha(ORDERED/"frozen_prediction.json") != read_json(path)["frozen_prediction_sha256"]:
        raise ValueError("Pinned ordered tracer plan changed")
    info, _, _ = archived(96)
    np.testing.assert_allclose(plan["period_length_m"],SIDE_M,rtol=1e-14,atol=0.)
    velocity = info["carrier_directional_flow_m3_s"][0]/SIDE_M**2/info["porosity"]
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")+"-"+uuid4().hex[:8]
    base = OUTPUT/"runs"/stamp
    base.mkdir(parents=True)
    centers, acceptance = disordered()
    controls = ordered(1)
    for label, points in (("disordered", centers), ("ordered", controls)):
        np.savetxt(base/(label+"_centres_m.txt"),points*SIDE_M,fmt="%.17g")
    comparable = dict(disordered=describe(centers,PERIODS), ordered=describe(ordered(PERIODS),PERIODS))
    np.testing.assert_equal(comparable["ordered"]["solid_volume_m3"],comparable["disordered"]["solid_volume_m3"])
    p = dict(geometry=comparable, seed=SEED, sweeps=SWEEPS, packing_acceptance=acceptance,
             period_length_m=SIDE_M, particle_diameter_m=DIAMETER_M, dry_carbon_mass_kg=plan["dry_carbon_mass_kg"],
             ordered_reference_sha256=ORDERED_SHA, reference_interstitial_m_s=velocity,
             velocity_scale=info["velocity_scale"], diameter_cases=plan["predictions"],
             ordered_reference_voxel_porosity=info["porosity"],
             centres_sha256={f:sha(base/f) for f in ("disordered_centres_m.txt","ordered_centres_m.txt")})
    write_json(base/"plan.json",p)
    write_json(OUTPUT/"latest_plan.json",dict(path=str(base.relative_to(ROOT))))
    print(base,flush=True)
    return base


def native(base, label, n, threads, pilot=False):
    plan = read_json(base/"plan.json")
    building = read_json(BUILD/"build_record/build.json")
    for name, digest in building["source_sha256"].items():
        if sha(HERE/"openlb"/name) != digest:
            raise ValueError("Stock TRT application source changed since external build")
    exe = Path(building["binary"])
    if sha(exe) != building["binary_sha256"]:
        raise ValueError("External binary changed")
    centers = base/(label+"_centres_m.txt")
    if sha(centers) != plan["centres_sha256"][centers.name]:
        raise ValueError("Frozen sphere centers changed")
    initializer = base/f"pilot_{label}_N{n}"
    folder = initializer if pilot else base/f"{label}_N{n}"
    folder.mkdir()
    shutil.copy2(BUILD/"build_record/build.json",folder/"build.json")
    tau = .8 if pilot else {48:.54,64:.553,96:.58}[n]
    warm = ()
    gradient = (2277.,0.,0.)
    if not pilot and (initializer/"raw.json").exists():
        old = read_json(initializer/"raw.json")
        gradient = old["gradient_Pa_m"]
        warm = (initializer/"field.f64",((tau-.5)/(old["tau_even"]-.5))**2)
    cmd = [str(exe),*map(str,(n,PERIODS if label=="disordered" else 1,SIDE_M,tau,
           plan["reference_interstitial_m_s"]/plan["velocity_scale"],*gradient,centers,*warm))]
    write_json(folder/"request.json",dict(command=cmd,omp_threads=threads,
               warm_field_sha256=sha(warm[0]) if warm else None))
    env = {**os.environ,"OMP_NUM_THREADS":str(threads),"OMP_DYNAMIC":"FALSE"}
    print("Stock TRT start",folder,flush=True)
    start = perf_counter()
    with (folder/"native.log").open("w") as log:
        result = subprocess.run(cmd,cwd=folder,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=14400)
    write_json(folder/"execution.json",dict(returncode=result.returncode,wall_s=perf_counter()-start))
    if result.returncode:
        raise RuntimeError(f"Stock TRT failed. See {folder/'native.log'}")
    print("Stock TRT completed",folder,read_json(folder/"raw.json"),flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action",required=True)
    sub.add_parser("prepare")
    sub.add_parser("build")
    run = sub.add_parser("run")
    run.add_argument("folder",type=Path)
    run.add_argument("label",choices=("ordered","disordered"))
    run.add_argument("grid",type=int,choices=(48,64,96))
    run.add_argument("--threads",type=int,default=10)
    run.add_argument("--pilot",action="store_true")
    args = parser.parse_args()
    if args.action == "build":
        print(build()["binary"],flush=True)
    elif args.action == "prepare":
        prepare()
    else:
        native(args.folder.resolve(),args.label,args.grid,args.threads,args.pilot)


if __name__ == "__main__":
    main()
