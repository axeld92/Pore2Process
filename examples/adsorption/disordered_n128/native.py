"""Run the stock OpenLB N128 refinement from the verified N96 field."""
import argparse
from datetime import datetime, timezone
import os
from pathlib import Path
import shutil
import subprocess
from time import perf_counter
from uuid import uuid4

import numpy as np

from examples.adsorption.disordered_packing.mesh import fluid_mask, load_field
from examples.adsorption.disordered_packing.native import BUILD, HERE as PACKING
from examples.adsorption.disordered_packing.packing import PERIODS, SIDE_M
from examples.adsorption.disordered_refinement.native import PARENT as N64, interpolate
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT

HERE = Path(__file__).resolve().parent
OUTPUT = ROOT/"results/adsorption/disordered_n128"
PARENT = ROOT/"results/adsorption/disordered_refinement/runs/20261112T171431Z-193be93e"
PARENT_SHA = "4cf9d60ca9471036681fe62d7f166414c470a3e7445e6ac9e7b476a8f697dd0f"
GRID = 128


def parent():
    path = PARENT/"summary.json"
    if sha(path)!=PARENT_SHA:
        raise ValueError("Pinned N96 refinement changed")
    record = read_json(path)
    paths = ("plan.json","disordered_centres_m.txt","disordered_N96/raw.json",
             "disordered_N96/field.f64","disordered_N96/tracer_summary.json",
             "disordered_N96/tracer_plan.json","disordered_N96/compact/summary.json",
             "projection_N64_parity/summary.json")
    for name in paths:
        if sha(PARENT/name)!=record["artifacts_sha256"][name]:
            raise ValueError(f"Pinned N96 input changed: {name}")
    return record


def initial_field(base):
    previous,ijk,_,info = load_field(PARENT,"disordered",96)
    field = np.memmap(previous/"field.f64",mode="r",dtype="<f8",shape=(len(ijk),7))
    shape = (PERIODS*GRID,GRID,GRID)
    mask = fluid_mask(np.loadtxt(base/"disordered_centres_m.txt"),shape,SIDE_M/GRID)
    fine = np.argwhere(mask)
    target = base/"initial_N128.f64"
    result = np.memmap(target,mode="w+",dtype="<f8",shape=(len(fine),7))
    result[:,:3] = fine
    for j in range(4):
        moments = np.zeros(info["shape"])
        baseline = 1000. if j==3 else 0.
        moments[tuple(ijk.T)] = field[:,3+j]-baseline
        result[:,3+j] = interpolate(moments,shape)[mask]+baseline
    result.flush()
    return target


def run(base,threads):
    config = read_json(base/"plan.json")
    record = parent()
    fine = record["analysis"]["native"]["raw"]
    coarse = read_json(N64/"disordered_N64/raw.json")
    expected = read_json(N64/"summary.json")["analysis"]["hydraulics"]["disordered_N64"]["info"]["raw"]
    if coarse!=expected:
        raise ValueError("N64 gradient reference changed")
    building = read_json(BUILD/"build_record/build.json")
    for name,digest in building["source_sha256"].items():
        if sha(PACKING/"openlb"/name)!=digest:
            raise ValueError("External OpenLB application source changed")
    binary = Path(building["binary"])
    prior_build = read_json(PARENT/"disordered_N96/build.json")
    if sha(binary)!=building["binary_sha256"] or building["source_sha256"]!=prior_build["source_sha256"]:
        raise ValueError("Stock TRT executable or its application changed")
    folder = base/f"disordered_N{GRID}"
    folder.mkdir()
    warm = initial_field(base)
    tau = .5+(fine["tau_even"]-.5)*GRID/fine["resolution"]
    gradient = np.array(fine["gradient_Pa_m"])
    slope = (gradient-coarse["gradient_Pa_m"])/(1/fine["resolution"]-1/coarse["resolution"])
    gradient += slope*(1/GRID-1/fine["resolution"])
    centres = base/"disordered_centres_m.txt"
    cmd = [str(binary),*map(str,(GRID,PERIODS,SIDE_M,tau,
        config["reference_interstitial_m_s"]/config["velocity_scale"],*gradient,centres,warm,1.))]
    shutil.copy2(BUILD/"build_record/build.json",folder/"build.json")
    write_json(folder/"request.json",dict(command=cmd,omp_threads=threads,parent_sha256=PARENT_SHA,
               warm_field_sha256=sha(warm),centres_sha256=sha(centres),
               protocol_sha256=sha(HERE/"protocol.md"),native_source_sha256=sha(Path(__file__))))
    env = {**os.environ,"OMP_NUM_THREADS":str(threads),"OMP_DYNAMIC":"FALSE"}
    print("Stock TRT N128 start",folder,flush=True)
    start = perf_counter()
    with (folder/"native.log").open("w") as log:
        result = subprocess.run(cmd,cwd=folder,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=21600)
    write_json(folder/"execution.json",dict(returncode=result.returncode,wall_s=perf_counter()-start))
    if result.returncode:
        raise RuntimeError(f"Stock TRT N128 failed. See {folder/'native.log'}")
    _,_,_,info = load_field(base,"disordered",GRID)
    write_json(folder/"verified_field.json",info)
    print("Stock TRT N128 accepted",info,flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--threads",type=int,default=10)
    args = p.parse_args()
    parent()
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")+"-"+uuid4().hex[:8]
    base = OUTPUT/"runs"/stamp
    base.mkdir(parents=True)
    for name in ("plan.json","disordered_centres_m.txt"):
        shutil.copy2(PARENT/name,base/name)
    write_json(OUTPUT/"latest_plan.json",dict(path=str(base.relative_to(ROOT))))
    print(base,flush=True)
    run(base,args.threads)


if __name__=="__main__":
    main()
