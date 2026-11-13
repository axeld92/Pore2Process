"""Solve the N160 lattice with the archived stock TRT application."""
import argparse
from datetime import datetime, timezone
import os
from pathlib import Path
import shutil
import subprocess
from time import perf_counter
from uuid import uuid4

import numpy as np

from examples.adsorption.disordered_n128.native import parent as n96_archive
from examples.adsorption.disordered_packing.mesh import fluid_mask, load_field
from examples.adsorption.disordered_packing.native import BUILD, HERE as PACKING
from examples.adsorption.disordered_packing.packing import PERIODS, SIDE_M
from examples.adsorption.disordered_refinement.native import interpolate
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT

HERE = Path(__file__).resolve().parent
OUTPUT = ROOT/"results/adsorption/disordered_n160"
PARENT = ROOT/"results/adsorption/disordered_n128/runs/20261112T214029Z-59d60ca7"
PARENT_SHA = "e28d69bb9dac4441edcd5ba20f5486463a3502552d0832a16d828de957414117"
GRID = 160
PARENT_GRID = 128


def parent():
    path = PARENT/"summary.json"
    if sha(path) != PARENT_SHA:
        raise ValueError("Pinned N128 refinement changed")
    record = read_json(path)
    previous = PARENT/f"disordered_N{PARENT_GRID}"
    paths = [PARENT/name for name in ("plan.json","disordered_centres_m.txt")]
    paths += [previous/name for name in ("raw.json","field.f64","build.json","tracer_summary.json",
                                        "tracer_plan.json","compact/summary.json")]
    for path in paths:
        name = str(path.relative_to(PARENT))
        if sha(path) != record["artifacts_sha256"][name]:
            raise ValueError(f"Pinned N128 input changed: {name}")
    return record


def initial_field(base, centres):
    previous = PARENT/f"disordered_N{PARENT_GRID}"
    old_shape = (PERIODS*PARENT_GRID,PARENT_GRID,PARENT_GRID)
    field = np.memmap(previous/"field.f64",mode="r",dtype="<f8").reshape(-1,7)
    old_coords = tuple(field[:,d].astype(np.int32) for d in range(3))
    shape = (PERIODS*GRID,GRID,GRID)
    mask = fluid_mask(np.loadtxt(centres),shape,SIDE_M/GRID)
    indices = np.flatnonzero(mask)
    target = base/f"initial_N{GRID}.f64"
    result = np.memmap(target,mode="w+",dtype="<f8",shape=(len(indices),7))
    result[:,0] = indices//(GRID*GRID)
    result[:,1] = indices//GRID%GRID
    result[:,2] = indices%GRID
    for j in range(4):
        baseline = 1000. if j==3 else 0.
        moments = np.zeros(old_shape)
        moments[old_coords] = field[:,3+j]-baseline
        result[:,3+j] = interpolate(moments,shape)[mask]+baseline
    result.flush()
    return target


def run(base, record, threads):
    config = read_json(base/"plan.json")
    fine = record["analysis"]["native"]["raw"]
    coarse = n96_archive()["analysis"]["native"]["raw"]
    building = read_json(BUILD/"build_record/build.json")
    prior = read_json(PARENT/f"disordered_N{PARENT_GRID}/build.json")
    binary = Path(building["binary"])
    if sha(binary)!=building["binary_sha256"] or building["source_sha256"]!=prior["source_sha256"]:
        raise ValueError("Stock TRT executable or its application changed")
    for name,digest in building["source_sha256"].items():
        if sha(PACKING/"openlb"/name)!=digest:
            raise ValueError("External OpenLB application source changed")
    folder = base/f"disordered_N{GRID}"
    folder.mkdir()
    centres = base/"disordered_centres_m.txt"
    warm = initial_field(base, centres)
    tau = .5+(fine["tau_even"]-.5)*GRID/PARENT_GRID
    gradient = np.array(fine["gradient_Pa_m"])
    slope = (gradient-coarse["gradient_Pa_m"])/(1/PARENT_GRID-1/coarse["resolution"])
    gradient += slope*(1/GRID-1/PARENT_GRID)
    cmd = [str(binary),*map(str,(GRID,PERIODS,SIDE_M,tau,
        config["reference_interstitial_m_s"]/config["velocity_scale"],*gradient,centres,warm,1.))]
    shutil.copy2(BUILD/"build_record/build.json",folder/"build.json")
    write_json(folder/"request.json",dict(command=cmd,omp_threads=threads,parent_sha256=PARENT_SHA,
               warm_field_sha256=sha(warm),centres_sha256=sha(centres),
               protocol_sha256=sha(HERE/"protocol.md"),native_source_sha256=sha(Path(__file__))))
    env = {**os.environ,"OMP_NUM_THREADS":str(threads),"OMP_DYNAMIC":"FALSE"}
    print("Stock TRT N160 start",folder,flush=True)
    start = perf_counter()
    with (folder/"native.log").open("w") as log:
        result = subprocess.run(cmd,cwd=folder,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=21600)
    write_json(folder/"execution.json",dict(returncode=result.returncode,wall_s=perf_counter()-start))
    if result.returncode:
        raise RuntimeError(f"Stock TRT N160 failed. See {folder/'native.log'}")
    _,_,_,info = load_field(base,"disordered",GRID)
    write_json(folder/"verified_field.json",info)
    print("Stock TRT N160 accepted",info,flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--threads",type=int,default=10)
    args = parser.parse_args()
    record = parent()
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")+"-"+uuid4().hex[:8]
    base = OUTPUT/"runs"/stamp
    base.mkdir(parents=True)
    for name in ("plan.json","disordered_centres_m.txt"):
        shutil.copy2(PARENT/name,base/name)
    write_json(OUTPUT/"latest_plan.json",dict(path=str(base.relative_to(ROOT))))
    print(base,flush=True)
    run(base,record,args.threads)


if __name__=="__main__":
    main()
