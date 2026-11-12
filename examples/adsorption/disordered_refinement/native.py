"""Continue the archived disordered packing on a stock TRT N96 lattice."""
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
from examples.adsorption.disordered_packing.native import BUILD
from examples.adsorption.disordered_packing.native import HERE as PACKING
from examples.adsorption.disordered_packing.packing import PERIODS, SIDE_M
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT

HERE = Path(__file__).resolve().parent
OUTPUT = ROOT/"results/adsorption/disordered_refinement"
PARENT = ROOT/"results/adsorption/disordered_packing/runs/20261112T161615Z-52964a01"
PARENT_SHA = "526ea6083ca472166363589f81ba9046a16e8765bd040382e8f047d7d054aac7"
GRID = 96
TAU = .58


def parent():
    path = PARENT/"summary.json"
    if sha(path) != PARENT_SHA:
        raise ValueError("Pinned disordered packing comparison changed")
    record = read_json(path)
    for name in ("plan.json","disordered_centres_m.txt","disordered_N64/raw.json",
                 "disordered_N64/field.f64","disordered_N64/tracer_summary.json"):
        if sha(PARENT/name) != record["artifacts_sha256"][name]:
            raise ValueError(f"Pinned N64 input changed: {name}")
    return record


def interpolate(values, target_shape):
    """Periodic trilinear interpolation between the two cell-centred meshes."""
    for axis,size in enumerate(target_shape):
        old = values.shape[axis]
        coordinate = (np.arange(size)+.5)*old/size-.5
        lo = np.floor(coordinate).astype(int)
        weight_shape = [1]*values.ndim
        weight_shape[axis] = size
        w = (coordinate-lo).reshape(weight_shape)
        values = ((1.-w)*np.take(values,lo%old,axis=axis)
                  +w*np.take(values,(lo+1)%old,axis=axis))
    return values


def initial_field(base):
    old, ijk, _, info = load_field(PARENT,"disordered",64)
    native = np.fromfile(old/"field.f64",dtype="<f8").reshape(-1,7)
    moments = np.zeros((*info["shape"],4))
    moments[tuple(ijk.T)] = native[:,3:]
    moments[tuple(ijk.T)+(3,)] -= 1000.
    shape = (PERIODS*GRID,GRID,GRID)
    centres = np.loadtxt(base/"disordered_centres_m.txt")
    fine = np.argwhere(fluid_mask(centres,shape,SIDE_M/GRID))
    mapped = interpolate(moments,shape)[tuple(fine.T)]
    mapped[:,3] += 1000.
    target = base/"initial_N96.f64"
    np.column_stack((fine,mapped)).astype("<f8",copy=False).tofile(target)
    return target


def run(base,threads):
    config = read_json(base/"plan.json")
    old = read_json(PARENT/"disordered_N64/raw.json")
    building = read_json(BUILD/"build_record/build.json")
    for name,digest in building["source_sha256"].items():
        if sha(PACKING/"openlb"/name)!=digest:
            raise ValueError("External OpenLB application source changed")
    binary = Path(building["binary"])
    if sha(binary)!=building["binary_sha256"]:
        raise ValueError("Stock TRT executable changed")
    folder = base/f"disordered_N{GRID}"
    folder.mkdir()
    warm = initial_field(base)
    density_scale = (((TAU-.5)/GRID)/((old["tau_even"]-.5)/64))**2
    centres = base/"disordered_centres_m.txt"
    cmd = [str(binary),*map(str,(GRID,PERIODS,SIDE_M,TAU,
        config["reference_interstitial_m_s"]/config["velocity_scale"],*old["gradient_Pa_m"],
        centres,warm,density_scale))]
    shutil.copy2(BUILD/"build_record/build.json",folder/"build.json")
    write_json(folder/"request.json",dict(command=cmd,omp_threads=threads,
        parent_sha256=PARENT_SHA,warm_field_sha256=sha(warm),centres_sha256=sha(centres),
        native_source_sha256=sha(Path(__file__)),protocol_sha256=sha(HERE/"protocol.md")))
    env = {**os.environ,"OMP_NUM_THREADS":str(threads),"OMP_DYNAMIC":"FALSE"}
    print("Stock TRT N96 start",folder,flush=True)
    start = perf_counter()
    with (folder/"native.log").open("w") as log:
        result = subprocess.run(cmd,cwd=folder,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=14400)
    write_json(folder/"execution.json",dict(returncode=result.returncode,wall_s=perf_counter()-start))
    if result.returncode:
        raise RuntimeError(f"Stock TRT N96 failed. See {folder/'native.log'}")
    _,_,_,info = load_field(base,"disordered",GRID)
    write_json(folder/"verified_field.json",info)
    print("Stock TRT N96 accepted",info,flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--threads",type=int,default=10)
    args = parser.parse_args()
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
