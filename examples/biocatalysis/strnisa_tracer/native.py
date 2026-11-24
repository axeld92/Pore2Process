"""Whole-length rectangular monolayer stock OpenLB TRT flow study."""
import argparse
import math
import os
from pathlib import Path
import shutil
import subprocess
import tomllib
from datetime import datetime, timezone
from uuid import uuid4

import numpy as np

from examples.adsorption.transport_root_cause.hydraulics.stock_omp import build_case_omp
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT

HERE = Path(__file__).resolve().parent
BUILD = ROOT / ".local/biocatalysis/strnisa_tracer"
OUTPUT = ROOT / "results/biocatalysis/strnisa_tracer"
SOURCE = ROOT / "literature/sources/strnisa2018/evidence.toml"
PROTOCOL = HERE / "protocol.toml"


def config():
    with PROTOCOL.open("rb") as f:
        return tomllib.load(f)


def source():
    with SOURCE.open("rb") as f:
        return tomllib.load(f)


def dimensions(reactor):
    row = next(r for r in source()["table1"]["rows"] if r[0] == reactor)
    return np.array(row[1:4], dtype=float) * 1e-3, float(row[4])


def bead_centers(extent, eps, diameter, seed, sweeps, seam_fraction, clearance_fraction):
    """Randomize a feasible staggered starting geometry without overlaps.

    Beads are confined to the measured rectangular channel, tangent to
    its lower plate. The axial hydraulic periodic seam is bead-free.
    All distances are nondimensionalized by the bead diameter here.
    """
    lx, ly, lz = np.asarray(extent) / diameter
    if lz < 1:
        raise ValueError("Representative beads do not fit in channel")
    n = round((1 - eps) * lx * ly * lz / (math.pi / 6))
    rng = np.random.default_rng(seed)
    # Center limits, including a prescribed clear slot at both ends.
    lower = np.array([.5 + seam_fraction, .5])
    upper = np.array([lx - .5 - seam_fraction, ly - .5])
    clearance = 1 + clearance_fraction
    points = None
    for rows in range(math.floor((ly-1)/(math.sqrt(3)/2*clearance))+1, 2, -1):
        cols = math.ceil(n/rows)
        step = (upper[0]-lower[0])/cols
        dy = (upper[1]-lower[1])/(rows-1)
        if min(step, math.hypot(step/2,dy), 2*dy) < clearance:
            continue
        points = np.array([
            [lower[0]+(i+(.25 if j%2==0 else .75))*step, lower[1]+j*dy]
            for j in range(rows) for i in range(cols)
        ])
        points = points[rng.choice(len(points), n, replace=False)]
        break
    if points is None:
        raise ValueError("Cannot initialize nonoverlapping monolayer")
    # Geometry-generation cell bins only; actual spherical indicators,
    # rasterization and every lattice operation come from stock OpenLB.
    bins = {}
    def bin_key(x):
        return tuple(np.floor(x / clearance).astype(int))
    for i,p in enumerate(points):
        bins.setdefault(bin_key(p),set()).add(i)
    accepted = 0
    for _ in range(sweeps):
        for i in rng.permutation(n):
            candidate = points[i] + rng.uniform(-.24,.24,2)
            if np.any(candidate < lower) or np.any(candidate > upper):
                continue
            key = bin_key(candidate)
            old = bin_key(points[i])
            conflict = False
            for a in range(-1,2):
                for b in range(-1,2):
                    for k in bins.get((key[0]+a,key[1]+b),()):
                        if k != i:
                            delta = points[k] - candidate
                            if delta @ delta < clearance**2:
                                conflict = True
                                break
                    if conflict:
                        break
                if conflict:
                    break
            if not conflict:
                if key != old:
                    bins[old].remove(i)
                    bins.setdefault(key,set()).add(i)
                points[i] = candidate
                accepted += 1
    minimum = min(
        np.min(np.sum((points[i+1:]-points[i])**2,axis=1))
        for i in range(n-1)
    )**.5
    if minimum < clearance - 1e-12:
        raise RuntimeError("Beads overlap")
    beads = np.column_stack((points,np.full(n,.5),np.full(n,.5)))*diameter
    return beads, dict(
        count=n, minimum_gap_m=(minimum-1)*diameter,
        continuous_porosity=1-n*math.pi*diameter**3/6/np.prod(extent),
        acceptance=accepted/(n*sweeps), seam_min_surface_clearance_m=seam_fraction*diameter)


def voxel_mask(beads, shape, h):
    """Independent verification of stock indicator geometry at fluid centers."""
    nx,ny,nz = shape
    mask = np.zeros((nx,ny+2,nz+2), dtype=bool)
    mask[:,1:-1,1:-1] = True
    offset = np.array([.5,-.5,-.5])
    padded_shape = np.array(mask.shape)
    for x,y,z,r in beads:
        center = np.array([x,y,z])
        lo = np.maximum(0,np.floor((center-r)/h-offset).astype(int))
        hi = np.minimum(padded_shape-1,np.ceil((center+r)/h-offset).astype(int))
        inds = [np.arange(a,b+1) for a,b in zip(lo,hi)]
        grid = np.meshgrid(*inds, indexing="ij", sparse=True)
        distance2 = sum(((q+off)*h-c)**2 for q,off,c in zip(grid,offset,center))
        mask[np.ix_(*inds)] &= distance2 > r*r
    return mask


def prepare():
    p = config()
    g = p["geometry"]
    base = OUTPUT / "runs" / (datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ") + "-" + uuid4().hex[:8])
    base.mkdir(parents=True)
    reactors = {}
    for name in ("A","C"):
        extent,eps = dimensions(name)
        beads, record = bead_centers(extent,eps,g["diameter_um"]*1e-6,g[f"seed_{name}"],
                          g["sweeps"],g["seam_buffer_fraction_of_diameter"],g["minimum_bead_gap_fraction"])
        f = base / f"{name}_beads_m.txt"
        np.savetxt(f,beads,fmt="%.17g")
        reactors[name] = dict(extent_m=extent.tolist(),observed_porosity=eps,
                              geometry=record,centres_sha256=sha(f))
    write_json(base/"plan.json",dict(reactors=reactors,protocol=p,
                    protocol_sha256=sha(PROTOCOL),source_sha256=sha(SOURCE),
                    source_doi=p["source_doi"]))
    write_json(OUTPUT/"latest_plan.json",dict(path=str(base.relative_to(ROOT))))
    print(base,flush=True)
    return base


def build():
    return build_case_omp(HERE/"openlb",BUILD/"build",BUILD/"record")


def native(base,reactor,nz,tau=None,warm=None):
    plan = read_json(base/"plan.json")
    p=plan["protocol"]
    r=plan["reactors"][reactor]
    if sha(PROTOCOL)!=plan["protocol_sha256"] or sha(SOURCE)!=plan["source_sha256"]:
        raise ValueError("Frozen protocol or source changed")
    file=base/f"{reactor}_beads_m.txt"
    if sha(file)!=r["centres_sha256"]:
        raise ValueError("Frozen packing changed")
    extent=np.asarray(r["extent_m"])
    h=extent[2]/nz
    shape=tuple(int(v) for v in np.rint(extent/h))
    np.testing.assert_allclose(np.array(shape)*h,extent,rtol=1e-12,atol=0)
    f=p["fluid"]
    n=p["numerics"]
    tau=n["tau_even"] if tau is None else tau
    q=f["flow_ul_min"]*1e-9/60
    record=read_json(BUILD/"record/build.json")
    for name,digest in record["source_sha256"].items():
        if sha(HERE/"openlb"/name)!=digest:
            raise ValueError("Application changed after build")
    if sha(Path(record["binary"]))!=record["binary_sha256"]:
        raise ValueError("Binary changed")
    directory=base/f"{reactor}_H{nz}"
    directory.mkdir()
    mask=voxel_mask(np.loadtxt(file),shape,h)
    porosity=mask.sum()/np.prod(shape)
    if abs(porosity-r["observed_porosity"])>.01:
        raise ValueError(f"Voxel porosity misses source by {porosity-r['observed_porosity']}")
    command=[record["binary"],*map(str,(*shape,h,tau,q,f["viscosity_Pa_s"],f["density_kg_m3"],file))]
    warm_hash=None
    if warm is not None:
        warm=warm.resolve()
        old=read_json(warm/"raw.json")
        old_request=read_json(warm/"request.json")
        if (old_request["centres_sha256"]!=sha(file) or
                old["shape_with_walls"]!=[shape[0],shape[1]+2,shape[2]+2] or
                old["fluid_nodes"]!=int(mask.sum()) or
                old["dx_m"]!=h or old["collision"]!="TRT" or not old["converged"]):
            raise ValueError("Warm start must be a converged field on this exact geometry")
        warm_hash=sha(warm/"field.f64")
        command += [str(warm/"field.f64"),str(old["tau_even"]),str(old["gradient_Pa_m"])]
    shutil.copyfile(BUILD/"record/build.json",directory/"build.json")
    write_json(directory/"request.json",dict(command=command,threads=n["threads"],
                   physical_shape=list(shape),voxel_porosity=porosity,
                   plan_sha256=sha(base/"plan.json"),centres_sha256=sha(file),
                   warm_field_sha256=warm_hash))
    print("Native stock TRT start:",directory,"cells including walls",list(mask.shape),flush=True)
    env={**os.environ,"OMP_NUM_THREADS":str(n["threads"]),"OMP_DYNAMIC":"FALSE"}
    with (directory/"native.log").open("w") as log:
        subprocess.run(command,cwd=directory,env=env,stdout=log,stderr=subprocess.STDOUT,
                       check=True,timeout=14400)
    raw=read_json(directory/"raw.json")
    if (raw["solver"],raw["collision"],raw["descriptor"],raw["converged"])!=("OpenLB 1.9r0","TRT","D3Q19",True):
        raise ValueError("Native stock TRT provenance failed")
    if raw["max_Mach"]>n["max_Mach"] or raw["max_density_departure"]>n["max_density_departure"]:
        raise ValueError("Native low-Mach/density departure limit failed")
    if abs(raw["flow_m3_s"]/q-1)>n["max_native_flow_relative_error"]:
        raise ValueError("Native carrier throughput mismatch")
    data=np.memmap(directory/"field.f64",dtype="<f8",mode="r",shape=(raw["fluid_nodes"],7))
    if not np.isfinite(data).all():
        raise ValueError("Nonfinite native field")
    ijk=data[:,:3].astype(np.int32)
    np.testing.assert_array_equal(ijk,data[:,:3])
    np.testing.assert_array_equal(ijk,np.argwhere(mask))
    write_json(directory/"hydraulic.json",dict(native=raw,shape_with_walls=list(mask.shape),
                  target_flow_m3_s=q,voxel_porosity=porosity,field_sha256=sha(directory/"field.f64")))
    print("Native TRT complete:",directory,raw,flush=True)
    return directory


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    s=parser.add_subparsers(dest="command",required=True)
    s.add_parser("build")
    s.add_parser("prepare")
    run=s.add_parser("run")
    run.add_argument("base",type=Path)
    run.add_argument("reactor",choices=("A","C"))
    run.add_argument("depth_cells",type=int,choices=(8,12,16,20,24,28,32))
    run.add_argument("--tau",type=float)
    run.add_argument("--warm",type=Path)
    a=parser.parse_args()
    if a.command=="prepare":
        prepare()
    elif a.command=="build":
        print(build()["binary"])
    else:
        native(a.base.resolve(),a.reactor,a.depth_cells,a.tau,a.warm)


if __name__=="__main__":
    main()
