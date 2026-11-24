"""Conservative face handoff at the measured throughput, then inert first passage."""
import argparse
import math
from pathlib import Path
import subprocess
from time import perf_counter

import numpy as np

from examples.adsorption.breakthrough_hypotheses.hydrodynamics.handoff import (
    conservation, connected_pores, julia,
)
from examples.adsorption.disordered_packing.mesh import FLOW_FILES,shared_faces
from pore2process_demo.io import read_json,sha,write_json
from pore2process_demo.paths import ROOT

from .native import HERE,voxel_mask


def handoff(base,reactor,nz,diagnostic=False):
    p=read_json(base/"plan.json")["protocol"]
    n=p["numerics"]
    directory=base/f"{reactor}_H{nz}"
    info=read_json(directory/"hydraulic.json")
    raw=info["native"]
    if sha(directory/"field.f64")!=info["field_sha256"]:
        raise ValueError("Carrier archive changed")
    shape=tuple(info["shape_with_walls"])
    h=raw["dx_m"]
    rows=np.memmap(directory/"field.f64",dtype="<f8",mode="r",shape=(raw["fluid_nodes"],7))
    ijk=rows[:,:3].astype(np.int32)
    velocity=np.array(rows[:,3:6])
    material=voxel_mask(np.loadtxt(base/f"{reactor}_beads_m.txt"),
                        (shape[0],shape[1]-2,shape[2]-2),h)
    np.testing.assert_array_equal(ijk,np.argwhere(material))
    # The existing connectivity audit also uses roll; lateral solid
    # padding ensures that no transverse seam can carry flow.
    active=connected_pores(ijk,material,velocity)
    removed=int(np.count_nonzero(~active))
    qnative=raw["flow_m3_s"]
    active_flow_fraction=float(velocity[active,0].sum()/velocity[:,0].sum())
    ijk=ijk[active]
    velocity=velocity[active]
    faces=shared_faces(ijk,velocity,shape,h)
    d=directory/"faces"
    d.mkdir()
    for name,arr in faces.items():
        arr.tofile(d/name)
    ijk.astype("<i4").tofile(d/"coords.i32")
    command=[*julia(1),str(HERE/"balance.jl"),str(d),str(len(ijk)),str(shape[0]*qnative)]
    print("Projecting shared faces",directory,flush=True)
    start=perf_counter()
    with (d/"balance.log").open("w") as log:
        subprocess.run(command,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=5400)
    q=np.fromfile(d/"balanced.f64",dtype="<f8")
    checks=conservation(ijk,dict(shape=shape,directional_flow_m3_s=[qnative,0.,0.]),faces,q)
    write_json(d/"screen.json",dict(checks=checks,
                field_sha256=info["field_sha256"],
                files={f:sha(d/f) for f in (*FLOW_FILES,"coords.i32","balanced.f64")}))
    if not checks["qualifies_local_correction"] and not diagnostic:
        raise ValueError("The face projection changes the native flow by more than 5% RMS")
    a,b,axis=(faces[f]-1 for f in ("left.i32","right.i32","axis.u8"))
    exterior=(axis==0)&(ijk[a,0]==shape[0]-1)
    if not np.any(exterior) or np.any(ijk[b[exterior],0]!=0):
        raise ValueError("No periodic seam")
    if np.any(q[exterior]<=0):
        raise ValueError("Tracer inlet/outlet seam is not wholly forward-flowing")
    # Every shared face contributes identically to the adjacent voxels.
    # Remove only the seam neighbor IDs, keeping positive open face rates.
    before=np.zeros((len(ijk),3),dtype="<i4")
    after=np.zeros_like(before)
    low=np.zeros((len(ijk),3),dtype="<f8")
    high=np.zeros_like(low)
    low[b,axis]=q
    high[a,axis]=q
    interior=~exterior
    before[b[interior],axis[interior]]=a[interior]+1
    after[a[interior],axis[interior]]=b[interior]+1
    # Check all zero-neighbor fluxes are wall-zero or the axial open faces.
    valid_low=(before!=0) | (low==0)
    valid_high=(after!=0) | (high==0)
    valid_low[b[exterior],0]=True
    valid_high[a[exterior],0]=True
    if not np.all(valid_low) or not np.all(valid_high):
        raise ValueError("Nonzero flow into a solid wall")
    midpoint=int(np.count_nonzero(ijk[:,0]<shape[0]//2))
    d2=directory/"compact"
    d2.mkdir()
    filenames=("before.i32","after.i32","incoming.f64","outgoing.f64")
    for name,arr in zip(filenames,(before,after,low,high)):
        arr.tofile(d2/name)
    record=dict(checks=checks,carrier_qualified=checks["qualifies_local_correction"],
                diagnostic_only=not checks["qualifies_local_correction"],native_throughput_m3_s=qnative,
                pore_volume_m3=len(ijk)*h**3,ideal_mean_s=len(ijk)*h**3/qnative,
                removed_corner_only_cells=removed,retained_native_axial_flow_fraction=active_flow_fraction,
                seam_faces=int(exterior.sum()),seam_min_m3_s=float(q[exterior].min()),
                seam_flow_m3_s=float(q[exterior].sum()),midpoint_cell_count=midpoint,
                dx_m=h,cells=len(ijk),face_projection_s=perf_counter()-start,
                field_sha256=info["field_sha256"],command=command,
                face_files_sha256={f:sha(d/f) for f in (*FLOW_FILES,"coords.i32","balanced.f64")},
                compact_files_sha256={f:sha(d2/f) for f in filenames})
    write_json(directory/"handoff.json",record)
    print("Conservative handoff complete",directory,checks,flush=True)


def matched_dt(protocol,finest_depth_cells,diffusivity_factor=1.2):
    """Match time steps across grids while staying within 12-sigma criterion."""
    h=.0004/finest_depth_cells
    d=protocol["fluid"]["glucose_diffusivity_m2_s"]*diffusivity_factor
    limit=.98*h*h/(288*d)
    cap=protocol["numerics"]["particle_dt_s"]
    return cap/2**max(0,math.ceil(math.log2(cap/limit)))


def tracer(base,reactor,nz,dt,factor=1.,count=None,seed=None,diagnostic=False):
    p=read_json(base/"plan.json")["protocol"]
    n=p["numerics"]
    directory=base/f"{reactor}_H{nz}"
    info=read_json(directory/"handoff.json")
    if not info["carrier_qualified"] and not diagnostic:
        raise ValueError("Carrier failed the projection screen; tracer requires explicit diagnostic mode")
    d=directory/"compact"
    for file,digest in info["compact_files_sha256"].items():
        if sha(d/file)!=digest:
            raise ValueError("Balanced field archive changed")
    count=n["particles"] if count is None else count
    seed=n["particle_seed"] if seed is None else seed
    stop=n[f"marker_stop_{reactor}_s"]
    diffusion=factor*p["fluid"]["glucose_diffusivity_m2_s"]
    tag=f"dt{dt:.8g}_D{factor:.4g}_n{count}_seed{seed}"
    path=directory/(tag+".f64")
    if path.exists():
        raise ValueError("Marker result already exists")
    command=[*julia(n["threads"]),str(HERE/"track.jl"),str(d),str(path),
             *map(str,(info["midpoint_cell_count"],info["dx_m"],diffusion,dt,stop,count,seed))]
    start=perf_counter()
    print("Tracer first passage",reactor,nz,tag,flush=True)
    with (directory/(tag+".log")).open("w") as log:
        subprocess.run(command,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=5400)
    times=np.fromfile(path,dtype="<f8")
    if len(times)!=count or np.any(times[np.isfinite(times)]<=0):
        raise ValueError("Invalid exit time ledger")
    arrived=times[np.isfinite(times)]
    recovery=len(arrived)/count
    mean=float(np.mean(arrived))
    variance=float(np.var(arrived))
    # Conditional moments are not presented as total RTD when recovery
    # fails; unresolved tracer remains explicit.
    block_means=[]; block_variances=[]
    for block in np.array_split(times,8):
        vals=block[np.isfinite(block)]
        block_means.append(float(np.mean(vals)))
        block_variances.append(float(np.var(vals)))
    record=dict(reactor=reactor,depth_cells=nz,dt_s=dt,
                factor=factor,D_m2_s=diffusion,particle_count=count,seed=seed,stop_s=stop,
                recovery=recovery,conditional_mean_s=mean,conditional_variance_s2=variance,
                conditional_skewness=float(np.mean((arrived-mean)**3)/variance**1.5),
                block_mean_sem_s=float(np.std(block_means,ddof=1)/np.sqrt(8)),
                block_variance_sem_s2=float(np.std(block_variances,ddof=1)/np.sqrt(8)),
                mean_vs_volume_relative=mean/info["ideal_mean_s"]-1,
                time_quantiles_s=np.quantile(arrived,[.01,.05,.1,.5,.9,.95,.99]).tolist(),
                numerical_mass_recovery_passed=recovery>=n["min_recovery"],
                carrier_qualified=info["carrier_qualified"],
                marker_sha256=sha(path),handoff_sha256=sha(directory/"handoff.json"),
                command=command,wall_s=perf_counter()-start)
    write_json(directory/(tag+".json"),record)
    print("Tracer result",tag,record,flush=True)
    return record


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("base",type=Path)
    parser.add_argument("reactor",choices=("A","C"))
    parser.add_argument("depth_cells",type=int)
    parser.add_argument("action",choices=("handoff","trace"))
    parser.add_argument("--dt",type=float,default=.0025)
    parser.add_argument("--Dfactor",type=float,default=1.)
    parser.add_argument("--count",type=int)
    parser.add_argument("--seed",type=int)
    parser.add_argument("--diagnostic",action="store_true")
    a=parser.parse_args()
    if a.action=="handoff":
        handoff(a.base.resolve(),a.reactor,a.depth_cells,a.diagnostic)
    else:
        tracer(a.base.resolve(),a.reactor,a.depth_cells,a.dt,a.Dfactor,a.count,a.seed,a.diagnostic)


if __name__=="__main__":
    main()
