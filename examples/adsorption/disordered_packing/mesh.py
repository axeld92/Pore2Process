"""Check the stock OpenLB export and prepare conservative rectangular face flows."""
import argparse
from pathlib import Path

import numpy as np

from examples.adsorption.dispersion_holdout.run import execute
from pore2process_demo.io import read_json, sha, write_json

from .native import HERE, ORDERED
from .packing import DIAMETER_M, PERIODS

COMPACT = dict(before="<i4", after="<i4", incoming="<f8", outgoing="<f8")
FLOW_FILES = dict(zip(("left.i32","right.i32","axis.u8","unbalanced.f64"),
                      ("<i4","<i4","u1","<f8"),strict=True))


def fluid_mask(centres, shape, h):
    solid = np.zeros(shape,dtype=bool)
    r = DIAMETER_M/2
    for centre in centres:
        lower = np.floor((centre-r)/h-.5).astype(int)
        upper = np.ceil((centre+r)/h-.5).astype(int)
        axes = [np.arange(a,b+1) for a,b in zip(lower,upper,strict=True)]
        grid = np.meshgrid(*axes,indexing="ij",sparse=True)
        test = sum(((index+.5)*h-c)**2 for index,c in zip(grid,centre,strict=True)) <= r**2
        coordinates = np.ix_(*(index%n for index,n in zip(axes,shape,strict=True)))
        solid[coordinates] |= test
    return ~solid


def load_field(base,label,n):
    plan = read_json(base/"plan.json")
    folder = base/f"{label}_N{n}"
    raw = read_json(folder/"raw.json")
    if (raw["solver"],raw["collision"],raw["resolution"],raw["converged"]) != ("OpenLB 1.9r0","TRT",n,True):
        raise ValueError("Expected a converged stock OpenLB TRT field")
    if raw["max_Mach"] > .04 or raw["max_density_departure"] > .002:
        raise ValueError("Native flow exceeds the Mach or density limits")
    data = np.fromfile(folder/"field.f64",dtype="<f8").reshape(raw["fluid_nodes"],7)
    ijk = data[:,:3].astype(np.int32)
    np.testing.assert_array_equal(ijk,data[:,:3])
    h = raw["dx_m"]
    shape = (raw["periods"]*n,n,n)
    centres = np.loadtxt(base/(label+"_centres_m.txt"))
    np.testing.assert_array_equal(ijk,np.argwhere(fluid_mask(centres,shape,h)))
    velocity = data[:,3:6]*plan["velocity_scale"]
    if not np.isfinite(data).all():
        raise ValueError("Non-finite OpenLB field")
    u = velocity.mean(axis=0)
    relative = np.linalg.norm(u-[plan["reference_interstitial_m_s"],0.,0.])/plan["reference_interstitial_m_s"]
    if relative > .002:
        raise ValueError("Native interstitial velocity outside target")
    porosity = len(ijk)/np.prod(shape)
    if abs(porosity-.6) > .001:
        raise ValueError("Voxel carbon inventory departs from the continuous solid fraction")
    directional = velocity.sum(axis=0)*h**2/np.array(shape)
    return folder, ijk, velocity, dict(shape=list(shape), raw=raw, porosity=float(porosity),
                directional_flow_m3_s=directional.tolist(), interstitial_velocity_m_s=u.tolist(),
                velocity_match_relative=float(relative), field_sha256=sha(folder/"field.f64"),
                voxel_solid_volume_m3=float((np.prod(shape)-len(ijk))*h**3))


def shared_faces(ijk,velocity,shape,h):
    ids = np.zeros(shape,dtype=np.int32)
    ids[tuple(ijk.T)] = np.arange(1,len(ijk)+1,dtype=np.int32)
    parts = []
    for axis in range(3):
        neighbor = np.roll(ids,-1,axis=axis)[tuple(ijk.T)]
        owner = np.flatnonzero(neighbor)
        downstream = neighbor[owner]-1
        q = .5*h**2*(velocity[owner,axis]+velocity[downstream,axis])
        parts.append((owner+1,downstream+1,np.full(len(owner),axis+1),q))
    return {name:np.concatenate(col).astype(dtype) for (name,dtype),col in
            zip(FLOW_FILES.items(),zip(*parts,strict=True),strict=True)}


def flow_checks(ijk,shape,info,faces,q):
    a,b,axis = (faces[name]-1 for name in FLOW_FILES if name!="unbalanced.f64")
    ref = np.array(info["directional_flow_m3_s"])
    raw = faces["unbalanced.f64"]
    if q.shape!=raw.shape or not np.isfinite(q).all():
        raise ValueError("Invalid balanced shared faces")
    net = np.bincount(a,weights=q,minlength=len(ijk))-np.bincount(b,weights=q,minlength=len(ijk))
    max_net = np.max(abs(net))/np.max(abs(q))
    max_plane = 0.
    for d,n in enumerate(shape):
        use = axis==d
        rates = np.bincount(ijk[a[use],d],weights=q[use],minlength=n)
        max_plane = max(max_plane,np.max(abs(rates-ref[d]))/ref[0])
    rms_change = np.linalg.norm(q-raw)/np.linalg.norm(raw)
    if max_net > 1e-10 or max_plane > 1e-9 or rms_change > .05:
        raise ValueError(f"Carrier face handoff outside limits: {max_net:g}, {max_plane:g}, {rms_change:g}")
    return dict(max_cell_imbalance_relative=float(max_net), max_plane_error_relative=float(max_plane),
                face_rms_correction_relative=float(rms_change))


def project(base,label,n):
    folder,ijk,velocity,info = load_field(base,label,n)
    d = folder/"faces"
    d.mkdir(exist_ok=True)
    shape = info["shape"]
    faces = shared_faces(ijk,velocity,shape,info["raw"]["dx_m"])
    for name, value in faces.items():
        value.tofile(d/name)
    ijk.astype("<i4").tofile(d/"coords.i32")
    runtime = read_json(ORDERED/"frozen_prediction.json")["parity_command"][:4]
    total = np.array(shape)*info["directional_flow_m3_s"]
    cmd = [*runtime,str(HERE/"balance.jl"),str(d),str(len(ijk)),*map(str,total)]
    print("Face projection start",folder,flush=True)
    elapsed = execute(cmd,d/"balance.log",timeout=7200)
    q = np.fromfile(d/"balanced.f64",dtype="<f8")
    checks = flow_checks(ijk,shape,info,faces,q)
    write_json(d/"summary.json",dict(info=info,checks=checks,elapsed=elapsed,command=cmd,
                                    files={f:sha(d/f) for f in (*FLOW_FILES,"coords.i32","balanced.f64")}))
    print("Face projection complete",folder,checks,flush=True)


def load_faces(folder):
    d = folder/"faces"
    metadata = read_json(d/"summary.json")
    for file,digest in metadata["files"].items():
        if sha(d/file)!=digest:
            raise ValueError("Frozen shared face export changed")
    ijk = np.fromfile(d/"coords.i32",dtype="<i4").reshape(-1,3)
    faces = {name:np.fromfile(d/name,dtype=dtype) for name,dtype in FLOW_FILES.items()}
    return ijk, faces, np.fromfile(d/"balanced.f64",dtype="<f8"), metadata


def forward_cuts(ijk,faces,q,nx):
    left = faces["left.i32"]-1
    axial = faces["axis.u8"]==1
    plane = ijk[left[axial],0]
    minima = np.full(nx,np.inf)
    np.minimum.at(minima,plane,q[axial])
    return minima


def choose_cut(base):
    """Pick the first common physical forward-flowing cut on a 1/16-period grid."""
    cuts = np.arange(8*16)
    for n in (48,64):
        ijk,faces,q,m = load_faces(base/f"disordered_N{n}")
        nx = m["info"]["shape"][0]
        outgoing = (cuts*(n//16)-1)%nx
        minimum = forward_cuts(ijk,faces,q,nx)[outgoing]
        cuts = cuts[np.isfinite(minimum)&(minimum>0)]
    if len(cuts)==0:
        raise ValueError("No common forward-flowing cut in both disordered meshes")
    return int(cuts[0])


def compact(folder,cut,*,write=True):
    ijk,faces,q,m = load_faces(folder)
    shape = m["info"]["shape"]
    nx,n,_ = shape
    left,right,axis = (faces[f]-1 for f in FLOW_FILES if f!="unbalanced.f64")
    offsets = (ijk[:,0]-cut)%nx
    order = np.argsort(offsets,kind="stable")
    index = np.empty(len(ijk),dtype=np.int32)
    index[order] = np.arange(len(ijk),dtype=np.int32)
    exterior = (axis==0)&(offsets[left]==nx-1)
    np.testing.assert_array_equal(offsets[right[exterior]],0)
    if not np.any(exterior) or np.any(q[exterior]<=0):
        raise ValueError("The open tracer cut must be forward-flowing on all fluid faces")
    before,after = (np.zeros((len(ijk),3),dtype="<i4") for _ in range(2))
    low,high = (np.zeros((len(ijk),3),dtype="<f8") for _ in range(2))
    low[index[right],axis] = q
    high[index[left],axis] = q
    inside = ~exterior
    before[index[right[inside]],axis[inside]] = index[left[inside]]+1
    after[index[left[inside]],axis[inside]] = index[right[inside]]+1
    counts = np.searchsorted(offsets[order],np.arange(n,nx+1,n))
    if nx == n:
        outlets,inlets = (index[ids[exterior]] for ids in (left,right))
        count = len(ijk)
        before,after,low,high = (np.tile(a,(PERIODS,1)) for a in (before,after,low,high))
        for period in range(PERIODS):
            start,end = period*count,(period+1)*count
            for topology in (before,after):
                topology[start:end] += start*(topology[start:end]>0)
        for period in range(PERIODS-1):
            start,end = period*count,(period+1)*count
            before[end+inlets,0] = start+outlets+1
            after[start+outlets,0] = end+inlets+1
        counts = count*np.arange(1,PERIODS+1)
    for axis in range(3):
        linked = np.flatnonzero(after[:,axis])
        downstream = after[linked,axis]-1
        np.testing.assert_array_equal(before[downstream,axis],linked+1)
        np.testing.assert_array_equal(low[downstream,axis],high[linked,axis])
    np.testing.assert_allclose(q[exterior].sum(),m["info"]["directional_flow_m3_s"][0],rtol=1e-9,atol=0)
    d = folder/"compact"
    if write:
        d.mkdir(exist_ok=True)
    for (name,dtype),value in zip(COMPACT.items(),(before,after,low,high),strict=True):
        file = d/(name+(".i32" if dtype=="<i4" else ".f64"))
        if write:
            value.tofile(file)
        else:
            np.testing.assert_array_equal(np.fromfile(file,dtype=dtype),value.ravel())
    record = dict(cut_index=int(cut), cut_offset_m=cut*m["info"]["raw"]["dx_m"],
                  period_end_counts=counts.tolist(), inlet_faces=int(exterior.sum()),
                  minimum_inlet_flow_m3_s=float(min(q[exterior])), throughput_m3_s=float(q[exterior].sum()),
                  resolution=n, face_summary_sha256=sha(folder/"faces/summary.json"))
    if write:
        write_json(d/"summary.json",record)
    return record


if __name__=="__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("base",type=Path)
    p.add_argument("label",choices=("ordered","disordered"))
    p.add_argument("grid",type=int)
    a = p.parse_args()
    project(a.base.resolve(),a.label,a.grid)
