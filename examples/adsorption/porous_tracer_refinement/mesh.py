"""Subdivide the saved pore voxels while preserving their face fluxes."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from contextlib import ExitStack
from pathlib import Path

import numpy as np

from examples.adsorption.porous_face_handoff.run import archived
from examples.adsorption.stationary_curved_wall import reference as ref
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT

from .time_half import CENTER, run_solver, records

N = 96
RATIO = 2
FINE_N = N*RATIO
FILES = {'left.i32':'<i4','right.i32':'<i4','axis.u8':'u1','balanced.f64':'<f8'}


def coarse_faces(source):
    return {name:np.fromfile(source/name,dtype=dtype) for name,dtype in FILES.items()}


def coarse_flow(ijk, faces):
    left, right = faces['left.i32']-1, faces['right.i32']-1
    axis, flux = faces['axis.u8']-1, faces['balanced.f64']
    incoming = np.zeros((len(ijk),3))
    outgoing = np.zeros_like(incoming)
    incoming[right,axis] = flux
    outgoing[left,axis] = flux
    edge_ids = np.full(outgoing.shape,-1,dtype=np.int32)
    edge_ids[left,axis] = np.arange(len(flux),dtype=np.int32)
    return incoming,outgoing,edge_ids


def timestep_counts(cap):
    return [int(np.ceil((end-start)/cap)) for start,end in zip((0.,*ref.TIMES[:-1]),ref.TIMES)]


def add_outgoing(total, left, right, q, dx):
    diffusion = ref.D/dx**2
    rate = q/dx**3
    total += np.bincount(left,weights=diffusion+np.maximum(rate,0.),minlength=len(total))
    total += np.bincount(right,weights=diffusion+np.maximum(-rate,0.),minlength=len(total))


def fine_coordinates(ijk):
    offsets = np.indices((RATIO,)*3,dtype=np.int32).reshape(3,-1).T
    return (RATIO*ijk[:,None,:]+offsets).reshape(-1,3)


def refined_mesh(directory, source, info, ijk, *, write):
    faces = coarse_faces(source)
    q0 = faces['balanced.f64']
    incoming,outgoing,edge_ids = coarse_flow(ijk,faces)
    fine = fine_coordinates(ijk)
    count = len(fine)
    ids = np.zeros((FINE_N,)*3,dtype=np.int32)
    ids[tuple(fine.T)] = np.arange(1,count+1,dtype=np.int32)
    dx = info['dx_m']/RATIO
    target = np.asarray(info['carrier_directional_flow_m3_s'])
    max_flux = 0.
    net = np.zeros(count)
    rate_sum = np.zeros(count)
    planes = []
    initial = (.5+ref.A*np.sin(2*np.pi*(fine[:,0]+.5)/FINE_N)).astype('<f8')
    if write:
        directory.mkdir(parents=True)
        initial.tofile(directory/'varying_initial.f64')
    else:
        np.testing.assert_array_equal(initial,np.fromfile(directory/'varying_initial.f64',dtype='<f8'))
    pos = 0
    with ExitStack() as stack:
        streams = {name:stack.enter_context((directory/name).open('wb' if write else 'rb'))
                   for name in FILES}
        for axis in range(3):
            dst = np.roll(ids,-1,axis=axis)[tuple(fine.T)]
            left = np.flatnonzero(dst).astype(np.int32)
            right = dst[left]-1
            parent = left//RATIO**3
            outer = fine[left,axis]%RATIO == RATIO-1
            # Constant normal velocity on each face of a voxel, linear between them.
            q = (incoming[parent,axis]+outgoing[parent,axis])/(2*RATIO**2)
            q[outer] = outgoing[parent[outer],axis]/RATIO**2
            parent_edges = edge_ids[parent[outer],axis]
            if np.any(parent_edges<0):
                raise ValueError('Refinement placed a face through a solid wall')
            coarse_axis = faces['axis.u8']==axis+1
            np.testing.assert_array_equal(np.bincount(parent_edges,minlength=len(q0))[coarse_axis],RATIO**2)
            recovered = np.bincount(parent_edges,weights=q[outer],minlength=len(q0))
            np.testing.assert_array_equal(recovered[coarse_axis],q0[coarse_axis])
            values = (left+1,right+1,np.full(len(left),axis+1,dtype=np.uint8),q)
            for (name,dtype),array in zip(FILES.items(),values,strict=True):
                values_on_disk = array.astype(dtype,copy=False)
                if write:
                    values_on_disk.tofile(streams[name])
                else:
                    np.testing.assert_array_equal(np.fromfile(streams[name],dtype=dtype,count=len(array)),values_on_disk)
            net += np.bincount(left,weights=q,minlength=count)-np.bincount(right,weights=q,minlength=count)
            add_outgoing(rate_sum,left,right,q,dx)
            planes.append(np.bincount(fine[left,axis],weights=q,minlength=FINE_N))
            max_flux = max(max_flux,float(np.max(abs(q))))
            pos += len(left)
        if not write and any(f.read(1) for f in streams.values()):
            raise ValueError('Trailing data in the refined face arrays')
    residual = float(np.max(abs(net))/max_flux)
    plane_error = float(np.max(abs(np.asarray(planes)-target[:,None]))/abs(target[0]))
    volume = count*dx**3
    np.testing.assert_allclose(volume,len(ijk)*info['dx_m']**3,rtol=1e-14,atol=0)
    if residual > 1e-10 or plane_error > 1e-9:
        raise ValueError('Refined mesh fails the flow preservation checks')
    coarse_rates = np.zeros(len(ijk))
    add_outgoing(coarse_rates,faces['left.i32']-1,faces['right.i32']-1,q0,info['dx_m'])
    cap_coarse = .8/float(np.max(coarse_rates))
    cap_fine = .8/float(np.max(rate_sum))
    # A rounding margin keeps the shared physical cap within both stability bounds.
    cap_common = min(.5*cap_coarse,cap_fine)*(1-1e-12)
    matched_case = ('N96_half_step' if timestep_counts(.5*cap_coarse)==timestep_counts(cap_common)
                    else 'N96_matched_step')
    return dict(fine_cells=count,fine_faces=pos,dx_m=dx,liquid_volume_m3=volume,
                maximum_local_imbalance_relative=residual,max_plane_error_relative_to_axial=plane_error,
                coarse_face_flux_preservation='exact',coarse_stability_cap_s=cap_coarse,
                fine_stability_cap_s=cap_fine,common_time_cap_s=cap_common,
                matched_case=matched_case,common_step_counts=timestep_counts(cap_common),
                fine_time_factor=cap_common/cap_fine,matched_coarse_time_factor=cap_common/cap_coarse)


def execute(base):
    base = base.resolve()
    center,_ = records()
    half = base/'N96_half_step'
    if read_json(half/'execution.json')['returncode']:
        raise ValueError('Half-step run failed')
    source = CENTER.parent/'N96'
    info,ijk,_ = archived(N)
    directory = base/'N192_faces'
    audit = refined_mesh(directory,source,info,ijk,write=True)
    write_json(base/'mesh_plan.json',dict(mesh=audit,mesh_source_sha256=sha(Path(__file__))))
    print('Refined face audit',audit,flush=True)
    args = [(base/'N192',directory,audit['dx_m'],center['julia_runtime'],audit['fine_time_factor'])]
    if audit['matched_case']=='N96_matched_step':
        args.append((base/audit['matched_case'],source,info['dx_m'],center['julia_runtime'],
                     audit['matched_coarse_time_factor']))
    with ThreadPoolExecutor(max_workers=2) as pool:
        for future in [pool.submit(run_solver,*case) for case in args]:
            future.result()


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('study',type=Path)
    execute(parser.parse_args().study)
