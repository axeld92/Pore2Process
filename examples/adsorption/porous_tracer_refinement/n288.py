"""Refine the fixed archived N96 flow to N288 and match the N192 time grid."""
import argparse
from contextlib import ExitStack
from datetime import datetime, timezone
from pathlib import Path
import shutil
import subprocess
from time import perf_counter
from uuid import uuid4

import numpy as np

from examples.adsorption.porous_face_handoff.run import SNAPSHOTS, archived, scalar_metrics
from examples.adsorption.stationary_curved_wall import reference as ref
from examples.adsorption.transport_accuracy.run import validate_record
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT
from pore2process_demo.provenance import source_manifest

from .analyze import variance_balance
from .mesh import FILES, N, add_outgoing, coarse_faces, coarse_flow
from .n192_time import PARENT, RESULTS as N192_RESULTS, inputs
from .time_half import CENTER, OUTPUT

HERE = Path(__file__).parent
RESULTS = OUTPUT/'N288_space'
CONTROL = N192_RESULTS/'20261111T162208Z-42f74625/summary.json'
CONTROL_SHA = '9733d80c58f58e903ae6ba692a9418f0b48494545f33f7b0a5d68c966a0861cf'
RATIO = 3
FINE_N = RATIO*N
THREADS = 8


def records():
    if sha(CONTROL)!=CONTROL_SHA:
        raise ValueError('Archived N192 half-step summary changed')
    old = validate_record(CONTROL)
    parent,runtime,plan = inputs()
    if old['plan']!=plan or not old['analysis']['checks']['exact_half_steps']:
        raise ValueError('N192 control time grid differs from the archived plan')
    return old,parent,runtime,plan


def child_coordinates(ijk):
    offsets = np.indices((RATIO,)*3,dtype=np.int32).reshape(3,-1).T
    return (RATIO*ijk[:,None,:]+offsets).reshape(-1,3)


def refine(directory,info,coarse,fine,*,write):
    faces = coarse_faces(CENTER.parent/'N96')
    q0 = faces['balanced.f64']
    incoming,outgoing,edge_ids = coarse_flow(coarse,faces)
    count = len(fine)
    ids = np.zeros((FINE_N,)*3,dtype=np.int32)
    ids[tuple(fine.T)] = np.arange(1,count+1,dtype=np.int32)
    dx = info['dx_m']/RATIO
    target = np.asarray(info['carrier_directional_flow_m3_s'])
    local = np.zeros(count)
    rates = np.zeros(count)
    planes = []
    total_edges = 0
    maximum_flux = 0.
    parent_error = 0.
    initial = (.5+ref.A*np.sin(2*np.pi*(fine[:,0]+.5)/FINE_N)).astype('<f8')
    if write:
        directory.mkdir(parents=True)
        initial.tofile(directory/'varying_initial.f64')
    else:
        np.testing.assert_array_equal(initial,np.fromfile(directory/'varying_initial.f64',dtype='<f8'))
    with ExitStack() as stack:
        streams = {name:stack.enter_context((directory/name).open('wb' if write else 'rb'))
                   for name in FILES}
        for axis in range(3):
            destination = np.roll(ids,-1,axis=axis)[tuple(fine.T)]
            left = np.flatnonzero(destination).astype(np.int32)
            right = destination[left]-1
            parent = left//RATIO**3
            offset = fine[left,axis]%RATIO+1
            s = offset/RATIO
            # The positive normal component is linear from the negative face
            # to the positive face of each original voxel.
            q = ((1-s)*incoming[parent,axis]+s*outgoing[parent,axis])/RATIO**2
            outer = offset==RATIO
            parent_edges = edge_ids[parent[outer],axis]
            if np.any(parent_edges<0):
                raise ValueError('Unexpected subface through a wall')
            coarse_axis = faces['axis.u8']==axis+1
            np.testing.assert_array_equal(np.bincount(parent_edges,minlength=len(q0))[coarse_axis],RATIO**2)
            recovered = np.bincount(parent_edges,weights=q[outer],minlength=len(q0))[coarse_axis]
            np.testing.assert_allclose(recovered,q0[coarse_axis],rtol=3e-15,atol=0)
            parent_error = max(parent_error,float(np.max(abs(recovered-q0[coarse_axis]))))
            values = (left+1,right+1,np.full(len(left),axis+1,dtype=np.uint8),q)
            for (name,dtype),data in zip(FILES.items(),values,strict=True):
                data = data.astype(dtype,copy=False)
                if write:
                    data.tofile(streams[name])
                else:
                    np.testing.assert_array_equal(np.fromfile(streams[name],dtype=dtype,count=len(data)),data)
            local += np.bincount(left,weights=q,minlength=count)-np.bincount(right,weights=q,minlength=count)
            add_outgoing(rates,left,right,q,dx)
            planes.append(np.bincount(fine[left,axis],weights=q,minlength=FINE_N))
            maximum_flux = max(maximum_flux,float(np.max(abs(q))))
            total_edges += len(q)
        if not write and any(f.read(1) for f in streams.values()):
            raise ValueError('Trailing data in refined faces')
    volume = count*dx**3
    np.testing.assert_allclose(volume,len(coarse)*info['dx_m']**3,rtol=1e-14,atol=0)
    imbalance = float(np.max(abs(local))/maximum_flux)
    plane_error = float(np.max(abs(np.asarray(planes)-target[:,None]))/abs(target[0]))
    if imbalance>1e-10 or plane_error>1e-9:
        raise ValueError('N288 face reconstruction failed balance or throughput')
    return dict(resolution=FINE_N,cells=count,faces=total_edges,dx_m=dx,liquid_volume_m3=volume,
                archived_axial_flow_m3_s=info['carrier_directional_flow_m3_s'][0],
                max_plane_error_relative_to_axial=plane_error,maximum_local_imbalance_relative=imbalance,
                maximum_parent_face_error_relative_to_largest_parent_flow=parent_error/float(np.max(abs(q0))),
                stable_time_cap_s=.8/float(np.max(rates)))


def command(runtime,mode,source,directory,dx,steps):
    args = [mode,source,directory,dx,ref.D,*ref.TIMES,*steps,THREADS]
    return [runtime[0],f'--threads={THREADS}',*runtime[1:],str(HERE/'solve_n288.jl'),*map(str,args)]


def execute(cmd,log_path,*,timeout):
    start = perf_counter()
    with log_path.open('w') as log:
        process = subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,timeout=timeout)
    wall_s = perf_counter()-start
    if process.returncode:
        raise RuntimeError(f'Julia process failed. See {log_path}')
    return dict(returncode=process.returncode,wall_s=wall_s)


def kernel_comparison(path):
    lines = path.read_text().splitlines()
    if len(lines)!=4 or lines[0]!='field,max_concentration_difference_kg_m3,relative_dissipation_difference':
        raise ValueError('Unexpected serial-to-threaded verification output')
    audit = {}
    for line in lines[1:]:
        name,difference,dissipation = line.split(',')
        delta,relative = float(difference),float(dissipation)
        if delta!=0 or not 0<=relative<=1e-10:
            raise ValueError('Threaded computation differs from archived serial scheme')
        audit[name] = dict(max_concentration_difference_kg_m3=delta,relative_dissipation_difference=relative)
    if set(audit)!={'initial','middle','final'}:
        raise ValueError('Missing serial-to-threaded verification case')
    return audit


def request(base,runtime,plan,mesh):
    n192 = CONTROL.parent/'N192_half_step'
    steps = plan['new_step_counts']
    return dict(verify=command(runtime,'verify',PARENT.parent/'N192_faces',n192,plan['dx_m'],steps),
                solve=command(runtime,'solve',base/'N288_faces',base/'N288',mesh['dx_m'],steps),
                step_counts=steps,effective_dt_s=plan['new_dt_s'],
                flow_sha256=sha(base/'N288_faces/balanced.f64'),
                initial_sha256=sha(base/'N288_faces/varying_initial.f64'))


def compare(base,old,parent,runtime,plan,mesh,ijk):
    expected = request(base,runtime,plan,mesh)
    if read_json(base/'request.json')!=expected or read_json(base/'execution.json')['returncode']:
        raise ValueError('Saved Julia run differs from the request')
    directory = base/'N288'
    initial = np.fromfile(base/'N288_faces/varying_initial.f64',dtype='<f8')
    metrics,_ = scalar_metrics(directory,'varying',initial,ijk,FINE_N)
    times = np.genfromtxt(directory/'varying_history.csv',delimiter=',',names=True)['time_s']
    old_history = CONTROL.parent/'N192_half_step'/'varying_history.csv'
    np.testing.assert_array_equal(times,np.genfromtxt(old_history,delimiter=',',names=True)['time_s'])
    if metrics['history_steps']!=sum(plan['new_step_counts']):
        raise ValueError('N288 and N192 time schedules differ')
    losses = variance_balance(metrics,float(np.var(initial)))
    change = {}
    for t in SNAPSHOTS:
        baseline = old['analysis']['half_step']['snapshots'][t]
        current = metrics['snapshots'][t]
        v = current['variance_kg2_m6']/baseline['variance_kg2_m6']-1
        a = current['fundamental_amplitude_kg_m3']/baseline['fundamental_amplitude_kg_m3']-1
        phase = float(np.angle(np.exp(1j*(current['phase_rad']-baseline['phase_rad']))))
        change[t] = dict(relative_variance_change=v,relative_fundamental_amplitude_change=a,
                         fundamental_phase_change_rad=phase,
                         variance_passed=abs(v)<=.05,amplitude_passed=abs(a)<=.05,phase_passed=abs(phase)<=.05)
        metrics['checks'].pop(t+'_variance_budget')
    checks = dict(same_time_grid=True,flow_preserved=True,
                  bounds_and_inventory=metrics['checks']['mass'] and metrics['checks']['physical_bounds'],
                  spatial_variance=all(t['variance_passed'] for t in change.values()),
                  spatial_amplitude=all(t['amplitude_passed'] for t in change.values()),
                  spatial_phase=all(t['phase_passed'] for t in change.values()),
                  numerical_share=all(x['within_numerical_limit'] for x in losses.values()))
    return dict(mesh=mesh,serial_equivalence=kernel_comparison(base/'serial_equivalence.csv'),
                N192=old['analysis']['half_step'],N288=metrics,
                N192_variance_budget=old['analysis']['variance_budget']['half_step'],
                N288_variance_budget=losses,changes=change,checks=checks,
                prior_N96_N192=parent['analysis']['spatial_change'])


def run():
    old,parent,runtime,plan = records()
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid4().hex[:8]
    base = RESULTS/stamp
    base.mkdir(parents=True)
    files = [HERE/k for k in ('n288.py','solve_n288.jl','n288_protocol.md','mesh.py','analyze.py','n192_time.py')]
    files += [ROOT/'julia/adsorption'/k for k in ('BoundedFaceTransport.jl','ThreadedBoundedFaceTransport.jl')]
    hashes = source_manifest(files)
    for name in hashes:
        target = base/'source_snapshot'/name
        target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(ROOT/name,target)
    record = dict(schema='pore2process.N288-space-check.v1',run_id=stamp,source_sha256=hashes,
                  control_summary=str(CONTROL.relative_to(ROOT)),control_summary_sha256=CONTROL_SHA,
                  new_openlb_runs=0,new_bed_runs=0,new_DWSIM_runs=0)
    write_json(base/'plan.json',record)
    print('N288 tracer study',base,flush=True)
    info,ijk,_ = archived(N)
    fine = child_coordinates(ijk)
    mesh = refine(base/'N288_faces',info,ijk,fine,write=True)
    if max(plan['new_dt_s'])>mesh['stable_time_cap_s']:
        raise ValueError('Archived N192 time steps exceed the N288 stability cap')
    req = request(base,runtime,plan,mesh)
    write_json(base/'mesh.json',mesh)
    write_json(base/'request.json',req)
    print('Refined flow audit',mesh,flush=True)
    verify = execute(req['verify'],base/'serial_equivalence.csv',timeout=300)
    print('Threaded equivalence',kernel_comparison(base/'serial_equivalence.csv'),flush=True)
    directory = base/'N288'
    directory.mkdir()
    execution = execute(req['solve'],directory/'stdout.log',timeout=7200)
    write_json(base/'execution.json',dict(execution,serial_equivalence_wall_s=verify['wall_s']))
    print('N288 Julia complete',execution,flush=True)
    record['analysis'] = compare(base,old,parent,runtime,plan,mesh,fine)
    record['artifacts_sha256'] = {str(p.relative_to(base)):sha(p) for p in base.rglob('*') if p.is_file()}
    path = base/'summary.json'
    write_json(path,record)
    write_json(RESULTS/'latest.json',dict(summary=str(path.relative_to(ROOT))))
    print('Completed N288 tracer study',path,record['analysis']['checks'],flush=True)


def check(path):
    path = path.resolve()
    record = validate_record(path)
    old,parent,runtime,plan = records()
    if record['control_summary_sha256']!=CONTROL_SHA:
        raise ValueError('Wrong N192 control')
    info,ijk,_ = archived(N)
    fine = child_coordinates(ijk)
    mesh = refine(path.parent/'N288_faces',info,ijk,fine,write=False)
    if mesh!=read_json(path.parent/'mesh.json') or record['analysis']!=compare(
            path.parent,old,parent,runtime,plan,mesh,fine):
        raise ValueError('N288 result differs on replay')
    result = dict(summary=str(path.relative_to(ROOT)),summary_sha256=sha(path),sources=len(record['source_sha256']),
                  artifacts=len(record['artifacts_sha256']),checks=record['analysis']['checks'])
    write_json(RESULTS/'verification.json',result)
    print(result)


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check',type=Path)
    args = parser.parse_args()
    check(args.check) if args.check else run()
