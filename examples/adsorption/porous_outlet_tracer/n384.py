"""Compare the archived N288 outlet tracer with N384 at the same step."""
import argparse
from datetime import datetime, timezone
from pathlib import Path
import shutil
import subprocess
from time import perf_counter
from uuid import uuid4

import numpy as np

from examples.adsorption.porous_tracer_refinement.mesh import coarse_faces, coarse_flow
from examples.adsorption.porous_tracer_refinement.time_half import CENTER
from examples.adsorption.stationary_curved_wall import reference as ref
from examples.adsorption.transport_accuracy.run import validate_record
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT
from pore2process_demo.provenance import source_manifest

from .run import FACES, HERE, OUTPUT, PULSE_S, controls, endpoint, history

PARENT = OUTPUT/'runs/20261111T194841Z-0489ef5c/summary.json'
PARENT_SHA = '186f2899c43dcbb371701814412304199a090dbc2e667c2c00276acc7ac831e2'
RESULTS = OUTPUT/'N384'
DT = 8e-5
RATIO = 4
THREADS = 10


def inputs():
    if sha(PARENT)!=PARENT_SHA:
        raise ValueError('Archived outlet study changed')
    old = validate_record(PARENT)
    runtime,info,coarse = controls()
    return old,runtime,info,coarse


def compact_mesh(directory,info,coarse,*,write):
    face = coarse_faces(CENTER.parent/'N96')
    incoming,outgoing,_ = coarse_flow(coarse,face)
    left,right,axis = (face[k]-1 for k in ('left.i32','right.i32','axis.u8'))
    seam = (axis==0)&(coarse[left,0]==95)
    q = face['balanced.f64']
    if np.any(q[seam]<=0):
        raise ValueError('The archived axial seam must be forward flowing')
    np.testing.assert_array_equal(coarse[right[seam],0],0)
    np.testing.assert_array_equal(coarse[left[seam],1:],coarse[right[seam],1:])
    before = np.zeros(incoming.shape,dtype='<i4')
    after = np.zeros_like(before)
    internal = ~seam
    before[right[internal],axis[internal]] = left[internal]+1
    after[left[internal],axis[internal]] = right[internal]+1
    arrays = {'before.i32':before,'after.i32':after,'incoming.f64':incoming,'outgoing.f64':outgoing}
    if write:
        directory.mkdir()
    for name,a in arrays.items():
        dtype = '<i4' if name.endswith('.i32') else '<f8'
        encoded = a.astype(dtype,copy=False).ravel(order='C')
        if write:
            encoded.tofile(directory/name)
        else:
            np.testing.assert_array_equal(np.fromfile(directory/name,dtype=dtype),encoded)
    dx = info['dx_m']/RATIO
    volume = dx**3
    count = len(coarse)*RATIO**3
    throughput = float(q[seam].sum())
    np.testing.assert_allclose(throughput,info['carrier_directional_flow_m3_s'][0],rtol=1e-12,atol=0)
    np.testing.assert_allclose(count*volume,len(coarse)*info['dx_m']**3,rtol=1e-14,atol=0)
    restored = q[seam]/RATIO**2*RATIO**2
    np.testing.assert_array_equal(restored,q[seam])
    # Check all 64 possible subcell offsets in every archived fluid voxel.
    rates,divergences = [],[]
    plane_error = 0.
    for a in range(3):
        lo,hi = incoming[:,a],outgoing[:,a]
        rate,net = [],[]
        for s in range(RATIO):
            f0,f1 = s/RATIO,(s+1)/RATIO
            q0 = ((1-f0)*lo+f0*hi)/RATIO**2
            q1 = ((1-f1)*lo+f1*hi)/RATIO**2
            incoming_diff = before[:,a]!=0 if s==0 else True
            outgoing_diff = after[:,a]!=0 if s==RATIO-1 else True
            degree = np.asarray(incoming_diff,dtype=int)+np.asarray(outgoing_diff,dtype=int)
            rate.append((np.maximum(-q0,0.)+np.maximum(q1,0.))/volume+ref.D/dx**2*degree)
            net.append(q1-q0)
            plane = np.bincount(coarse[:,a],weights=q1*RATIO**2,minlength=96)
            plane_error = max(plane_error,float(np.max(abs(plane-info['carrier_directional_flow_m3_s'][a]))))
        rates.append(rate)
        divergences.append(net)
    max_rate = max(float(np.max(rates[0][i]+rates[1][j]+rates[2][k])) for i,j,k in np.ndindex((RATIO,)*3))
    max_net = max(float(np.max(abs(divergences[0][i]+divergences[1][j]+divergences[2][k])))
                  for i,j,k in np.ndindex((RATIO,)*3))
    face_scale = float(max(abs(q)))/RATIO**2
    plane_relative = plane_error/throughput
    if max_net/face_scale>1e-10 or plane_relative>1e-9 or not 0<DT*max_rate<1:
        raise ValueError('N384 flow or requested monotone time step failed')
    return dict(n=384,cells=count,interior_faces=int(np.count_nonzero(internal))*RATIO**2+
        len(coarse)*3*(RATIO-1)*RATIO**2,inlet_faces=len(q[seam])*RATIO**2,
        dx_m=dx,volume_m3=count*volume,flow_m3_s=throughput,transit_s=count*volume/throughput,
        max_local_imbalance_relative=max_net/face_scale,max_plane_throughput_error_relative=plane_relative,
        max_depletion_s_inverse=max_rate,depletion_number=DT*max_rate,
        former_safety_cap_s=.8/max_rate,monotone_cap_s=1/max_rate)


def command(runtime,mode,src,dst,ratio,dx,steps,extra=()):
    return [runtime[0],f'--threads={THREADS}',*runtime[1:],str(HERE/'solve_brick.jl'),
            *map(str,[mode,src,dst,ratio,dx,ref.D,DT,steps,*extra])]


def execute(command,log,timeout):
    start = perf_counter()
    with log.open('w') as stream:
        p = subprocess.run(command,stdout=stream,stderr=subprocess.STDOUT,timeout=timeout)
    if p.returncode:
        raise RuntimeError(f'Julia failed. See {log}')
    return dict(returncode=p.returncode,wall_s=perf_counter()-start)


def comparison_test(path):
    lines = path.read_text().splitlines()
    if len(lines)!=4 or not lines[0].startswith('field,max_concentration'):
        raise ValueError('Unexpected outlet equivalence output')
    result = {}
    for line in lines[1:]:
        label,delta,relative = line.split(',')
        delta,relative = float(delta),float(relative)
        if delta!=0 or not 0<=relative<=1e-12:
            raise ValueError('Brick update did not match archived kernel')
        result[label] = dict(max_concentration_difference_kg_m3=delta,relative_exit_mass_difference=relative)
    if set(result)!={'clean','smooth','late'}:
        raise ValueError('Incomplete equivalence cases')
    return result


def summarize(base,old,mesh):
    path = base/'N384'
    h = history(path/'history.csv')
    prior = history(PARENT.parent/'N288/history.csv')
    np.testing.assert_array_equal(h['time_s'],prior['time_s'])
    np.testing.assert_allclose(h['injected_kg'],prior['injected_kg'],rtol=1e-12,atol=0)
    if not np.isfinite(h.view(np.float64)).all():
        raise ValueError('Non-finite outlet history')
    entered,out,liquid = (h[k] for k in ('injected_kg','recovered_kg','liquid_kg'))
    target = mesh['flow_m3_s']*PULSE_S
    np.testing.assert_allclose(entered,mesh['flow_m3_s']*np.minimum(h['time_s'],PULSE_S),rtol=1e-12,atol=0)
    dm = np.diff(out,prepend=0.)
    np.testing.assert_allclose(dm,h['interval_exit_kg'],rtol=0,atol=1e-14*target)
    if np.any(dm< -1e-14*target):
        raise ValueError('Unphysical negative outlet transfer')
    for key,power in (('first_moment_kg_s',1),('second_moment_kg_s2',2)):
        increments = np.diff(h[key])
        time = h['time_s']
        tolerance = 1e-14*target*time[-1]**power
        if np.any(increments<time[:-1]**power*dm[1:]-tolerance) or np.any(increments>time[1:]**power*dm[1:]+tolerance):
            raise ValueError('Inconsistent exit moment increment')
    final = np.fromfile(path/'final.f64',dtype='<f8')
    if len(final)!=mesh['cells'] or not np.isfinite(final).all():
        raise ValueError('Wrong saved final tracer field')
    np.testing.assert_allclose(final.sum()*mesh['dx_m']**3,liquid[-1],rtol=1e-12,atol=0)
    np.testing.assert_allclose([final.min(),final.max()],[h['min_kg_m3'][-1],h['max_kg_m3'][-1]],rtol=0,atol=1e-14)
    residual = float(np.max(abs(liquid+out-entered))/target)
    low,high = float(h['min_kg_m3'].min()),float(h['max_kg_m3'].max())
    fine = endpoint(h,mesh['transit_s'])
    fine.update(mass_balance_relative=residual,range_kg_m3=[low,high],remaining_kg=float(liquid[-1]),
                injected_kg=float(entered[-1]),recovered_kg=float(out[-1]))
    control = old['analysis']['cases']['288']
    tails = {}
    for n,data,stats in [('288',prior,control),('384',h,fine)]:
        tau = stats['mean_exit_s']
        recovered = data['recovered_kg']
        recovery = recovered/data['injected_kg'][-1]
        t = data['time_s']
        j = int(np.searchsorted(t,.30))
        central = (data['second_moment_kg_s2']-2*tau*data['first_moment_kg_s']+tau**2*recovered)/recovered[-1]
        tails[n] = dict(time_to_99pct_s=float(np.interp(.99,recovery,t)) if recovery[-1]>=.99 else None,
            unrecovered_fraction_at_0_30_s=float(1-recovery[j]),
            variance_from_after_0_30_s2=float(central[-1]-central[j]),
            variance_fraction_from_after_0_30_s=float((central[-1]-central[j])/central[-1]))
    mean = abs(control['mean_exit_s']/fine['mean_exit_s']-1)
    variance = abs(control['variance_s2']/fine['variance_s2']-1)
    checks = dict(same_time_grid=True,mass=residual<=1e-8,bounds=low>=-1e-8 and high<=1+1e-8,
                  recovery_and_tail=fine['tail_passed'] and control['tail_passed'],mean=mean<=.05,variance=variance<=.10)
    return dict(N288=control,N384=fine,tail=tails,mean_relative_difference=mean,
                variance_relative_difference=variance,checks=checks,short_segment_qualified=all(checks.values()))


def plan(base,runtime,info,old,mesh):
    steps = old['steps']['288']
    coarse_dx = info['dx_m']
    verify = command(runtime,'verify',base/'compact',base,3,coarse_dx/3,4,
        extra=(PARENT.parent/'N288_faces',FACES[288],PARENT.parent/'N288'))
    solve = command(runtime,'solve',base/'compact',base/'N384',RATIO,mesh['dx_m'],steps)
    return dict(verify=verify,solve=solve,steps=steps,dt_s=DT)


def run():
    old,runtime,info,coarse = inputs()
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid4().hex[:8]
    base = RESULTS/stamp
    base.mkdir(parents=True)
    files = [HERE/k for k in ('n384.py','solve_brick.jl','n384_protocol.md','run.py','solve.jl')]
    files += [ROOT/'julia/adsorption'/k for k in ('BrickOpenTracer.jl','OpenBoundedFaceTransport.jl',
              'ThreadedBoundedFaceTransport.jl','BoundedFaceTransport.jl')]
    hashes = source_manifest(files)
    for name in hashes:
        target = base/'source_snapshot'/name
        target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(ROOT/name,target)
    record = dict(schema='pore2process.N384-outlet.v1',run_id=stamp,source_sha256=hashes,
                  parent_summary=str(PARENT.relative_to(ROOT)),parent_summary_sha256=PARENT_SHA,
                  new_openlb_runs=0,new_bed_runs=0,new_DWSIM_runs=0)
    print('N384 outlet tracer study',base,flush=True)
    mesh = compact_mesh(base/'compact',info,coarse,write=True)
    requests = plan(base,runtime,info,old,mesh)
    write_json(base/'plan.json',dict(record,mesh=mesh,commands=requests))
    print('N384 face and stability audit',mesh,flush=True)
    test = execute(requests['verify'],base/'kernel_equivalence.csv',timeout=600)
    equivalence = comparison_test(base/'kernel_equivalence.csv')
    print('Archived N288 kernel equivalence',equivalence,flush=True)
    (base/'N384').mkdir()
    print('N384 outlet tracer start',requests['steps'],'steps',flush=True)
    result = execute(requests['solve'],base/'N384/stdout.log',timeout=18000)
    write_json(base/'execution.json',dict(verification=test,solve=result))
    print('N384 outlet tracer complete',result,flush=True)
    record.update(mesh=mesh,commands=requests,kernel_equivalence=equivalence,analysis=summarize(base,old,mesh))
    record['artifacts_sha256'] = {str(p.relative_to(base)):sha(p) for p in base.rglob('*') if p.is_file()}
    path = base/'summary.json'
    write_json(path,record)
    write_json(RESULTS/'latest.json',dict(summary=str(path.relative_to(ROOT))))
    print('Completed N384 comparison',path,record['analysis']['checks'],flush=True)


def check(path):
    path = path.resolve()
    record = validate_record(path)
    old,runtime,info,coarse = inputs()
    if record['parent_summary_sha256']!=PARENT_SHA:
        raise ValueError('Wrong N288 comparator')
    mesh = compact_mesh(path.parent/'compact',info,coarse,write=False)
    execution = read_json(path.parent/'execution.json')
    if (record['mesh']!=mesh or record['commands']!=plan(path.parent,runtime,info,old,mesh)
            or any(e['returncode'] for e in execution.values())
            or record['kernel_equivalence']!=comparison_test(path.parent/'kernel_equivalence.csv')
            or record['analysis']!=summarize(path.parent,old,mesh)):
        raise ValueError('N384 comparison differs on replay')
    result = dict(summary=str(path.relative_to(ROOT)),summary_sha256=sha(path),sources=len(record['source_sha256']),
                  artifacts=len(record['artifacts_sha256']),checks=record['analysis']['checks'])
    write_json(RESULTS/'verification.json',result)
    print(result)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check',type=Path)
    args=parser.parse_args()
    check(args.check) if args.check else run()
