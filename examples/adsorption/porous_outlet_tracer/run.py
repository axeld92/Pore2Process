"""Run and replay a one-cell pulse on the archived TRT-derived face flows."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
import shutil
import subprocess
from time import perf_counter
from uuid import uuid4

import numpy as np

from examples.adsorption.porous_face_handoff.run import archived
from examples.adsorption.porous_tracer_refinement.mesh import FILES, add_outgoing, coarse_faces
from examples.adsorption.porous_tracer_refinement.n288 import RESULTS as N288_RESULTS
from examples.adsorption.porous_tracer_refinement.n192_time import PARENT, inputs
from examples.adsorption.porous_tracer_refinement.time_half import CENTER
from examples.adsorption.stationary_curved_wall import reference as ref
from examples.adsorption.transport_accuracy.run import validate_record
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT
from pore2process_demo.provenance import source_manifest

HERE = Path(__file__).parent
OUTPUT = ROOT/'results/adsorption/porous_outlet_tracer'
N288 = N288_RESULTS/'20261111T164622Z-a501da76/summary.json'
N288_SHA = 'b254edb0c5e995a10e24bcaf81a69a3253e92cb3520fbadb85095e8fda80513e'
FACES = {96:CENTER.parent/'N96',192:PARENT.parent/'N192_faces',288:N288.parent/'N288_faces'}
DT = {96:2e-4,192:8e-5,288:8e-5}
THREADS = {96:4,192:4,288:8}
PULSE_S = .004


def controls():
    if sha(N288)!=N288_SHA:
        raise ValueError('Saved N288 flow record differs')
    old = validate_record(N288)
    if not old['analysis']['checks']['flow_preserved']:
        raise ValueError('Saved N288 flow failed conservation')
    _,runtime,_ = inputs()
    info,coarse,_ = archived(96)
    return runtime,info,coarse


def coordinates(index,coarse,ratio):
    parent,child = np.divmod(index,ratio**3)
    offsets = np.indices((ratio,)*3,dtype=np.int32).reshape(3,-1).T
    return ratio*coarse[parent]+offsets[child]


def open_faces(n,directory,info,coarse,*,write):
    source = FACES[n]
    ratio = n//96
    count = len(coarse)*ratio**3
    dx = info['dx_m']/ratio
    faces = coarse_faces(source)
    left,right,axis,q = (faces[k] for k in FILES)
    axial = np.flatnonzero(axis==1)
    cell = left[axial]-1
    x = ratio*coarse[cell//ratio**3,0]+(cell//ratio**2)%ratio
    seam = axial[x==n-1]
    a,b = coordinates(left[seam]-1,coarse,ratio),coordinates(right[seam]-1,coarse,ratio)
    np.testing.assert_array_equal(a[:,1:],b[:,1:])
    np.testing.assert_array_equal(a[:,0],n-1)
    np.testing.assert_array_equal(b[:,0],0)
    if np.any(q[seam]<=0) or len(np.unique(left[seam]))!=len(seam) or len(np.unique(right[seam]))!=len(seam):
        raise ValueError('Expected positive flow through distinct axial faces')
    target = info['carrier_directional_flow_m3_s'][0]
    throughput = float(q[seam].sum())
    np.testing.assert_allclose(throughput,target,rtol=1e-12,atol=0)
    boundary = np.zeros(count)
    boundary[right[seam]-1] = q[seam]/dx**3
    boundary[left[seam]-1] = -q[seam]/dx**3
    keep = np.ones(len(q),dtype=bool)
    keep[seam] = False
    internal = {name:value[keep] for name,value in faces.items()}
    internal['boundary_rate.f64'] = boundary
    if write:
        directory.mkdir(parents=True)
    for name,values in internal.items():
        if write:
            values.tofile(directory/name)
        else:
            np.testing.assert_array_equal(np.fromfile(directory/name,dtype=values.dtype),values)
    l,r,f = (internal[k] for k in ('left.i32','right.i32','balanced.f64'))
    l,r = l-1,r-1
    net = np.bincount(l,weights=f,minlength=count)-np.bincount(r,weights=f,minlength=count)
    net -= boundary*dx**3
    residual = float(max(abs(net))/max(abs(f)))
    if residual>1e-10:
        raise ValueError('Opening the seam changed local carrier continuity')
    outgoing = np.maximum(-boundary,0.)
    add_outgoing(outgoing,l,r,f,dx)
    cap = .8/float(max(outgoing))
    if DT[n]>cap:
        raise ValueError('Time step exceeds the stable cap')
    volume = count*dx**3
    return dict(n=n,cells=count,interior_faces=len(f),inlet_faces=len(seam),outlet_faces=len(seam),
                dx_m=dx,volume_m3=volume,flow_m3_s=throughput,source_flow_m3_s=target,
                minimum_boundary_flow_m3_s=float(q[seam].min()),max_local_imbalance_relative=residual,
                stability_cap_s=cap,transit_s=volume/throughput,
                saved_face_flow_sha256=sha(source/'balanced.f64'))


def cmd(base,n,mesh,runtime,steps):
    dt = DT[n]
    pulse_steps = round(PULSE_S/dt)
    if abs(pulse_steps*dt-PULSE_S)>1e-16:
        raise ValueError('Pulse must end exactly on a time step')
    return [runtime[0],f'--threads={THREADS[n]}',*runtime[1:],str(HERE/'solve.jl'),
            *map(str,[base/f'N{n}_faces',base/f'N{n}',mesh['dx_m'],ref.D,dt,pulse_steps,steps,THREADS[n]])]


def execute(base,n,mesh,runtime,steps):
    directory = base/f'N{n}'
    directory.mkdir()
    command = cmd(base,n,mesh,runtime,steps)
    write_json(directory/'request.json',dict(command=command,steps=steps,dt_s=DT[n],
        open_face_sha256=sha(base/f'N{n}_faces/balanced.f64'),
        boundary_sha256=sha(base/f'N{n}_faces/boundary_rate.f64')))
    print('Outlet tracer start',n,steps,'steps',flush=True)
    start = perf_counter()
    with (directory/'stdout.log').open('w') as log:
        p = subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,timeout=9000)
    write_json(directory/'execution.json',dict(returncode=p.returncode,wall_s=perf_counter()-start))
    if p.returncode:
        raise RuntimeError(f'Julia outlet case failed. See {directory / "stdout.log"}')
    print('Outlet tracer complete',n,round(perf_counter()-start,1),'s',flush=True)


def history(path):
    return np.genfromtxt(path,delimiter=',',names=True)


def endpoint(data,transit_s):
    time = data['time_s']
    out = data['recovered_kg']
    injected = float(data['injected_kg'][-1])
    if injected<=0 or out[-1]<=0:
        raise ValueError('No tracer entered or exited')
    k = np.searchsorted(time,time[-1]-transit_s)
    m1,m2 = data['first_moment_kg_s'],data['second_moment_kg_s2']
    tail = [float((v[-1]-v[k])/v[-1]) for v in (m1,m2)]
    recovery = float(out[-1]/injected)
    mean = float(m1[-1]/out[-1])
    variance = float(m2[-1]/out[-1]-mean**2)
    early = {str(f):float(np.interp(f,out/injected,time)) if recovery>=f else None for f in (.05,.1)}
    return dict(time_s=float(time[-1]),recovery=recovery,mean_exit_s=mean,variance_s2=variance,
                early_arrival_s=early,last_transit_raw_moment_fractions=tail,
                tail_passed=recovery>=.999 and max(tail)<=.01)


def checkpoints(tau,dt):
    return [int(np.ceil(i*tau/dt)) for i in range(1,11)]


def observation(scout,mesh):
    cp = checkpoints(mesh['transit_s'],DT[96])
    for i,k in enumerate(cp[3:],start=4):
        if endpoint(scout[:k+1],mesh['transit_s'])['tail_passed']:
            return min(10,i+1)
    return 10


def analyze_case(base,n,mesh,runtime,steps):
    path = base/f'N{n}'
    req = dict(command=cmd(base,n,mesh,runtime,steps),steps=steps,dt_s=DT[n],
               open_face_sha256=sha(base/f'N{n}_faces/balanced.f64'),
               boundary_sha256=sha(base/f'N{n}_faces/boundary_rate.f64'))
    if read_json(path/'request.json')!=req or read_json(path/'execution.json')['returncode']:
        raise ValueError('Archived simulation request differs')
    h = history(path/'history.csv')
    if not np.isfinite(h.view(np.float64)).all():
        raise ValueError('Non-finite outlet history')
    np.testing.assert_allclose(h['time_s'],np.arange(steps+1)*DT[n],rtol=1e-14,atol=0)
    injected,recovered,inventory = (h[k] for k in ('injected_kg','recovered_kg','liquid_kg'))
    target = mesh['flow_m3_s']*PULSE_S
    np.testing.assert_allclose(injected,mesh['flow_m3_s']*np.minimum(h['time_s'],PULSE_S),rtol=1e-12,atol=0)
    dm = np.diff(recovered,prepend=0.)
    np.testing.assert_allclose(dm,h['interval_exit_kg'],rtol=0,atol=1e-14*target)
    if np.any(dm < -1e-14*target):
        raise ValueError('Unphysical negative outlet transfer')
    for key,power in (('first_moment_kg_s',1),('second_moment_kg_s2',2)):
        increment = np.diff(h[key])
        scale = target*h['time_s'][-1]**power
        lo = h['time_s'][:-1]**power*dm[1:]
        hi = h['time_s'][1:]**power*dm[1:]
        if np.any(increment<lo-1e-14*scale) or np.any(increment>hi+1e-14*scale):
            raise ValueError('Exit moment increment inconsistent with mass transfer')
    residual = float(max(abs(inventory+recovered-injected))/target)
    final = np.fromfile(path/'final.f64',dtype='<f8')
    if len(final)!=mesh['cells'] or not np.isfinite(final).all():
        raise ValueError('Invalid final tracer field')
    np.testing.assert_allclose([final.min(),final.max()],
        [h['min_kg_m3'][-1],h['max_kg_m3'][-1]],rtol=0,atol=1e-14)
    np.testing.assert_allclose(final.sum()*mesh['dx_m']**3,inventory[-1],rtol=1e-12,atol=1e-20*target)
    low,high = float(h['min_kg_m3'].min()),float(h['max_kg_m3'].max())
    metrics = endpoint(h,mesh['transit_s'])
    metrics.update(injected_kg=float(injected[-1]),recovered_kg=float(recovered[-1]),
        remaining_kg=float(inventory[-1]),mass_balance_relative=residual,range_kg_m3=[low,high],
        mass_passed=residual<=1e-8,bounds_passed=low>=-1e-8 and high<=1+1e-8)
    return metrics,h


def compare(base,meshes,runtime,steps):
    cases = {str(n):analyze_case(base,n,meshes[str(n)],runtime,steps[str(n)])[0] for n in FACES}
    a,b = cases['192'],cases['288']
    if a['time_s']!=b['time_s']:
        raise ValueError('Comparison windows differ')
    mean_change = abs(a['mean_exit_s']/b['mean_exit_s']-1)
    variance_change = abs(a['variance_s2']/b['variance_s2']-1)
    early_change = {key:a['early_arrival_s'][key]-value if a['early_arrival_s'][key] is not None and value is not None
                    else None for key,value in b['early_arrival_s'].items()}
    checks = dict(mass=all(c['mass_passed'] for c in (a,b)),bounds=all(c['bounds_passed'] for c in (a,b)),
        recovery_and_tail=all(c['tail_passed'] for c in (a,b)),mean=mean_change<=.05,variance=variance_change<=.10)
    return dict(cases=cases,mean_relative_difference=mean_change,variance_relative_difference=variance_change,
        N192_minus_N288_early_arrival_s=early_change,checks=checks,short_segment_qualified=all(checks.values()))


def run():
    runtime,info,coarse = controls()
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid4().hex[:8]
    base = OUTPUT/'runs'/stamp
    base.mkdir(parents=True)
    files = [HERE/k for k in ('run.py','solve.jl','test_open.jl','protocol.md')]
    files += [ROOT/'julia/adsorption'/k for k in
              ('BoundedFaceTransport.jl','ThreadedBoundedFaceTransport.jl','OpenBoundedFaceTransport.jl')]
    hashes = source_manifest(files)
    for name in hashes:
        target = base/'source_snapshot'/name
        target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(ROOT/name,target)
    record = dict(schema='pore2process.short-outlet-tracer.v1',run_id=stamp,source_sha256=hashes,
        parent_N288_sha256=N288_SHA,new_openlb_runs=0,new_bed_runs=0,new_DWSIM_runs=0)
    write_json(base/'plan.json',record)
    print('Short outlet tracer study',base,flush=True)
    test = [runtime[0],'--threads=4',*runtime[1:],str(HERE/'test_open.jl')]
    with (base/'boundary_test.log').open('w') as log:
        subprocess.run(test,stdout=log,stderr=subprocess.STDOUT,timeout=120,check=True)
    write_json(base/'boundary_test.json',dict(command=test,returncode=0))
    meshes = {'96':open_faces(96,base/'N96_faces',info,coarse,write=True)}
    steps = {'96':checkpoints(meshes['96']['transit_s'],DT[96])[-1]}
    with ThreadPoolExecutor(max_workers=2) as pool:
        scout = pool.submit(execute,base,96,meshes['96'],runtime,steps['96'])
        for n in (192,288):
            meshes[str(n)] = open_faces(n,base/f'N{n}_faces',info,coarse,write=True)
            print('Open mesh audit',meshes[str(n)],flush=True)
        scout.result()
        scout_metrics,h = analyze_case(base,96,meshes['96'],runtime,steps['96'])
        if not scout_metrics['mass_passed'] or not scout_metrics['bounds_passed']:
            raise ValueError('Scout failed mass or concentration checks')
        windows = observation(h,meshes['96'])
        for n in (192,288):
            steps[str(n)] = checkpoints(meshes[str(n)]['transit_s'],DT[n])[windows-1]
        if steps['192']!=steps['288']:
            raise ValueError('The grids must use a matched window')
        write_json(base/'schedule.json',dict(nominal_transits=windows,meshes=meshes,steps=steps))
        print('Matched window',windows,'transits',steps['288'],'steps',flush=True)
        for future in [pool.submit(execute,base,n,meshes[str(n)],runtime,steps[str(n)]) for n in (192,288)]:
            future.result()
    record.update(meshes=meshes,steps=steps,nominal_transits=windows,analysis=compare(base,meshes,runtime,steps))
    record['artifacts_sha256'] = {str(p.relative_to(base)):sha(p) for p in base.rglob('*') if p.is_file()}
    path = base/'summary.json'
    write_json(path,record)
    write_json(OUTPUT/'latest.json',dict(summary=str(path.relative_to(ROOT))))
    print('Completed outlet tracer study',path,record['analysis']['checks'],flush=True)


def check(path):
    path = path.resolve()
    saved = validate_record(path)
    runtime,info,coarse = controls()
    if saved['parent_N288_sha256']!=N288_SHA:
        raise ValueError('Parent record differs')
    base = path.parent
    meshes = {str(n):open_faces(n,base/f'N{n}_faces',info,coarse,write=False) for n in FACES}
    steps = saved['steps']
    _,h = analyze_case(base,96,meshes['96'],runtime,steps['96'])
    windows = observation(h,meshes['96'])
    expected = {str(n):checkpoints(meshes[str(n)]['transit_s'],DT[n])[9 if n==96 else windows-1] for n in FACES}
    if (meshes!=saved['meshes'] or steps!=expected or windows!=saved['nominal_transits']
            or compare(base,meshes,runtime,steps)!=saved['analysis']):
        raise ValueError('Outlet tracer replay differs')
    result = dict(summary=str(path.relative_to(ROOT)),summary_sha256=sha(path),sources=len(saved['source_sha256']),
                  artifacts=len(saved['artifacts_sha256']),checks=saved['analysis']['checks'])
    write_json(OUTPUT/'verification.json',result)
    print(result)


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check',type=Path)
    args = parser.parse_args()
    check(args.check) if args.check else run()
