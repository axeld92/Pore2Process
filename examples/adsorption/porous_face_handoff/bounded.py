"""Test bounded tracer spreading using the exact saved balanced face flows."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
import shutil
import subprocess
from time import perf_counter
from uuid import uuid4

import numpy as np

from examples.adsorption.stationary_curved_wall import reference as ref
from examples.adsorption.transport_accuracy.run import validate_record
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT
from pore2process_demo.provenance import source_manifest

from .run import HERE, OUTPUT, SNAPSHOTS, archived, scalar_metrics

JOBS = ((64,'varying',1.),(64,'half_step',.5),(96,'varying',1.),(96,'uniform',1.))


def execute(base, source, parent, job):
    n, name, factor = job
    directory = base/f'N{n}'/name
    directory.mkdir(parents=True)
    src = source.parent/f'N{n}'
    dx = parent['analysis']['cases'][str(n)]['source']['dx_m']
    args = [dx,ref.D,*ref.TIMES,name,factor]
    cmd = [*parent['julia_runtime'],str(HERE/'solve_bounded.jl'),str(src),str(directory),*map(str,args)]
    write_json(directory/'request.json',dict(command=cmd,balanced_sha256=sha(src/'balanced.f64')))
    print('Bounded Julia start',n,name,flush=True)
    start = perf_counter()
    with (directory/'stdout.log').open('w') as log:
        p = subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,timeout=1800)
    write_json(directory/'execution.json',dict(returncode=p.returncode,wall_s=perf_counter()-start))
    if p.returncode:
        raise RuntimeError(f'Bounded Julia transport failed. See {directory / "stdout.log"}')
    print('Bounded Julia complete',n,name,round(perf_counter()-start,1),'s',flush=True)


def analyze(base, source, parent):
    cases, fields = {}, {}
    for n in (64,96):
        info, ijk, _ = archived(n)
        initial = np.fromfile(source.parent/f'N{n}'/'varying_initial.f64',dtype='<f8')
        for _, name, factor in (job for job in JOBS if job[0]==n):
            directory = base/f'N{n}'/name
            src = source.parent/f'N{n}'
            args = [info['dx_m'],ref.D,*ref.TIMES,name,factor]
            cmd = [*parent['julia_runtime'],str(HERE/'solve_bounded.jl'),str(src),str(directory),*map(str,args)]
            if (read_json(directory/'request.json') != dict(command=cmd,balanced_sha256=sha(src/'balanced.f64'))
                    or read_json(directory/'execution.json')['returncode']):
                raise ValueError('Saved execution differs')
            u = initial if name!='uniform' else np.full(len(initial),.5)
            s, f = scalar_metrics(directory,name,u,ijk,n)
            for time, metrics in s['snapshots'].items():
                loss = .5*(float(np.var(initial))-metrics['variance_kg2_m6'])
                physical = metrics['integrated_dissipation']
                share = (loss-physical)/loss if loss>0 else None
                metrics['net_discretization_fraction_of_variance_loss'] = share
                s['checks'].pop(time+'_variance_budget')
                s['checks'][time+'_numerical_spreading'] = share is not None and abs(share)<=.15
            key = f'N{n}_{name}'
            cases[key], fields[key] = s, f
    space = {}
    time = {}
    for i,t in enumerate(SNAPSHOTS):
        a,b = (cases[f'N{n}_varying']['snapshots'][t] for n in (64,96))
        variance = abs(a['variance_kg2_m6']/b['variance_kg2_m6']-1)
        amplitude = abs(a['fundamental_amplitude_kg_m3']/b['fundamental_amplitude_kg_m3']-1)
        phase = abs(float(np.angle(np.exp(1j*(a['phase_rad']-b['phase_rad'])))))
        space[t] = dict(variance_relative_change=variance,amplitude_relative_change=amplitude,phase_change_rad=phase,
                        passed=variance<=.05 and amplitude<=.05 and phase<=.05)
        delta = (fields['N64_varying'][i]-fields['N64_half_step'][i])/ref.A
        time[t] = float(np.sqrt(np.mean(delta**2)))
    checks = dict(handoff=parent['analysis']['handoff_qualified'],
                  uniform=cases['N96_uniform']['checks']['constant'],
                  bounded_and_conserved=all(c['checks']['physical_bounds'] and c['checks']['mass'] for c in cases.values()),
                  spatial_spreading=all(s['passed'] for s in space.values()),
                  time_refinement=max(time.values())<=.001,
                  numerical_spreading=all(v for c in cases.values() for k,v in c['checks'].items() if k.endswith('_numerical_spreading')))
    return dict(cases=cases,spatial_change=space,time_change_RMS_over_amplitude=time,
                checks=checks,periodic_control_qualified=all(checks.values()))


def run(source):
    source = source.resolve()
    parent = validate_record(source)
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid4().hex[:8]
    base = OUTPUT/'bounded'/stamp
    base.mkdir(parents=True)
    files = [HERE/k for k in ('bounded.py','solve_bounded.jl','bounded_protocol.md','test_bounded.jl')]
    files.append(ROOT/'julia/adsorption/BoundedFaceTransport.jl')
    hashes = source_manifest(files)
    for name in hashes:
        target = base/'source_snapshot'/name
        target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(ROOT/name,target)
    record=dict(schema='pore2process.porous-face-bounded.v1',run_id=stamp,source_sha256=hashes,
                centered_summary=str(source.relative_to(ROOT)),centered_summary_sha256=sha(source),
                new_openlb_runs=0,new_bed_runs=0,new_DWSIM_runs=0)
    write_json(base/'plan.json',record)
    with ThreadPoolExecutor(max_workers=2) as pool:
        for future in [pool.submit(execute,base,source,parent,job) for job in JOBS]:
            future.result()
    record['analysis']=analyze(base,source,parent)
    record['artifacts_sha256']={str(p.relative_to(base)):sha(p) for p in base.rglob('*') if p.is_file()}
    path=base/'summary.json'
    write_json(path,record)
    write_json(OUTPUT/'latest_bounded.json',dict(summary=str(path.relative_to(ROOT))))
    print('Bounded comparison',path,record['analysis']['checks'],flush=True)


def check(path):
    path = path.resolve()
    record = validate_record(path)
    source = ROOT/record['centered_summary']
    if sha(source)!=record['centered_summary_sha256']:
        raise ValueError('Parent centered archive differs')
    parent = validate_record(source)
    analysis = analyze(path.parent,source,parent)
    if analysis!=record['analysis']:
        raise ValueError('Replay differs')
    report=dict(summary=str(path.relative_to(ROOT)),summary_sha256=sha(path),sources=len(record['source_sha256']),
                artifacts=len(record['artifacts_sha256']),checks=analysis['checks'])
    write_json(OUTPUT/'bounded_verification.json',report)
    print(report)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--centered',type=Path)
    parser.add_argument('--check',type=Path)
    args=parser.parse_args()
    check(args.check) if args.check else run(args.centered or ROOT/read_json(OUTPUT/'latest.json')['summary'])
