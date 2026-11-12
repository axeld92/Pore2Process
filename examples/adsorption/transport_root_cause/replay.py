"""Deterministic peak or half-time-step trial using an archived native parent."""
import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import shutil
import subprocess
from time import perf_counter
from uuid import uuid4

import numpy as np

from pore2process_demo.paths import ROOT
from pore2process_demo.io import sha,write_json
from pore2process_demo.provenance import source_manifest
from pore2process_demo.runtime.openlb import verify_vendor
from examples.adsorption.transport_accuracy.run import validate_record
from .study import HERE,OUTPUT,analyze,check


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--parent',type=Path)
    p.add_argument('--check',type=Path)
    mode=p.add_mutually_exclusive_group()
    mode.add_argument('--peak-audit',action='store_true')
    mode.add_argument('--half-dt',action='store_true')
    args=p.parse_args()
    if args.check:check(args.check);return
    if args.parent is None or not (args.peak_audit or args.half_dt):p.error('Require --parent and either --peak-audit or --half-dt')
    parent_path=args.parent.resolve();parent=validate_record(parent_path)
    if len(parent['specifications'])!=1:raise ValueError('Choose a single-case native parent')
    previous=parent['specifications'][0];name=previous['id'];old=parent['cases'][name]
    request=json.loads((parent_path.parent/'native'/name/'request.json').read_text())
    command=list(request['command']);raw=old['raw'];spec=dict(previous)
    if args.half_dt:
        dt=raw['dt_s']/2
        command[5]=repr(.5+3e-9*dt/raw['dx_m']**2)
        command[6]=repr(raw['steps']*raw['dt_s']-.01*dt)
        command[13]='-1';spec.update(id=name+'_half_dt',audit_step=-1,end_s=float(command[6]),requested_dt_s=dt)
    else:
        step=old['peak']['step']
        if step<1:raise ValueError('No evolved peak to audit')
        command[6]=repr((step-.01)*raw['dt_s']);command[13]=str(step)
        spec.update(id=name+f'_audit{step}',audit_step=step,end_s=float(command[6]))
    rid=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid4().hex[:8]
    base=OUTPUT/'replay/runs'/rid;base.mkdir(parents=True)
    binary=base/'binary/capillary3d';binary.parent.mkdir()
    source_binary=Path(request['command'][0])
    if sha(source_binary)!=request['binary_sha256']:raise ValueError('Native parent executable changed; no silent rebuild')
    shutil.copy2(source_binary,binary);command[0]=str(binary)
    files=[Path(__file__),HERE/'study.py',HERE/'protocol.md',ROOT/'examples/adsorption/transport_accuracy/profile.py',
           ROOT/'examples/adsorption/transport_accuracy/run.py',ROOT/'examples/adsorption/column_transport/run.py']
    sources=source_manifest(files)
    for file in sources:
        target=base/'source_snapshot'/file;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/file,target)
    record=dict(schema='pore2process.transport-root-native-replay.v1',run_id=rid,source_sha256=sources,
        specifications=[spec],parent_summary=str(parent_path.relative_to(ROOT)),parent_summary_sha256=sha(parent_path),
        hydraulic_summary=parent['hydraulic_summary'],hydraulic_summary_sha256=parent['hydraulic_summary_sha256'],
        mode='half_dt' if args.half_dt else 'peak_audit',new_native_runs=1,new_DWSIM_runs=0,
        new_adsorption_predictions=0,service_prediction_qualified=False,
        binary_sha256=sha(binary),OMP_NUM_THREADS=parent.get('OMP_NUM_THREADS',1),
        note='Same preserved executable and source-point snapshot; only explicit time/diagnostic inputs change.')
    if parent.get('refined_summary'):
        record.update(refined_summary=parent['refined_summary'],refined_summary_sha256=parent['refined_summary_sha256'])
        shutil.copy2(parent_path.parent/'source_point.json',base/'source_point.json')
    dest=base/'native'/spec['id'];dest.mkdir(parents=True)
    write_json(dest/'request.json',dict(spec=spec,command=command,binary_sha256=sha(binary),native_field_sha256=sha(Path(spec['field']))))
    write_json(base/'plan.json',record)
    env=dict(os.environ,OMP_NUM_THREADS=str(record['OMP_NUM_THREADS']),OMP_DYNAMIC='FALSE')
    print('Native replay starts',base,flush=True);started=perf_counter()
    with (dest/'stdout.log').open('w') as stream:
        result=subprocess.run(command,cwd=dest,env=env,stdout=stream,stderr=subprocess.STDOUT,timeout=7200)
    write_json(dest/'execution.json',dict(returncode=result.returncode,elapsed_s=perf_counter()-started))
    if result.returncode:raise RuntimeError('Native replay failed: '+str(dest/'stdout.log'))
    record['cases']={spec['id']:analyze(dest)}
    if args.peak_audit:
        old_field=np.genfromtxt(parent_path.parent/'native'/name/'peak.csv',delimiter=',',names=True)
        replay=np.genfromtxt(dest/'final.csv',delimiter=',',names=True)
        np.testing.assert_array_equal(old_field[['i','j','k']],replay[['i','j','k']])
        error=float(np.max(abs(old_field['concentration_kg_m3']-replay['concentration_kg_m3'])))
        if error>1e-12:raise ValueError('Prefix field does not reproduce the parent peak')
        record['peak_field_maximum_absolute_reproduction_error_kg_m3']=error
    else:
        previous_field=np.genfromtxt(parent_path.parent/'native'/name/'final.csv',delimiter=',',names=True)
        current_field=np.genfromtxt(dest/'final.csv',delimiter=',',names=True)
        fluid=current_field['material']!=2
        difference=current_field['concentration_kg_m3'][fluid]-previous_field['concentration_kg_m3'][fluid]
        newraw=record['cases'][spec['id']]['raw']
        if abs(newraw['steps']*newraw['dt_s']/(raw['steps']*raw['dt_s'])-1)>1e-9:raise ValueError('Time-control physical window differs')
        record['time_comparison']=dict(final_fluid_RMS_difference_over_initial_amplitude=float(np.sqrt(np.mean(difference**2))/.25),
             final_fluid_max_difference_over_initial_amplitude=float(np.max(abs(difference))/.25),
             actual_duration_relative_difference=abs(newraw['steps']*newraw['dt_s']/(raw['steps']*raw['dt_s'])-1))
    record['vendor']=verify_vendor()
    record['artifacts_sha256']={str(x.relative_to(base)):sha(x) for x in base.rglob('*') if x.is_file()}
    write_json(base/'summary.json',record)
    write_json(OUTPUT/(record['mode']+'_latest.json'),dict(summary=str((base/'summary.json').relative_to(ROOT))))
    print('Native replay completed',base,round(perf_counter()-started,2),'s',flush=True)


if __name__=='__main__':main()
