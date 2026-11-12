"""Matched constant-background execution using the exact archived sine binary."""
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
from pore2process_demo.io import sha, write_json
from pore2process_demo.provenance import source_manifest
from pore2process_demo.runtime.openlb import verify_vendor
from examples.adsorption.transport_accuracy.run import validate_record
from .study import HERE, OUTPUT, analyze, check


def comparison(parent_path, native):
    parent=validate_record(parent_path);spec=parent['specifications'][0]
    previous=parent['cases'][spec['id']];step=previous['peak']['step']
    before=np.genfromtxt(parent_path.parent/'native'/spec['id']/'peak.csv',delimiter=',',names=True)
    control=np.genfromtxt(native/f'audit_{step}.csv',delimiter=',',names=True)
    np.testing.assert_array_equal(before[['i','j','k']],control[['i','j','k']])
    np.testing.assert_array_equal(before['material'],control['material'])
    fluid=before['material']!=2
    background=control['concentration_after_stream'][fluid]-.5
    wave=before['concentration_kg_m3'][fluid]
    diagnostic=wave-background
    a=analyze(native)
    maximum=max(abs(a['minimum_concentration_kg_m3']-.5),abs(a['maximum_concentration_kg_m3']-.5))
    return dict(sine_peak_step=step,sine_peak_time_s=previous['peak']['time_s'],
        constant_whole_history_maximum_departure_kg_m3=maximum,
        whole_history_background_bound_over_sine_excursion=maximum/previous['peak']['maximum_excursion_kg_m3'],
        constant_at_sine_peak_maximum_departure_kg_m3=float(np.max(abs(background))),
        constant_at_sine_peak_RMS_departure_kg_m3=float(np.sqrt(np.mean(background**2))),
        diagnostic_background_removed_minimum_kg_m3=float(np.min(diagnostic)),
        diagnostic_background_removed_maximum_kg_m3=float(np.max(diagnostic)),
        production_correction_applied=False)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--parent',type=Path)
    parser.add_argument('--check',type=Path)
    args=parser.parse_args()
    if args.check:
        s=check(args.check);spec=s['specifications'][0]
        if comparison(ROOT/s['parent_summary'],args.check.resolve().parent/'native'/spec['id'])!=s['background_comparison']:
            raise ValueError('Constant-background comparison does not replay')
        return
    if args.parent is None:parser.error('Require --parent')
    parent_path=args.parent.resolve();parent=validate_record(parent_path)
    if len(parent['specifications'])!=1:raise ValueError('Require a single-case parent')
    old=parent['specifications'][0];name=old['id'];previous=parent['cases'][name]
    if previous['peak']['initial_amplitude_kg_m3']!=.25:raise ValueError('Require the original sinusoidal parent')
    request=json.loads((parent_path.parent/'native'/name/'request.json').read_text())
    command=list(request['command']);binary_parent=Path(command[0])
    if sha(binary_parent)!=request['binary_sha256']:raise ValueError('Parent executable changed')
    spec=dict(old,id=f"n{old['n']}_constant_full",amplitude=0.,audit_step=previous['peak']['step'])
    command[13]=str(spec['audit_step']);command[14]='0'
    rid=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid4().hex[:8]
    base=OUTPUT/'constant/runs'/rid;base.mkdir(parents=True)
    binary=base/'binary/capillary3d';binary.parent.mkdir();shutil.copy2(binary_parent,binary);command[0]=str(binary)
    sources=source_manifest([Path(__file__),HERE/'constant_protocol.md',HERE/'study.py',
        ROOT/'examples/adsorption/transport_accuracy/profile.py',ROOT/'examples/adsorption/transport_accuracy/run.py'])
    for name in sources:
        target=base/'source_snapshot'/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/name,target)
    record=dict(schema='pore2process.matched-constant-background.v1',run_id=rid,
        source_sha256=sources,specifications=[spec],parent_summary=str(parent_path.relative_to(ROOT)),parent_summary_sha256=sha(parent_path),
        hydraulic_summary=parent['hydraulic_summary'],hydraulic_summary_sha256=parent['hydraulic_summary_sha256'],
        new_native_runs=1,new_DWSIM_runs=0,new_adsorption_predictions=0,service_prediction_qualified=False,
        binary_sha256=sha(binary),OMP_NUM_THREADS=parent.get('OMP_NUM_THREADS',1))
    if parent.get('refined_summary'):
        record.update(refined_summary=parent['refined_summary'],refined_summary_sha256=parent['refined_summary_sha256'])
        shutil.copy2(parent_path.parent/'source_point.json',base/'source_point.json')
    dest=base/'native'/spec['id'];dest.mkdir(parents=True)
    write_json(dest/'request.json',dict(spec=spec,command=command,binary_sha256=sha(binary),native_field_sha256=sha(Path(spec['field']))))
    write_json(base/'plan.json',record)
    env=dict(os.environ,OMP_NUM_THREADS=str(record['OMP_NUM_THREADS']),OMP_DYNAMIC='FALSE')
    print('Matched constant starts',base,flush=True);started=perf_counter()
    with (dest/'stdout.log').open('w') as stream:
        execution=subprocess.run(command,cwd=dest,env=env,stdout=stream,stderr=subprocess.STDOUT,timeout=7200)
    write_json(dest/'execution.json',dict(returncode=execution.returncode,elapsed_s=perf_counter()-started))
    if execution.returncode:raise RuntimeError('Native constant control failed')
    record['cases']={spec['id']:analyze(dest)}
    r=record['cases'][spec['id']]['raw'];oldraw=previous['raw']
    if r['dt_s']!=oldraw['dt_s'] or r['steps']!=oldraw['steps']:raise ValueError('Time step/window changed')
    record['background_comparison']=comparison(parent_path,dest)
    record['vendor']=verify_vendor()
    record['artifacts_sha256']={str(p.relative_to(base)):sha(p) for p in base.rglob('*') if p.is_file()}
    write_json(base/'summary.json',record)
    write_json(OUTPUT/'constant_latest.json',dict(summary=str((base/'summary.json').relative_to(ROOT))))
    print('Matched constant completed',base,round(perf_counter()-started,2),'s',flush=True)


if __name__=='__main__':main()
