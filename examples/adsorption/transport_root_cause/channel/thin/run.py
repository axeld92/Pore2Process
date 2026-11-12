"""Cheaper higher-resolution shear controls after a full-cube equivalence check."""
from datetime import datetime,timezone
import argparse
import json
from pathlib import Path
import shutil
import subprocess
from time import perf_counter
from uuid import uuid4
import numpy as np
from pore2process_demo.paths import ROOT
from pore2process_demo.io import sha,write_json
from pore2process_demo.provenance import source_manifest
from pore2process_demo.runtime.openlb import build_case,verify_vendor
from ..run import PHYSICAL,OUTPUT as PARENT_OUTPUT,analyze

HERE=Path(__file__).parent
OUTPUT=PARENT_OUTPUT/'thin'
DURATION=11294*PHYSICAL['dt_s']


def specifications():
    return [dict(id='thin32',n=32,dt_factor=1.),dict(id='thin96',n=96,dt_factor=1.),
            dict(id='thin128',n=128,dt_factor=1.),dict(id='thin128_half_dt',n=128,dt_factor=.5)]


def execute(base,spec,binary):
    dest=base/'native'/spec['id'];dest.mkdir(parents=True);p=PHYSICAL
    command=[str(binary),str(spec['n']),str(p['side_m']),str(p['diffusivity_m2_s']),str(p['dt_s']*spec['dt_factor']),
             str(DURATION),str(p['speed_m_s']),'.75','.25']
    write_json(dest/'request.json',dict(spec=spec,physical={**p,'end_s':DURATION,'dt_s':p['dt_s']*spec['dt_factor']},actual_target_duration_s=DURATION,command=command,binary_sha256=sha(binary)))
    print('Native thin start',spec['id'],flush=True);started=perf_counter()
    with (dest/'stdout.log').open('w') as stream:
        result=subprocess.run(command,cwd=dest,stdout=stream,stderr=subprocess.STDOUT,timeout=3600)
    write_json(dest/'execution.json',dict(returncode=result.returncode,elapsed_s=perf_counter()-started))
    if result.returncode:raise RuntimeError('Native thin run failed: '+str(dest/'stdout.log'))
    print('Native thin finished',spec['id'],f'{perf_counter()-started:.1f}s',flush=True)


def equivalence(full,thin):
    full=Path(full);thin=Path(thin);comparisons={}
    raw_full=json.loads((full/'raw.json').read_text());raw_thin=json.loads((thin/'raw.json').read_text())
    for name in ('initial.csv','peak.csv','final.csv'):
        a=np.genfromtxt(full/name,delimiter=',',names=True);b=np.genfromtxt(thin/name,delimiter=',',names=True)
        if a.shape!=b.shape:raise ValueError('Thin equivalence shape mismatch.')
        np.testing.assert_array_equal(a['i'],b['i']);np.testing.assert_array_equal(a['j'],b['j'])
        for key in ('x_m','y_m','ux_m_s'):np.testing.assert_allclose(a[key],b[key],atol=1e-15,rtol=1e-13)
        comparisons[name]=float(np.max(abs(a['concentration_kg_m3']-b['concentration_kg_m3'])))
        comparisons[name+'_maximum_z_variation_kg_m3']=float(np.max(b['z_variation_kg_m3']))
    a=np.genfromtxt(full/'inventory.csv',delimiter=',',names=True);b=np.genfromtxt(thin/'inventory.csv',delimiter=',',names=True)
    np.testing.assert_array_equal(a['step'],b['step']);np.testing.assert_allclose(a['time_s'],b['time_s'],atol=1e-15,rtol=1e-13)
    comparisons['all_step_minimum_error_kg_m3']=float(np.max(abs(a['min_concentration_kg_m3']-b['min_concentration_kg_m3'])))
    comparisons['all_step_maximum_error_kg_m3']=float(np.max(abs(a['max_concentration_kg_m3']-b['max_concentration_kg_m3'])))
    comparisons['normalized_inventory_history_max_difference']=float(np.max(abs(a['inventory_kg']/a['inventory_kg'][0]-b['inventory_kg']/b['inventory_kg'][0])))
    comparisons['initial_inventory_depth_scaling_relative_error']=abs(float(b['inventory_kg'][0]/a['inventory_kg'][0])/(4/raw_full['n'])-1)
    comparisons['peak_step_matches']=raw_full['peak_step']==raw_thin['peak_step']
    comparisons['passed']=all(value<=1e-11 for key,value in comparisons.items() if key!='peak_step_matches') and comparisons['peak_step_matches']
    return comparisons


def temporal_comparison(full,half):
    a=np.genfromtxt(Path(full)/'final.csv',delimiter=',',names=True);b=np.genfromtxt(Path(half)/'final.csv',delimiter=',',names=True)
    for key in ('i','j','x_m','y_m','ux_m_s'):np.testing.assert_array_equal(a[key],b[key])
    raw_a=json.loads((Path(full)/'raw.json').read_text());raw_b=json.loads((Path(half)/'raw.json').read_text())
    if abs(raw_a['duration_s']-raw_b['duration_s'])>1e-14:raise ValueError('Time-control duration differs.')
    difference=a['concentration_kg_m3']-b['concentration_kg_m3']
    return dict(final_RMS_difference_over_amplitude=float(np.sqrt(np.mean(difference*difference)))/.25,
                final_maximum_difference_over_amplitude=float(np.max(abs(difference)))/.25,
                RMS_tolerance=.01,time_sensitivity_small=float(np.sqrt(np.mean(difference*difference)))/.25<=.01)


def validate(path):
    p=Path(path);s=json.loads(p.read_text());base=p.parent
    for name,digest in s['source_sha256'].items():
        if sha(ROOT/name)!=digest or sha(base/'source_snapshot'/name)!=digest:raise ValueError('Thin source mismatch: '+name)
    for name,digest in s['artifacts_sha256'].items():
        if sha(base/name)!=digest:raise ValueError('Thin artifact mismatch: '+name)
    parent=ROOT/s['parent_summary']
    if sha(parent)!=s['parent_summary_sha256']:raise ValueError('Thin parent changed.')
    parent_record=json.loads(parent.read_text())
    for name,digest in parent_record['source_sha256'].items():
        if sha(ROOT/name)!=digest:raise ValueError('Parent current source changed: '+name)
    for name,digest in parent_record['artifacts_sha256'].items():
        if sha(parent.parent/name)!=digest:raise ValueError('Parent artifact changed: '+name)
    for key,value in s['cases'].items():
        if analyze(base/'native'/key)!=value:raise ValueError('Thin analysis changed: '+key)
    if equivalence(parent.parent/'native/shear32',base/'native/thin32')!=s['full_cube_equivalence']:raise ValueError('Thin equivalence replay differs.')
    if temporal_comparison(base/'native/thin128',base/'native/thin128_half_dt')!=s['time_comparison']:raise ValueError('Thin temporal comparison differs.')
    return s


def check():
    p=ROOT/json.loads((OUTPUT/'latest.json').read_text())['summary'];s=validate(p)
    result=dict(run_id=s['run_id'],summary=str(p.relative_to(ROOT)),summary_sha256=sha(p),
                sources_checked=len(s['source_sha256']),artifacts_checked=len(s['artifacts_sha256']),
                new_native_runs_in_check=0,new_DWSIM_runs_in_check=0,vendor=verify_vendor(),
                full_cube_equivalence=s['full_cube_equivalence'],time_comparison=s['time_comparison'],cases=s['cases'])
    write_json(OUTPUT/'verification.json',result);print(json.dumps(result,indent=2));return result


def run():
    from ..run import check as check_parent
    parent_verification=check_parent();parent=ROOT/parent_verification['source_summary']
    rid=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid4().hex[:8];base=OUTPUT/'runs'/rid;base.mkdir(parents=True)
    sources=source_manifest([HERE/'run.py',HERE/'protocol.md',*list((HERE/'openlb').iterdir()),HERE.parent/'run.py',HERE.parent/'reference.py'])
    for name in sources:
        dest=base/'source_snapshot'/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/name,dest)
    builddir=ROOT/'.local/adsorption/transport_root_cause_channel_thin'
    build=build_case(HERE/'openlb',builddir,base/'build');jobs=specifications()
    execute(base,jobs[0],builddir/'capillary3d')
    agreement=equivalence(parent.parent/'native/shear32',base/'native/thin32')
    write_json(base/'full_cube_equivalence.json',agreement)
    if not agreement['passed']:raise ValueError('Thin reduction failed full-cube equivalence; no finer runs authorized by protocol.')
    for spec in jobs[1:]:execute(base,spec,builddir/'capillary3d')
    record=dict(schema='pore2process.periodic-shear-thin-diagnosis.v1',run_id=rid,source_sha256=sources,
                parent_summary=str(parent.relative_to(ROOT)),parent_summary_sha256=sha(parent),build=build,
                specifications=jobs,physical_inputs={**PHYSICAL,'end_s':DURATION},actual_target_duration_s=DURATION,
                full_cube_equivalence=agreement,cases={s['id']:analyze(base/'native'/s['id']) for s in jobs},
                time_comparison=temporal_comparison(base/'native/thin128',base/'native/thin128_half_dt'),
                vendor=verify_vendor(),new_native_runs=len(jobs),new_DWSIM_runs=0,application_transport_qualified=False)
    record['artifacts_sha256']={str(p.relative_to(base)):sha(p) for p in base.rglob('*') if p.is_file()}
    write_json(base/'summary.json',record);write_json(OUTPUT/'latest.json',dict(summary=str((base/'summary.json').relative_to(ROOT))))
    print('Thin shear study complete',base,flush=True);return base/'summary.json'


def main():
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');args=p.parse_args()
    check() if args.check else run()

if __name__=='__main__':main()
