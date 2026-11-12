"""Execute the two independent diagnostic stages and replay saved evidence."""
import argparse
from datetime import datetime,timezone
import json
import os
from pathlib import Path
import shutil
import subprocess
from time import perf_counter
from uuid import uuid4

import numpy as np

from pore2process_demo.io import sha,write_json
from pore2process_demo.paths import ROOT
from pore2process_demo.provenance import source_manifest
from pore2process_demo.runtime.openlb import verify_vendor
from examples.adsorption.transport_accuracy.run import validate_record
from examples.adsorption.transport_root_cause.hydraulics.stock_omp import build_case_omp
from . import reference as ref
from .analysis import analyze,compare_depth

HERE=Path(__file__).parent
OUTPUT=ROOT/'results/adsorption/transport_verification'
THREADS=4


def cases(stage):
    def spec(name,kind,n,nz=None,speed=ref.U,amplitude=ref.A):
        return dict(name=name,kind=kind,n=n,nz=n if nz is None else nz,speed=speed,
                    amplitude=amplitude,dt=ref.DT,steps=ref.STEPS)
    if stage=='stretching':
        return [spec('uniform32','uniform',32),spec('constant32','stretch',32,amplitude=0.),
                *[spec(f'stretch{n}','stretch',n) for n in (32,64,96)]]
    return [spec('rotation32','disk',32,4),spec('rotation32_full_depth','disk',32),
            *[spec(f'rotation{n}','disk',n,4) for n in (64,96,128)],
            spec('diffusion32','disk',32,4,speed=0.),spec('diffusion128','disk',128,4,speed=0.),
            spec('constant128','disk',128,4,amplitude=0.)]


def snapshots(base):
    files=[HERE/'run.py',HERE/'analysis.py',HERE/'reference.py',HERE/'protocol.md',
           HERE.parent/'transport_root_cause/hydraulics/stock_omp.py',*list((HERE/'openlb').iterdir())]
    sources=source_manifest(files)
    for k in sources:
        p=base/'source_snapshot'/k;p.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/k,p)
    return sources


def execute(base,binary,spec,modes):
    dest=base/'native'/spec['name'];dest.mkdir(parents=True)
    command=[str(binary),*map(str,[spec['kind'],spec['n'],spec['nz'],spec['steps'],ref.L,spec['dt'],ref.D,spec['speed'],spec['amplitude']])]
    write_json(dest/'request.json',dict(spec=spec,command=command,binary_sha256=sha(binary),OMP_NUM_THREADS=THREADS))
    print('Native start',spec['name'],flush=True);start=perf_counter()
    env=dict(os.environ,OMP_NUM_THREADS=str(THREADS),OMP_DYNAMIC='FALSE')
    with (dest/'stdout.log').open('w') as log:
        result=subprocess.run(command,cwd=dest,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=5400)
    write_json(dest/'execution.json',dict(returncode=result.returncode,elapsed_s=perf_counter()-start))
    if result.returncode:raise RuntimeError('Native run failed: '+str(dest/'stdout.log'))
    stats=analyze(dest,modes=modes if spec['kind']=='stretch' and spec['amplitude'] else None)
    write_json(dest/'analysis.json',stats)
    print('Native complete',spec['name'],round(perf_counter()-start,2),'s',stats['range_kg_m3'],stats.get('fields',{}).get('final'),flush=True)
    return stats


def run(stage,parent=None):
    rid=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid4().hex[:8]
    base=OUTPUT/stage/rid;base.mkdir(parents=True);specs=cases(stage)
    s=dict(schema='pore2process.prescribed-transport-verification.v1',stage=stage,run_id=rid,
           source_sha256=snapshots(base),specifications=specs,threads=THREADS,new_DWSIM_runs=0,
           new_bed_runs=0,application_transport_qualified=False)
    write_json(base/'plan.json',s);print('Study',stage,base,flush=True)
    modes=None;binary=base/'binary/capillary3d';binary.parent.mkdir()
    if stage=='stretching':
        print('Qualifying independent stretching reference',flush=True)
        modes,info=ref.qualify_stretch();np.save(base/'reference_modes.npy',modes)
        s['reference']=dict(info,gradient_scales=[ref.gradients(f) for f in modes]);write_json(base/'reference.json',s['reference'])
        print('Reference passed',info,flush=True)
        build=build_case_omp(HERE/'openlb',ROOT/'.local/adsorption/transport_verification',base/'build')
        shutil.copy2(build['binary'],binary);s['build']=build
    else:
        parent=Path(parent).resolve();up=validate_record(parent)
        if up['stage']!='stretching' or not up['reference']['passed'] or len(up['cases'])!=5:
            raise ValueError('Complete independently referenced stretching stage first')
        shutil.copy2(parent.parent/'binary/capillary3d',binary)
        if sha(binary)!=up['binary_sha256']:raise ValueError('Changed archived stock executable')
        s.update(parent_summary=str(parent.relative_to(ROOT)),parent_summary_sha256=sha(parent),
                 reference={'kind':'exact Neumann Bessel mode with rigid rotation','root_J2_prime':ref.ROOT_J2_PRIME})
    s['binary_sha256']=sha(binary)
    write_json(base/'plan.json',s)
    s['cases']={job['name']:execute(base,binary,job,modes) for job in specs}
    if stage=='curved':s['depth_comparison']=compare_depth(base/'native/rotation32',base/'native/rotation32_full_depth')
    s.update(new_native_runs=len(specs),vendor=verify_vendor())
    s['artifacts_sha256']={str(p.relative_to(base)):sha(p) for p in base.rglob('*') if p.is_file()}
    p=base/'summary.json';write_json(p,s);write_json(OUTPUT/f'{stage}_latest.json',{'summary':str(p.relative_to(ROOT))})
    print('Completed',stage,p,flush=True);return p


def check(path):
    path=Path(path).resolve();s=validate_record(path);base=path.parent
    if 'parent_summary' in s:
        parent=ROOT/s['parent_summary']
        if sha(parent)!=s['parent_summary_sha256']:raise ValueError('Changed stretching parent')
        validate_record(parent)
    modes=np.load(base/'reference_modes.npy',allow_pickle=False) if s['stage']=='stretching' else None
    for spec in s['specifications']:
        result=analyze(base/'native'/spec['name'],modes if spec['kind']=='stretch' and spec['amplitude'] else None)
        if result!=s['cases'][spec['name']]:raise ValueError('Analysis changed: '+spec['name'])
    if s['stage']=='curved' and compare_depth(base/'native/rotation32',base/'native/rotation32_full_depth')!=s['depth_comparison']:
        raise ValueError('Depth comparison changed')
    result=dict(summary=str(path.relative_to(ROOT)),summary_sha256=sha(path),source_count=len(s['source_sha256']),
                artifact_count=len(s['artifacts_sha256']),native_runs=s['new_native_runs'],new_native_runs_in_check=0,
                new_DWSIM_runs_in_check=0,vendor=verify_vendor())
    write_json(OUTPUT/f"{s['stage']}_verification.json",result);print(json.dumps(result,indent=2),flush=True)
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--stage',choices=('stretching','curved','both'))
    p.add_argument('--stretching-summary',type=Path);p.add_argument('--check',type=Path);a=p.parse_args()
    if a.check:check(a.check);return
    if a.stage in ('stretching','both'):
        parent=run('stretching')
        if a.stage=='both':run('curved',parent)
    elif a.stage=='curved' and a.stretching_summary:run('curved',a.stretching_summary)
    else:p.error('Choose --stage stretching/both, --stage curved with --stretching-summary, or --check')


if __name__=='__main__':main()
