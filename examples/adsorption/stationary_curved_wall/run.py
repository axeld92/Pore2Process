"""Run or replay the stationary wall benchmark with the independent Julia oracle."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime,timezone
import json
import os
from pathlib import Path
import shutil
import subprocess
from time import perf_counter
from uuid import uuid4

from pore2process_demo.io import sha,write_json
from pore2process_demo.paths import ROOT
from pore2process_demo.provenance import source_manifest
from pore2process_demo.runtime.openlb import verify_vendor
from examples.adsorption.transport_accuracy.run import validate_record
from examples.adsorption.transport_root_cause.hydraulics.stock_omp import build_case_omp
from examples.adsorption.transport_verification.analysis import compare_depth,load
from . import reference as ref
from .analysis import analyze,compare_background,geometry

HERE=Path(__file__).parent
OUTPUT=ROOT/'results/adsorption/stationary_curved_wall'
THREADS=2
WORKERS=2


def specifications():
    return [dict(name=name,n=n,nz=nz,amplitude=amplitude) for name,n,nz,amplitude in [
        ('shear32',32,4,ref.A),('shear128',128,4,ref.A),('constant128',128,4,0.),
        ('shear96',96,4,ref.A),('shear32_full_depth',32,32,ref.A),
        ('shear64',64,4,ref.A),('constant64',64,4,0.)]]


def sources(base):
    files=[HERE/k for k in ('run.py','analysis.py','reference.py','radial_reference.jl','test_reference.jl','protocol.md')]
    files.extend((HERE/'openlb').iterdir())
    files.extend(HERE.parent/k for k in ('transport_verification/analysis.py','transport_verification/reference.py',
                                       'transport_root_cause/hydraulics/stock_omp.py'))
    hashes=source_manifest(files)
    for p in hashes:
        target=base/'source_snapshot'/p;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/p,target)
    return hashes


def execute(base,binary,spec,radial):
    directory=base/'native'/spec['name'];directory.mkdir(parents=True)
    cmd=[str(binary),*map(str,['stationary_disk',spec['n'],spec['nz'],ref.STEPS,ref.L,ref.DT,ref.D,ref.U,spec['amplitude']])]
    write_json(directory/'request.json',dict(command=cmd,spec=spec,binary_sha256=sha(binary),OMP_NUM_THREADS=THREADS))
    print('Native start',spec['name'],flush=True);start=perf_counter()
    env=dict(os.environ,OMP_NUM_THREADS=str(THREADS),OMP_DYNAMIC='FALSE')
    with (directory/'stdout.log').open('w') as stream:
        p=subprocess.run(cmd,cwd=directory,env=env,stdout=stream,stderr=subprocess.STDOUT,timeout=3600)
    write_json(directory/'execution.json',dict(returncode=p.returncode,elapsed_s=perf_counter()-start))
    if p.returncode:raise RuntimeError('OpenLB failed; see '+str(directory/'stdout.log'))
    result=analyze(directory,radial);write_json(directory/'analysis.json',result)
    print('Native complete',spec['name'],round(perf_counter()-start,2),'s',result['range_kg_m3'],result.get('fields',{}).get('final'),flush=True)
    return spec['name'],result


def rigid_comparison(parent,new_cases):
    parent=Path(parent).resolve();old=validate_record(parent)
    if old['stage']!='curved':raise ValueError('Expected archived rotating-wall cohort')
    source=old['cases']['constant128'];target=new_cases['constant128']
    for k in ('n','nz','dt_s','steps','side_m','radius_m','speed_m_s','diffusivity_m2_s','amplitude_kg_m3','tau_odd','tau_even'):
        if source['raw'][k]!=target['raw'][k]:raise ValueError('Rotating control differs in '+k)
    dest=parent.parent/'native/constant128';r,m,u,f,_,_=load(dest)
    _,_,_,link_stats=geometry(dest,m==1,u,f['initial'],r['dx_m'])
    return dict(parent_summary=str(parent.relative_to(ROOT)),parent_summary_sha256=sha(parent),
                rotating_constant_range_kg_m3=source['range_kg_m3'],stationary_constant_range_kg_m3=target['range_kg_m3'],
                rotating_maximum_constant_departure_kg_m3=source['maximum_excursion_kg_m3'],
                stationary_maximum_constant_departure_kg_m3=target['maximum_excursion_kg_m3'],
                rotating_link_velocities=link_stats,stationary_link_velocities={k:target[k] for k in link_stats},
                changed_velocity_role='Equal peak speed, different velocity distribution and wall speed. Compare uniform preservation only.')


def run(parent=None):
    stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid4().hex[:8]
    base=OUTPUT/'runs'/stamp;base.mkdir(parents=True)
    specs=specifications();record=dict(schema='pore2process.stationary-curved-wall.v1',run_id=stamp,source_sha256=sources(base),
        specifications=specs,threads_per_process=THREADS,maximum_concurrent_processes=WORKERS,
        new_bed_runs=0,new_DWSIM_runs=0,application_transport_qualified=False)
    write_json(base/'plan.json',record);print('Stationary wall study',base,flush=True)
    s,f,qualification=ref.qualify(base/'reference');radial=(s,f);record['reference']=qualification
    print('Reference qualified',qualification['radial_points'],qualification['steps_per_window'],qualification['comparisons'],flush=True)
    metadata=build_case_omp(HERE/'openlb',ROOT/'.local/adsorption/stationary_curved_wall',base/'build')
    binary=base/'binary/capillary3d';binary.parent.mkdir();shutil.copy2(metadata['binary'],binary)
    record.update(build=metadata,binary_sha256=sha(binary));write_json(base/'plan.json',record)
    result=dict([execute(base,binary,specs[0],radial)])
    with ThreadPoolExecutor(max_workers=WORKERS) as pool:
        result.update(pool.map(lambda spec:execute(base,binary,spec,radial),specs[1:]))
    record.update(cases=result,new_native_runs=len(result),depth_comparison=compare_depth(base/'native/shear32',base/'native/shear32_full_depth'),
                  background=compare_background(base/'native/shear128',base/'native/constant128',radial),vendor=verify_vendor())
    if parent:record['rigid_comparison']=rigid_comparison(parent,result)
    record['artifacts_sha256']={str(x.relative_to(base)):sha(x) for x in base.rglob('*') if x.is_file()}
    path=base/'summary.json';write_json(path,record);write_json(OUTPUT/'latest.json',dict(summary=str(path.relative_to(ROOT))))
    print('Completed stationary wall study',path,flush=True)
    return path


def check(path):
    path=Path(path).resolve();s=validate_record(path);base=path.parent
    if not s['reference']['passed']:raise ValueError('Radial oracle not qualified')
    mode=base/'reference'/s['reference']['reference_file']
    if sha(mode)!=s['reference']['reference_sha256']:raise ValueError('Changed radial oracle')
    radial=ref.read_mode(mode)
    for job in s['specifications']:
        name=job['name'];d=base/'native'/name
        if analyze(d,radial)!=s['cases'][name]:raise ValueError('Analysis differs: '+name)
        if json.loads((d/'execution.json').read_text())['returncode']!=0:raise ValueError('Native exit failure')
    if compare_depth(base/'native/shear32',base/'native/shear32_full_depth')!=s['depth_comparison']:
        raise ValueError('Depth comparison differs')
    if compare_background(base/'native/shear128',base/'native/constant128',radial)!=s['background']:
        raise ValueError('Background comparison differs')
    if 'rigid_comparison' in s and rigid_comparison(ROOT/s['rigid_comparison']['parent_summary'],s['cases'])!=s['rigid_comparison']:
        raise ValueError('Rotating-wall comparison differs')
    result=dict(summary=str(path.relative_to(ROOT)),summary_sha256=sha(path),sources=len(s['source_sha256']),
        artifacts=len(s['artifacts_sha256']),native_runs=s['new_native_runs'],new_native_runs_in_check=0,vendor=verify_vendor())
    write_json(OUTPUT/'verification.json',result);print(json.dumps(result,indent=2));return result


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--rotation-summary',type=Path);p.add_argument('--check',type=Path);a=p.parse_args()
    check(a.check) if a.check else run(a.rotation_summary)


if __name__=='__main__':main()
