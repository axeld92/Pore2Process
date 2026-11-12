"""Finite staged stock half-way candidate; never changes a native source mid-run."""
from datetime import datetime,timezone
import argparse
import json
from pathlib import Path
import shutil
import subprocess
from time import perf_counter
from uuid import uuid4
from pore2process_demo.paths import ROOT
from pore2process_demo.io import sha,write_json
from pore2process_demo.provenance import source_manifest
from pore2process_demo.runtime.openlb import build_case,verify_vendor
from .analysis import analyze,stage_one_gate
from ..flat_channel.run import PHYSICAL

HERE=Path(__file__).parent
OUTPUT=ROOT/'results/adsorption/transport_root_cause/half_way'
PARENTS={32:ROOT/'results/adsorption/transport_root_cause/runs/20261110T205424Z-ed34b563/native/n32_sine',
         64:ROOT/'results/adsorption/transport_root_cause/parallel/runs/20261110T212155Z-fb52fcc9/native/n64_sine'}


def native(base,name,binary,args,kind,parent=None):
    dest=base/'native'/name;dest.mkdir(parents=True)
    command=[str(binary),*map(str,args)]
    request=dict(command=command,binary_sha256=sha(binary),kind=kind)
    if parent:
        request.update(parent=str(parent.relative_to(ROOT)),parent_request_sha256=sha(parent/'request.json'),
                       parent_summary_raw_sha256=sha(parent/'raw.json'),native_field_sha256=sha(Path(args[6])))
    write_json(dest/'request.json',request)
    print('Half-way starts',name,base,flush=True);start=perf_counter()
    with (dest/'stdout.log').open('w') as stream:
        process=subprocess.Popen(command,cwd=dest,stdout=stream,stderr=subprocess.STDOUT)
        write_json(dest/'process.json',dict(pid=process.pid,command=command))
        code=process.wait(timeout=7200)
    write_json(dest/'execution.json',dict(returncode=code,wall_s=perf_counter()-start))
    if code:raise RuntimeError('Half-way native failed: '+str(dest/'stdout.log'))
    result=analyze(dest,kind);write_json(dest/'analysis.json',result)
    print('Half-way done',name,'bounds',result['maximum_bound_excursion_kg_m3'],'RMS',result.get('RMS_error_over_amplitude'),flush=True)
    return result


def run():
    rid=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid4().hex[:8]
    base=OUTPUT/'runs'/rid;base.mkdir(parents=True)
    files=[HERE/'run.py',HERE/'analysis.py',HERE/'protocol.md',HERE/'halfway_setup.h',
           *list((HERE/'flat_openlb').iterdir()),*list((HERE/'porous_openlb').iterdir()),
           HERE.parent/'flat_channel/run.py',HERE.parent/'flat_channel/reference.py',
           HERE.parent/'flat_channel/reference_qualified.py',HERE.parent.parent/'channel/reference.py']
    sources=source_manifest(files)
    for name in sources:
        dest=base/'source_snapshot'/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/name,dest)
    p=PHYSICAL;cases={};builds={};builddir=ROOT/'.local/adsorption/transport_root_cause_halfway_flat'
    builds['flat']=build_case(HERE/'flat_openlb',builddir,base/'flat_build')
    for name,n,speed,amplitude in [('constant32',32,1.,0.),('diffusion32',32,0.,.25),('shear32',32,1.,.25),('shear96',96,1.,.25)]:
        args=[n,4,p['height_m'],p['diffusivity_m2_s'],p['dt_s'],p['steps'],p['mean_speed_m_s']*speed,amplitude]
        cases[name]=native(base,name,builddir/'capillary3d',args,'flat')
    flat_passed=stage_one_gate(cases)
    write_json(base/'stage_one.json',dict(passed=flat_passed,cases=cases))
    if flat_passed:
        builddir=ROOT/'.local/adsorption/transport_root_cause_halfway_porous'
        builds['porous']=build_case(HERE/'porous_openlb',builddir,base/'porous_build')
        for n in (32,64):
            parent=PARENTS[n];request=json.loads((parent/'request.json').read_text());args=request['command'][1:]
            if sha(Path(args[6]))!=request['native_field_sha256']:raise ValueError('Archived flow field changed')
            cases[f'porous{n}']=native(base,f'porous{n}',builddir/'capillary3d',args,'porous',parent)
    for name,digest in sources.items():
        if sha(ROOT/name)!=digest:raise ValueError('Half-way source changed during run: '+name)
    record=dict(schema='pore2process.stock-halfway-candidate.v1',run_id=rid,source_sha256=sources,
                cases=cases,builds=builds,stage_one_passed=flat_passed,new_native_runs=len(cases),
                new_DWSIM_runs=0,application_transport_qualified=False,vendor=verify_vendor())
    record['artifacts_sha256']={str(p.relative_to(base)):sha(p) for p in base.rglob('*') if p.is_file()}
    write_json(base/'summary.json',record);write_json(OUTPUT/'latest.json',dict(summary=str((base/'summary.json').relative_to(ROOT))))
    print('Half-way complete',base,flush=True);return base/'summary.json'


def check():
    path=ROOT/json.loads((OUTPUT/'latest.json').read_text())['summary'];s=json.loads(path.read_text());base=path.parent
    for name,digest in s['source_sha256'].items():
        if sha(ROOT/name)!=digest or sha(base/'source_snapshot'/name)!=digest:raise ValueError('Changed source: '+name)
    for name,digest in s['artifacts_sha256'].items():
        if sha(base/name)!=digest:raise ValueError('Changed artifact: '+name)
    for name,value in s['cases'].items():
        if analyze(base/'native'/name,'porous' if name.startswith('porous') else 'flat')!=value:raise ValueError('Changed analysis: '+name)
        request=json.loads((base/'native'/name/'request.json').read_text())
        if 'parent' in request:
            parent=ROOT/request['parent']
            if sha(parent/'request.json')!=request['parent_request_sha256'] or sha(parent/'raw.json')!=request['parent_summary_raw_sha256']:raise ValueError('Changed parent')
            if sha(Path(request['command'][7]))!=request['native_field_sha256']:raise ValueError('Changed native field')
    result=dict(source_summary=str(path.relative_to(ROOT)),source_summary_sha256=sha(path),
                sources_checked=len(s['source_sha256']),artifacts_checked=len(s['artifacts_sha256']),
                completed_native_runs=s['new_native_runs'],new_native_runs_in_check=0,vendor=verify_vendor())
    write_json(OUTPUT/'verification.json',result);print(json.dumps(result,indent=2));return result


def main():
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();check() if a.check else run()


if __name__=='__main__':main()
