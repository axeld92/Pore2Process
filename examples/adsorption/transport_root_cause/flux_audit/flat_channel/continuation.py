"""Finish one missing native case and reassess preserved data with a finer oracle."""
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
from pore2process_demo.runtime.openlb import verify_vendor
from .run import HERE,OUTPUT,PHYSICAL,equivalence,specifications
from .analysis_qualified import analyze

PARENT=OUTPUT/'runs/20261110T212641Z-c0b585ff'


def source_snapshot(files,base):
    sources=source_manifest(files)
    for name in sources:
        target=base/'source_snapshot'/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/name,target)
    return sources


def run():
    # The earlier controller stopped at an unqualified oracle after native64
    # completed. Every parent file is frozen and guarded before continuation.
    parent_artifacts={str(p.relative_to(PARENT)):sha(p) for p in PARENT.rglob('*') if p.is_file()}
    parent_sources={str(p.relative_to(PARENT/'source_snapshot')):sha(p) for p in (PARENT/'source_snapshot').rglob('*') if p.is_file()}
    for name,digest in parent_sources.items():
        if sha(ROOT/name)!=digest:raise ValueError('Original flat-channel source changed: '+name)
    rid=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid4().hex[:8]
    base=OUTPUT/'qualified_runs'/rid;base.mkdir(parents=True)
    files=[HERE/name for name in ('continuation.py','continuation_protocol.md','analysis_qualified.py','reference_qualified.py','run.py','reference.py','protocol.md')]
    files += list((HERE/'openlb').iterdir())+[HERE.parent.parent/'channel/reference.py']
    sources=source_snapshot(files,base)
    interrupted=OUTPUT/'interrupted_attempts.json';shutil.copy2(interrupted,base/interrupted.name)
    write_json(base/'parent_manifest.json',dict(parent=str(PARENT.relative_to(ROOT)),artifacts_sha256=parent_artifacts,source_sha256=parent_sources))
    cases={};directories={}
    for name in ('full32','thin32','constant32','diffusion32','thin64'):
        dest=PARENT/'native'/name
        if not (dest/'execution.json').exists() or json.loads((dest/'execution.json').read_text())['returncode']!=0:
            raise ValueError('A completed parent native case is required: '+name)
        cases[name]=analyze(dest);directories[name]=str(dest.relative_to(ROOT))
    agreement=equivalence(PARENT/'native/full32',PARENT/'native/thin32')
    if not agreement['passed']:raise ValueError('Flat thin/full equivalence failed')
    spec=specifications()[-1];binary=ROOT/'.local/adsorption/transport_root_cause_flat_channel/capillary3d'
    parent_request=json.loads((PARENT/'native/thin64/request.json').read_text())
    if sha(binary)!=parent_request['binary_sha256']:raise ValueError('Native executable changed')
    dest=base/'native/thin96';dest.mkdir(parents=True);p=PHYSICAL
    command=[str(binary),'96','4',str(p['height_m']),str(p['diffusivity_m2_s']),str(p['dt_s']),str(p['steps']),str(p['mean_speed_m_s']),'.25']
    write_json(dest/'request.json',dict(spec=spec,physical_inputs=p,command=command,binary_sha256=sha(binary),parent_native_binary_sha256=parent_request['binary_sha256']))
    print('Flat channel continuation starts thin96',base,flush=True);start=perf_counter()
    with (dest/'stdout.log').open('w') as stream:
        result=subprocess.run(command,cwd=dest,stdout=stream,stderr=subprocess.STDOUT,timeout=3600)
    write_json(dest/'execution.json',dict(returncode=result.returncode,wall_s=perf_counter()-start))
    if result.returncode:raise RuntimeError('Native thin96 failed: '+str(dest/'stdout.log'))
    cases['thin96']=analyze(dest);directories['thin96']=str(dest.relative_to(ROOT))
    for name,digest in parent_artifacts.items():
        if sha(PARENT/name)!=digest:raise ValueError('Parent flat artifact changed: '+name)
    for name,digest in sources.items():
        if sha(ROOT/name)!=digest:raise ValueError('Qualified assessment source changed: '+name)
    record=dict(schema='pore2process.qualified-flat-wall-diagnostic.v1',run_id=rid,source_sha256=sources,
        parent_manifest='parent_manifest.json',parent_manifest_sha256=sha(base/'parent_manifest.json'),
        cases=cases,case_directories=directories,thin_equivalence=agreement,physical_inputs=PHYSICAL,
        parent_completed_native_runs=5,new_native_runs=1,total_completed_native_runs=6,interrupted_native_attempts=2,
        native_attempts_including_interrupted=8,new_DWSIM_runs=0,application_transport_qualified=False,vendor=verify_vendor())
    record['artifacts_sha256']={str(p.relative_to(base)):sha(p) for p in base.rglob('*') if p.is_file()}
    write_json(base/'summary.json',record);write_json(OUTPUT/'qualified_latest.json',dict(summary=str((base/'summary.json').relative_to(ROOT))))
    print('Qualified flat channel complete',base,flush=True);return base/'summary.json'


def check():
    path=ROOT/json.loads((OUTPUT/'qualified_latest.json').read_text())['summary'];s=json.loads(path.read_text());base=path.parent
    for name,digest in s['source_sha256'].items():
        if sha(ROOT/name)!=digest or sha(base/'source_snapshot'/name)!=digest:raise ValueError('Qualified source changed: '+name)
    for name,digest in s['artifacts_sha256'].items():
        if sha(base/name)!=digest:raise ValueError('Qualified artifact changed: '+name)
    manifest=base/s['parent_manifest']
    if sha(manifest)!=s['parent_manifest_sha256']:raise ValueError('Parent manifest changed')
    parent=json.loads(manifest.read_text());directory=ROOT/parent['parent']
    for name,digest in parent['artifacts_sha256'].items():
        if sha(directory/name)!=digest:raise ValueError('Original native artifact changed: '+name)
    for name,digest in parent['source_sha256'].items():
        if sha(ROOT/name)!=digest:raise ValueError('Original native source changed: '+name)
    for name,value in s['cases'].items():
        if analyze(ROOT/s['case_directories'][name])!=value:raise ValueError('Qualified analysis differs: '+name)
    if equivalence(ROOT/s['case_directories']['full32'],ROOT/s['case_directories']['thin32'])!=s['thin_equivalence']:
        raise ValueError('Thin equivalence differs')
    result=dict(source_summary=str(path.relative_to(ROOT)),source_summary_sha256=sha(path),
        sources_checked=len(s['source_sha256']),new_artifacts_checked=len(s['artifacts_sha256']),
        parent_artifacts_checked=len(parent['artifacts_sha256']),total_completed_native_runs=6,
        interrupted_native_attempts=2,new_native_runs_in_check=0,new_DWSIM_runs_in_check=0,vendor=verify_vendor())
    write_json(OUTPUT/'qualified_verification.json',result);print(json.dumps(result,indent=2));return result


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--check',action='store_true');args=parser.parse_args()
    check() if args.check else run()


if __name__=='__main__':main()
