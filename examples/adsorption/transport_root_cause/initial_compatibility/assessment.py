"""Typed assessment continuation preserving the original completed native run.

The first controller completed OpenLB but strict JSON rejected NumPy booleans.
Only serialization types change here; all initial-data and simulation sources,
thresholds, conditions and original native outputs remain untouched.
"""
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path
import shutil
import subprocess
from time import perf_counter
from uuid import uuid4

import numpy as np

from pore2process_demo.paths import ROOT
from pore2process_demo.io import sha,write_json
from pore2process_demo.provenance import source_manifest,verified_source_path
from pore2process_demo.runtime.openlb import verify_vendor
from examples.adsorption.transport_accuracy.run import validate_record
from examples.adsorption.transport_root_cause.study import DT,END
from .run import HERE,OUTPUT,analyze


def ordinary(value):
    if isinstance(value,np.generic):return value.item()
    if isinstance(value,dict):return {key:ordinary(item) for key,item in value.items()}
    if isinstance(value,(list,tuple)):return [ordinary(item) for item in value]
    return value


def preserve_controller_failure(parent):
    target=parent/'controller_failure.json'
    if target.exists():return target,validate_record(target)
    if (parent/'summary.json').exists():raise ValueError('Do not alter a completed controller summary.')
    plan=json.loads((parent/'plan.json').read_text())
    for name,digest in plan['source_sha256'].items():
        if sha(verified_source_path(name,digest))!=digest:raise ValueError('Changed original control source.')
    execution=json.loads((parent/'native/execution.json').read_text())
    if execution['returncode']!=0:raise ValueError('A successful original native execution is required.')
    record=dict(plan,status='native_completed_controller_serialization_failed',
        native_execution=execution,controller_error='TypeError: NumPy bool is not JSON serializable',
        native_outputs_unchanged=True,analysis_saved=False,
        failure_record_writer=str(Path(__file__).relative_to(ROOT)),failure_record_writer_sha256=sha(Path(__file__)))
    record['artifacts_sha256']={str(p.relative_to(parent)):sha(p) for p in parent.rglob('*') if p.is_file()}
    write_json(target,record)
    return target,record


def check(summary):
    summary=Path(summary).resolve();s=validate_record(summary)
    for role in ('parent','hydraulic','refined'):
        if not s.get(role+'_summary'):continue
        path=ROOT/s[role+'_summary']
        if sha(path)!=s[role+'_summary_sha256']:raise ValueError('Changed '+role+' parent.')
        validate_record(path)
    native=ROOT/s['native_directory']
    if ordinary(analyze(native))!=s['analysis']:raise ValueError('Typed assessment differs on replay.')
    print('Checked compatibleinitial assessment:',len(s['source_sha256']),'sources',len(s['artifacts_sha256']),
        'artifacts; 0newnative/DWSIM',flush=True)
    return s


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--native-parent',type=Path)
    parser.add_argument('--grid',type=int,choices=[32,64],default=32)
    parser.add_argument('--flow-summary',type=Path)
    parser.add_argument('--check',type=Path)
    args=parser.parse_args()
    if args.check:check(args.check/'summary.json');return
    if args.native_parent is None:raise ValueError('Name the preserved original native control.')
    parent=args.native_parent.resolve()
    parent_path,parent_record=preserve_controller_failure(parent)
    validate_record(parent_path)
    hydraulic=ROOT/parent_record['hydraulic_summary'];validate_record(hydraulic)
    if sha(hydraulic)!=parent_record['hydraulic_summary_sha256']:raise ValueError('Changed original hydraulic parent.')
    refined=None
    if args.grid==32:
        point=parent_record['specification']
    else:
        if args.flow_summary is None:raise ValueError('N64 needs its actual qualified refined flow summary.')
        refined_path=args.flow_summary.resolve();refined=validate_record(refined_path)
        if not refined['checks_passed'] or refined['resolution']!=64 or refined['hydraulic_summary_sha256']!=sha(hydraulic):
            raise ValueError('Wrong or unqualified refined operating point.')
        point=refined['source_point']
    rid=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid4().hex[:8]
    base=OUTPUT/'assessments'/rid;base.mkdir(parents=True)
    files=[Path(__file__),HERE/'reference.py',HERE/'run.py',HERE/'protocol.md',*list((HERE/'openlb').iterdir()),
        ROOT/'examples/adsorption/transport_root_cause/study.py',ROOT/'examples/adsorption/transport_accuracy/profile.py',
        ROOT/'examples/adsorption/transport_accuracy/run.py',ROOT/'examples/adsorption/transport_accuracy/analysis.py',
        ROOT/'examples/adsorption/column_transport/run.py']
    sources=source_manifest(files)
    for name in sources:
        dest=base/'source_snapshot'/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/name,dest)
    record=dict(schema='pore2process.initial-compatibility-typed-assessment.v1',run_id=rid,
        source_sha256=sources,parent_summary=str(parent_path.relative_to(ROOT)),parent_summary_sha256=sha(parent_path),
        hydraulic_summary=str(hydraulic.relative_to(ROOT)),hydraulic_summary_sha256=sha(hydraulic),
        refined_summary=None if refined is None else str(refined_path.relative_to(ROOT)),
        refined_summary_sha256=None if refined is None else sha(refined_path),
        original_profile_qualified=False,new_adsorption_predictions=0,new_DWSIM_runs=0,
        new_native_runs=0 if args.grid==32 else 1,reused_native_runs=1 if args.grid==32 else 0,
        serialization_change='NumPy scalar values converted to ordinary JSON scalar types; no threshold/model change.',
        specification=dict(point,n=args.grid))
    if args.grid==32:
        native=parent/'native'
    else:
        source_binary=parent/'binary/capillary3d'
        if sha(source_binary)!=json.loads((parent/'native/request.json').read_text())['binary_sha256']:
            raise ValueError('Original initial-control binary changed.')
        binary=base/'binary/capillary3d';binary.parent.mkdir();shutil.copy2(source_binary,binary)
        native=base/'native';native.mkdir()
        command=[str(binary),'64','1',repr(point['side_m']),'1e-9',
            repr(.5+3e-9*DT/(point['side_m']/64)**2),repr(END),point['field'],
            repr(point['speed_m_s']),'.5',repr(point['velocity_scale']),'1','zero','-1','.25']
        write_json(native/'request.json',dict(command=command,binary_sha256=sha(binary),
            native_field_sha256=sha(Path(point['field']))))
        record['binary_reused_from_parent']=str(source_binary.relative_to(ROOT))
        record['binary_sha256']=sha(binary)
        write_json(base/'plan.json',record)
        print('N64compatibleinitial control',base,flush=True)
        started=perf_counter()
        with (native/'stdout.log').open('w') as stream:
            result=subprocess.run(command,cwd=native,stdout=stream,stderr=subprocess.STDOUT,timeout=7200)
        write_json(native/'execution.json',dict(returncode=result.returncode,elapsed_s=perf_counter()-started))
        if result.returncode:raise RuntimeError('N64nativefailed: '+str(native/'stdout.log'))
    record['native_directory']=str(native.relative_to(ROOT))
    record['analysis']=ordinary(analyze(native))
    if not record['analysis']['initial_compatibility']['passed']:
        raise ValueError('Do not interpret an unverified compatibleinitial construction.')
    record['vendor']=verify_vendor()
    record['artifacts_sha256']={str(p.relative_to(base)):sha(p) for p in base.rglob('*') if p.is_file()}
    write_json(base/'summary.json',record)
    write_json(OUTPUT/f'grid{args.grid}_assessed_latest.json',dict(summary=str((base/'summary.json').relative_to(ROOT))))
    print('Typed compatibleinitial assessment',base,'range',record['analysis']['minimum_concentration_kg_m3'],
        record['analysis']['maximum_concentration_kg_m3'],flush=True)


if __name__=='__main__':main()
