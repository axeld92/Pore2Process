"""Execute a separate normal-compatible initial-data diagnostic."""
import argparse
from datetime import datetime,timezone
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
from examples.adsorption.transport_accuracy.run import validate_record,source_point
from examples.adsorption.transport_root_cause.study import DT,END,analyze as original_analysis
from .reference import validate_initial

HERE=Path(__file__).parent
OUTPUT=ROOT/'results/adsorption/transport_root_cause/initial_compatibility'


def analyze(directory):
    a=original_analysis(directory)
    a['initial_compatibility']=validate_initial(directory)
    a['passed']=a['passed'] and a['initial_compatibility']['passed']
    a['original_profile_qualified']=False
    return a


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--grid',type=int,choices=[32,64],default=32)
    p.add_argument('--flow-summary',type=Path)
    p.add_argument('--check',type=Path)
    args=p.parse_args()
    if args.check:
        summary=args.check.resolve()/'summary.json';s=validate_record(summary)
        if sha(ROOT/s['hydraulic_summary'])!=s['hydraulic_summary_sha256']:raise ValueError('Hydraulic parent changed.')
        validate_record(ROOT/s['hydraulic_summary'])
        if s.get('flow_summary'):
            if sha(ROOT/s['flow_summary'])!=s['flow_summary_sha256']:raise ValueError('Refined flow parent changed.')
            validate_record(ROOT/s['flow_summary'])
        if analyze(summary.parent/'native')!=s['analysis']:raise ValueError('Initial control replay differs.')
        print('Checked normal-compatible control:',len(s['source_sha256']),'sources',len(s['artifacts_sha256']),'artifacts; 0newnative/DWSIM',flush=True)
        return
    hydraulic_path=ROOT/json.loads((ROOT/'results/adsorption/column_transport/hydraulics/latest.json').read_text())['summary']
    parent=validate_record(hydraulic_path)
    if not parent['hydraulic_checks_passed']:raise ValueError('Qualified actual hydraulic parent required.')
    flow=None
    if args.grid==32:
        point=source_point(hydraulic_path.parent,'d800',32)
    else:
        if args.flow_summary is None:raise ValueError('N64 requires its qualified flow summary.')
        flow_path=args.flow_summary.resolve();flow=validate_record(flow_path)
        if not flow['checks_passed'] or flow['resolution']!=64 or flow['hydraulic_summary_sha256']!=sha(hydraulic_path):
            raise ValueError('Wrong or unqualified N64 physical flow.')
        point=flow['source_point']
    rid=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid4().hex[:8]
    base=OUTPUT/'runs'/rid;base.mkdir(parents=True)
    files=[Path(__file__),HERE/'reference.py',HERE/'protocol.md',*list((HERE/'openlb').iterdir()),
        ROOT/'examples/adsorption/transport_root_cause/study.py',
        ROOT/'examples/adsorption/transport_accuracy/profile.py',ROOT/'examples/adsorption/transport_accuracy/run.py',
        ROOT/'examples/adsorption/transport_accuracy/analysis.py',ROOT/'examples/adsorption/column_transport/run.py']
    sources=source_manifest(files)
    for name in sources:
        dest=base/'source_snapshot'/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/name,dest)
    metadata=build_case(HERE/'openlb',ROOT/'.local/adsorption/initial_compatibility',base/'build')
    binary=base/'binary/capillary3d';binary.parent.mkdir();shutil.copy2(ROOT/'.local/adsorption/initial_compatibility/capillary3d',binary)
    spec=dict(point,n=args.grid)
    record=dict(schema='pore2process.initial-normal-compatibility.v1',run_id=rid,
        source_sha256=sources,specification=spec,build=metadata,
        hydraulic_summary=str(hydraulic_path.relative_to(ROOT)),hydraulic_summary_sha256=sha(hydraulic_path),
        flow_summary=None if flow is None else str(flow_path.relative_to(ROOT)),flow_summary_sha256=None if flow is None else sha(flow_path),
        new_native_runs=1,new_DWSIM_runs=0,new_adsorption_predictions=0,original_profile_qualified=False,
        changed_input='Initial concentration only; surface-trace preserving normal extension and stock quintic blend.')
    write_json(base/'plan.json',record)
    dest=base/'native';dest.mkdir()
    command=[str(binary),str(args.grid),'1',repr(point['side_m']),'1e-9',
        repr(.5+3e-9*DT/(point['side_m']/args.grid)**2),repr(END),point['field'],
        repr(point['speed_m_s']),'.5',repr(point['velocity_scale']),'1','zero','-1','.25']
    write_json(dest/'request.json',dict(command=command,binary_sha256=sha(binary),native_field_sha256=sha(Path(point['field']))))
    print('Normal-compatible initial control',base,flush=True)
    start=perf_counter()
    with (dest/'stdout.log').open('w') as stream:
        result=subprocess.run(command,cwd=dest,stdout=stream,stderr=subprocess.STDOUT,timeout=7200)
    write_json(dest/'execution.json',dict(returncode=result.returncode,elapsed_s=perf_counter()-start))
    if result.returncode:raise RuntimeError('Native initial control failed: '+str(dest/'stdout.log'))
    record['analysis']=analyze(dest);record['vendor']=verify_vendor()
    record['artifacts_sha256']={str(x.relative_to(base)):sha(x) for x in base.rglob('*') if x.is_file()}
    write_json(base/'summary.json',record)
    write_json(OUTPUT/f'grid{args.grid}_latest.json',dict(summary=str((base/'summary.json').relative_to(ROOT))))
    print('Initial compatibility checks',record['analysis']['initial_compatibility']['checks'],flush=True)
    print('Initial control range',record['analysis']['minimum_concentration_kg_m3'],record['analysis']['maximum_concentration_kg_m3'],flush=True)


if __name__=='__main__':main()
