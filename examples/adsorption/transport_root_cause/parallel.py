"""Validated stock OpenMP execution of the same instrumented scalar case."""
import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import shutil
from uuid import uuid4

import numpy as np

from pore2process_demo.paths import ROOT
from pore2process_demo.io import sha, write_json
from pore2process_demo.provenance import source_manifest
from pore2process_demo.runtime.openlb import verify_vendor
from examples.adsorption.transport_accuracy.run import source_point, validate_record
from examples.adsorption.bed import analyze as analyze_flow
from .hydraulics.run import native_checks
from .hydraulics.stock_omp import build_case_omp
from .study import HERE, OUTPUT, END, DT, execute, analyze, check


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--grid',type=int,choices=[32,48,64,96],default=32)
    p.add_argument('--point-json',type=Path)
    p.add_argument('--flow-summary',type=Path)
    p.add_argument('--verify-serial',type=Path)
    p.add_argument('--threads',type=int,choices=[2,4],default=4)
    p.add_argument('--constant',action='store_true')
    p.add_argument('--audit-step',type=int,default=-1)
    p.add_argument('--end-s',type=float,default=END)
    p.add_argument('--check',type=Path)
    args=p.parse_args()
    if args.check:
        s=check(args.check)
        for role in ('refined','serial'):
            if role+'_summary' in s:
                parent=ROOT/s[role+'_summary']
                if sha(parent)!=s[role+'_summary_sha256']:raise ValueError('Changed '+role+' summary')
                validate_record(parent)
        return
    os.environ['OMP_NUM_THREADS']=str(args.threads)
    os.environ['OMP_DYNAMIC']='FALSE'
    parent_path=ROOT/json.loads((ROOT/'results/adsorption/transport_accuracy/diagnosis_latest.json').read_text())['summary']
    parent=validate_record(parent_path);hp=ROOT/parent['hydraulic_summary'];validate_record(hp)
    if sha(hp)!=parent['hydraulic_summary_sha256']:raise ValueError('Changed hydraulic parent')
    if args.grid in (32,48):point=source_point(hp.parent,'d800',args.grid)
    else:
        if args.point_json is None or args.flow_summary is None:raise ValueError('New grid requires qualified point and hydraulic summary')
        point=json.loads(args.point_json.read_text());refined=validate_record(args.flow_summary)
        if point['n']!=args.grid:raise ValueError('Source point grid mismatch')
        # Check the actual accepted native field again, rather than trusting a
        # flow declaration or interpolated initialization as a final solution.
        carrier=json.loads((hp.parent/'inputs.json').read_text())
        if not all(native_checks(analyze_flow(Path(point['field']).parent),carrier['native_target_velocity_m_s']).values()):
            raise ValueError('Refined hydraulic field has not passed the physical screens')
        if refined['hydraulic_summary_sha256']!=sha(hp):raise ValueError('Refined point uses a different hydraulic parent')
    if sha(Path(point['field']))!=point['field_sha256']:raise ValueError('Changed input velocity')
    name=f'n{args.grid}'+('_constant' if args.constant else '_sine')
    if args.audit_step>=0:name+=f'_audit{args.audit_step}'
    if args.verify_serial:name+='_parallel_validation'
    spec=dict(point,id=name,n=args.grid,amplitude=0. if args.constant else .25,audit_step=args.audit_step,end_s=args.end_s)
    rid=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid4().hex[:8]
    base=OUTPUT/'parallel/runs'/rid;base.mkdir(parents=True)
    files=[Path(__file__),HERE/'study.py',HERE/'protocol.md',HERE/'hydraulics/run.py',HERE/'hydraulics/stock_omp.py',
           HERE/'flux_audit/population_audit.h',HERE/'flux_audit/uniform_advection_audit.h',
           *list((HERE/'openlb').iterdir()),ROOT/'examples/adsorption/transport_accuracy/profile.py',
           ROOT/'examples/adsorption/transport_accuracy/run.py',ROOT/'examples/adsorption/column_transport/run.py',ROOT/'examples/adsorption/bed.py']
    sources=source_manifest(files)
    for name,digest in sources.items():
        dest=base/'source_snapshot'/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/name,dest)
    generated=base/'generated_source';shutil.copytree(HERE/'openlb',generated)
    for h in ('population_audit.h','uniform_advection_audit.h'):shutil.copy2(HERE/'flux_audit'/h,generated/h)
    build=build_case_omp(generated,ROOT/'.local/adsorption/transport_root_cause_scalar_omp',base/'build')
    binary=Path(build['binary'])
    record=dict(schema='pore2process.transport-root-cause-stock-omp.v1',run_id=rid,source_sha256=sources,
        specifications=[spec],build=build,OMP_NUM_THREADS=args.threads,OMP_DYNAMIC=False,
        parent_summary=str(parent_path.relative_to(ROOT)),parent_summary_sha256=sha(parent_path),
        hydraulic_summary=str(hp.relative_to(ROOT)),hydraulic_summary_sha256=sha(hp),
        new_native_runs=1,new_DWSIM_runs=0,new_adsorption_predictions=0,service_prediction_qualified=False)
    if args.flow_summary:
        record.update(refined_summary=str(args.flow_summary.resolve().relative_to(ROOT)),refined_summary_sha256=sha(args.flow_summary))
        shutil.copy2(args.point_json,base/'source_point.json')
    write_json(base/'plan.json',record);print('Parallel scalar run',base,flush=True)
    execute(base,spec,binary);record['cases']={spec['id']:analyze(base/'native'/spec['id'])}
    if args.verify_serial:
        serial=check(args.verify_serial)
        old_name=f'n{args.grid}_sine'
        old_dest=args.verify_serial.resolve().parent/'native'/old_name
        new_dest=base/'native'/spec['id']
        comparisons={}
        for filename in ('peak.csv','final.csv','inventory.csv','moments.csv'):
            a=np.genfromtxt(old_dest/filename,delimiter=',',names=True)
            b=np.genfromtxt(new_dest/filename,delimiter=',',names=True)
            if a.dtype.names!=b.dtype.names or a.shape!=b.shape:raise ValueError('Parallel diagnostic shape differs')
            for key in a.dtype.names:
                np.testing.assert_allclose(a[key],b[key],rtol=1e-12,atol=2e-12)
            comparisons[filename]=max(float(np.max(abs(a[key]-b[key]))) for key in a.dtype.names)
        record.update(serial_summary=str(args.verify_serial.resolve().relative_to(ROOT)),serial_summary_sha256=sha(args.verify_serial),
                      serial_equivalence_maximum_absolute_differences=comparisons)
    record['vendor']=verify_vendor()
    record['artifacts_sha256']={str(x.relative_to(base)):sha(x) for x in base.rglob('*') if x.is_file()}
    write_json(base/'summary.json',record)
    write_json(OUTPUT/'parallel_latest.json',dict(summary=str((base/'summary.json').relative_to(ROOT))))
    print('Parallel scalar completed',base,flush=True)


if __name__=='__main__':main()
