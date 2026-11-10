"""Replay Julia and inspect saved real DWSIM recycle evidence; no new solves."""
import argparse
import csv
import json
import math
from pathlib import Path

from examples.filtration.calibration.check import close_tree
from pore2process_demo.io import sha
from pore2process_demo.bridge import runtime
from examples.filtration.once_through.check import verify as verify_parent
from pore2process_demo.io import read_json as read; from pore2process_demo.paths import ROOT
from examples.filtration.crossflow_recycle.model import translate, adapter_inputs, validate, topology
from examples.filtration.crossflow_recycle.run import check_protocol, jobs, comparison_rows

def verify(path,*,verify_upstream=True):
    path=Path(path).resolve();base=path.parent;s=read(path)
    if s['status']!='passed':raise RuntimeError('Expected a completed study')
    files={str(f.relative_to(base)) for f in base.rglob('*') if f.is_file() and str(f.relative_to(base)) not in ('summary.json','status.json','verification.json')}
    if files!=set(s['artifact_sha256']):raise RuntimeError('Artifact coverage differs')
    for name,digest in s['artifact_sha256'].items():
        f=(base/name).resolve()
        if base not in f.parents or sha(f)!=digest:raise RuntimeError(f'Artifact hash mismatch: {name}')
    for key,value in read(base/'inputs.json').items():close_tree(s[key],value,key)
    p=s['protocol'];parent=read(base/'parent_summary.json');source=read(base/'source_summary.json');check_protocol(p,source)
    if sha(base/'parent_summary.json')!=s['parent_summary_sha256'] or sha(base/'source_summary.json')!=s['source_summary_sha256']:
        raise RuntimeError('Upstream snapshot differs')
    if parent['run_id']!=s['parent_run_id'] or source['run_id']!=s['source_run_id']:raise RuntimeError('Source run ID differs')
    close_tree(read(base/'parent_verification.json'),s['parent_verification'],'parent verification')
    if verify_upstream:
        original=ROOT/s['parent_summary_path']
        if sha(original)!=s['parent_summary_sha256']:raise RuntimeError('Original parent changed')
        close_tree(verify_parent(original),s['parent_verification'],'upstream replay')
    inventory=read(base/'adapter_build/installed_api.json');jl=runtime();expected=list(jobs(p,source))
    if len(s['cases'])!=len(expected):raise RuntimeError('Missing/extra recycle cases')
    for record,request in zip(s['cases'],expected,strict=True):
        close_tree(record['request'],request,'scenario request')
        support=next(v for v in source['supports'] if v['name']==request['support'])
        m,c=translate(jl,source,support,p,request['recovery'],request['module_ratio'],request['resistance_multiplier'])
        close_tree(record['julia_membrane'],m,'Julia membrane');close_tree(record['julia_process'],c,'Julia recycle')
        inputs=adapter_inputs(p,source,c,**request['seed_options']);d=base/'process'/request['case_id']
        actual=read(d/'summary.json');close_tree(record['dwsim'],actual,'saved DWSIM evidence')
        if actual['version']!=inventory['version']:raise RuntimeError('Installed API/solver version differs')
        close_tree(read(d/'inputs.json'),inputs,'adapter request')
        close_tree(record['validation'],validate(actual,inputs,c,p,source),'recomputed stream balances')
        close_tree(record['saved_connections'],topology(d/'crossflow.dwxmz'),'saved topology')
        baseline=next(v for v in s['cases'] if v['request']['condition']==request['condition'] and v['request']['support']=='straight' and not v['request']['validation_only'])
        actual_saving=1-record['validation']['total_shaft_power_W']/baseline['validation']['total_shaft_power_W']
        close_tree(record['energy_saving_vs_straight_fraction'],actual_saving,'actual geometry saving')
        if request['validation_only']:
            for k in ('feed_W','circulation_W','cooling_removed_W'):
                if not math.isclose(actual['powers'][k],baseline['dwsim']['powers'][k],rel_tol=p['acceptance']['seed_invariance_relative_power']):
                    raise RuntimeError('Seed/tolerance invariance failed')
        elif request['support']=='straight':
            nominal=s['cases'][0]
            # Same feed pressure and controlled circulation-inlet properties:
            # W_circulation scales with R*C^2 at constant reference flux.
            ratio=request['resistance_multiplier']*(request['module_ratio']/p['nominal_module_to_permeate_ratio'])**2
            if not math.isclose(actual['powers']['circulation_W']/nominal['dwsim']['powers']['circulation_W'],ratio,rel_tol=1e-8):
                raise RuntimeError('Channel resistance / circulation-square oracle failed')
            feed_ratio=p['nominal_recovery']/request['recovery']
            if not math.isclose(actual['powers']['feed_W']/nominal['dwsim']['powers']['feed_W'],feed_ratio,rel_tol=1e-9):
                raise RuntimeError('Fresh recovery pump oracle failed')
    rejected=read(base/'validation/nonconvergence/rejection.json')
    close_tree(rejected,s['nonconvergence_rejection'],'failed-solve evidence')
    log=(base/'validation/nonconvergence/stdout.log').read_text()
    if rejected['return_code']==0 or rejected['nonconvergence_rejected'] is not True or 'RECYCLE_NOT_CONVERGED' not in log:
        raise RuntimeError('Missing deliberate nonconvergence rejection')
    failed_input=adapter_inputs(p,source,s['cases'][0]['julia_process'],seed_multiplier=.1);failed_input['Maximum_iterations']=1
    close_tree(read(base/'validation/nonconvergence/inputs.json'),failed_input,'deliberate failed request')
    with (base/'comparison.csv').open(newline='') as f:rows=list(csv.DictReader(f))
    expected_rows=comparison_rows(s['cases'])
    if len(rows)!=len(expected_rows):raise RuntimeError('CSV row count differs')
    for row,expected_row in zip(rows,expected_rows,strict=True):
        for k,v in expected_row.items():close_tree(row[k] if isinstance(v,str) else float(row[k]),v,f'CSV.{k}')
    return {'status':'passed','actual_saved_DWSIM_recycle_cases_checked':len(s['cases']),
        'design_scenarios':sum(not v['request']['validation_only'] for v in s['cases']),
        'seed_and_tolerance_variants':sum(v['request']['validation_only'] for v in s['cases']),
        'saved_topologies_checked':len(s['cases']),'nonconvergence_rejection_checked':True,
        'artifacts_checked':len(files),'upstream_reverified':verify_upstream,'new_OpenLB_or_DWSIM_solves_performed_by_checker':False}

def main():
    a=argparse.ArgumentParser(description=__doc__);a.add_argument('summary',nargs='?',type=Path);args=a.parse_args()
    path=args.summary or Path(read(ROOT/'results/crossflow_recycle/latest.json')['summary'])
    print(json.dumps(verify(path),indent=2))

if __name__=='__main__':main()
