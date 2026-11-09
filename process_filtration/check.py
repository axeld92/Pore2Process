"""Hash-check saved actual DWSIM evidence and replay Julia; no external solves."""
import argparse
import csv
import json
import math
from pathlib import Path

from calibration.check import close_tree
from coupling.run import sha
from graded_support.check import verify as verify_source
from pore2process_demo.bridge import runtime
from process_filtration.model import translate,adapter_inputs,validate,topology
from process_filtration.run import ROOT,read,check_protocol,jobs,comparison_rows

def verify(path,*,verify_upstream=True):
    path=Path(path).resolve();base=path.parent;s=read(path)
    if s['status']!='passed': raise RuntimeError('Study did not pass')
    files={str(f.relative_to(base)) for f in base.rglob('*') if f.is_file() and str(f.relative_to(base)) not in ('summary.json','status.json','verification.json')}
    if files!=set(s['artifact_sha256']): raise RuntimeError('Artifact coverage differs')
    for name,digest in s['artifact_sha256'].items():
        f=(base/name).resolve()
        if base not in f.parents or sha(f)!=digest: raise RuntimeError(f'Artifact hash mismatch: {name}')
    inputs=read(base/'inputs.json')
    for key,value in inputs.items(): close_tree(s[key],value,key)
    source=read(base/'source_summary.json');p=s['protocol'];check_protocol(p,source)
    if sha(base/'source_summary.json')!=s['source_summary_sha256'] or source['run_id']!=s['source_run_id']:
        raise RuntimeError('Source boundary differs')
    inventory=read(base/'adapter_build/installed_api.json')
    if any(c['dwsim']['version']!=inventory['version'] for c in s['cases']):
        raise RuntimeError('Installed API/solver versions differ')
    if verify_upstream:
        source_path=ROOT/s['source_summary_path']
        if sha(source_path)!=s['source_summary_sha256']: raise RuntimeError('Upstream source summary changed')
        close_tree(verify_source(source_path,write_record=False),s['source_verification'],'upstream verification')
    close_tree(read(base/'source_verification.json'),s['source_verification'],'source evidence')
    jl=runtime();expected_jobs=list(jobs(p,source))
    if len(s['cases'])!=len(expected_jobs): raise RuntimeError('Missing/extra process cases')
    for c,(name,support,recovery,channel) in zip(s['cases'],expected_jobs,strict=True):
        for key,value in (('case_id',name),('support',support['name']),('recovery_fraction',recovery),('channel_pressure_drop_Pa',channel)):
            close_tree(c[key],value,key)
        load,closure=translate(jl,source,support,recovery,channel)
        close_tree(c['julia_membrane'],load,'membrane closure');close_tree(c['julia_process'],closure,'process closure')
        directory=base/'process'/name;actual=read(directory/'summary.json');snapshot=adapter_inputs(p,source,closure)
        close_tree(read(directory/'inputs.json'),snapshot,'adapter request')
        close_tree(c['dwsim'],actual,'actual saved DWSIM evidence')
        close_tree(c['validation'],validate(actual,snapshot,closure,p,source),'stream validation')
        close_tree(c['saved_connections'],topology(directory/'filtration.dwxmz'),'actual saved topology')
        group=[v for v in s['cases'] if v['recovery_fraction']==recovery and v['channel_pressure_drop_Pa']==channel]
        b=group[0];saving=1-actual['shaft_power_W']/b['dwsim']['shaft_power_W']
        # Independent reference: same recovery/efficiency; only TMP changes.
        oracle=(b['julia_membrane']['pressure_drop_Pa']-load['pressure_drop_Pa'])/(b['julia_membrane']['pressure_drop_Pa']+channel)
        if not math.isclose(saving,oracle,rel_tol=1e-10,abs_tol=1e-14): raise RuntimeError('Geometry savings oracle failed')
        close_tree(c['energy_saving_vs_straight_fraction'],saving,'geometry saving')
        if recovery!=p['nominal_recovery']:
            other=next(v for v in s['cases'] if v['recovery_fraction']==p['nominal_recovery'] and v['channel_pressure_drop_Pa']==channel and v['support']==support['name'])
            ratio=actual['shaft_power_W']/other['dwsim']['shaft_power_W']
            if not math.isclose(ratio,p['nominal_recovery']/recovery,rel_tol=1e-9): raise RuntimeError('Recovery sweep oracle failed')
    with (base/'comparison.csv').open(newline='') as f: rows=list(csv.DictReader(f))
    expected=comparison_rows(s['cases'])
    if len(rows)!=len(expected): raise RuntimeError('CSV row count differs')
    for row,e in zip(rows,expected,strict=True):
        for k,v in e.items(): close_tree(row[k] if isinstance(v,str) else float(row[k]),v,f'CSV.{k}')
    return {'status':'passed','actual_saved_DWSIM_filtration_cases_checked':len(expected),
        'actual_saved_topologies_checked':len(expected),'Julia_process_loads_replayed':len(expected),
        'artifacts_checked':len(files),'upstream_reverified':verify_upstream,
        'new_OpenLB_or_DWSIM_solves_performed_by_checker':False}

def main():
    a=argparse.ArgumentParser(description=__doc__);a.add_argument('summary',nargs='?',type=Path)
    args=a.parse_args();path=args.summary or Path(read(ROOT/'results/process_filtration/latest.json')['summary'])
    print(json.dumps(verify(path),indent=2))

if __name__=='__main__':main()
