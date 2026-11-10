"""Replay saved OpenLB fields, Julia contracts and actual DWSIM archives."""
import argparse
import csv
import math
import numpy as np
from pathlib import Path
from calibration.check import close_tree
from coupling.run import sha
from hydraulic.run import write_json
from hydraulic.build import ROOT,verify_vendor
from pore2process_demo.bridge import runtime
from process_filtration.run import read
from crossflow_recycle.model import validate,topology
from star_channel.geometry import polygon
from .native import inspect,checks
from .model import translate,inputs,capacity_decisions,comparison,requests,source_shape_diagnostics
from .selection import selections


def artifacts(path,s):
    base=path.parent.resolve()
    found={str(f.relative_to(base)) for f in base.rglob('*') if f.is_file() and f.resolve()!=path.resolve() and f.name not in ('verification.json','status.json')}
    if found!=set(s['artifacts_sha256']):raise ValueError('Changed artifact coverage.')
    for name,digest in s['artifacts_sha256'].items():
        f=(base/name).resolve()
        if base not in f.parents or sha(f)!=digest:raise ValueError('Changed artifact: '+name)
    return len(s['artifacts_sha256'])


def verify_native(path):
    path=Path(path).resolve();s=read(path);base=path.parent
    if s['status']!='passed':raise ValueError('Native envelope not accepted.')
    count=artifacts(path,s)
    for name,value in read(base/'inputs.json').items():close_tree(s[name],value,'Native input snapshot.'+name)
    for name,digest in s['source_sha256'].items():
        if sha(ROOT/name)!=digest:raise ValueError('Native application changed; rerun verification.')
    replay={name:inspect(base/name) for name in s['cases']}
    close_tree(replay,s['cases'],'OpenLB field replay')
    jl=runtime();v=jl.duct_stokes_reference(jl.seval('Matrix{Float64}')(polygon()),256)
    close_tree(float(v.conductance_geometry_m4),s['reference_C_m4'],'Independent Julia reference')
    close_tree(checks(replay,s['reference_C_m4'],s['protocol']),s['checks'],'Native numerical checks')
    if s['geometry_coefficient_m4']!=s['cases'][s['accepted_case']]['geometry_coefficient_m4']:
        raise ValueError('Process coefficient does not come from the accepted native case.')
    expected=[s['cases']['high1000_n64']['Re_mean'],s['cases']['high2000_n64']['Re_mean']]
    close_tree(expected,s['verified_Re_range'],'Native Reynolds envelope')
    seed=np.genfromtxt(base/'seed/section.csv',delimiter=',',names=True)
    nominal=np.genfromtxt(base/'high1388_n32/section.csv',delimiter=',',names=True)
    if seed.shape!=nominal.shape or not np.array_equal(seed['j'],nominal['j']) or not np.array_equal(seed['k'],nominal['k']):raise ValueError('Seed field indexing changed.')
    error=float(np.linalg.norm(seed['ux_m_s']-nominal['ux_m_s'])/np.linalg.norm(nominal['ux_m_s']))
    if error>1e-5:raise ValueError('Initial-estimate velocity profile invariance failed.')
    return {'status':'passed','native_cases':len(replay),'artifacts_checked':count,
            'seed_profile_relative_L2':error,'original_vendor':verify_vendor()}


def verify(path,upstream=True):
    path=Path(path).resolve();base=path.parent;s=read(path)
    if s['status']!='passed':raise ValueError('Incomplete map.')
    count=artifacts(path,s)
    for name,digest in s['implementation_sha256'].items():
        if sha(ROOT/name)!=digest:raise ValueError('Map implementation changed: '+name)
    for name,original,digest in [('native_summary.json','native_path','native_sha256'),('support_summary.json','support_path','support_sha256'),('parent_summary.json','parent_path','parent_sha256')]:
        if sha(base/name)!=s[digest]:raise ValueError('Changed source snapshot.')
        if upstream and sha(Path(s[original]))!=s[digest]:raise ValueError('Changed upstream summary.')
    native,source=read(base/'native_summary.json'),read(base/'support_summary.json')
    if upstream:
        close_tree(verify_native(s['native_path']),s['native_verification'],'Native ancestry')
        from graded_support.check import verify as verify_support
        close_tree(verify_support(Path(s['support_path']),write_record=False),s['support_verification'],'TRT support ancestry')
    from .evidence import extract
    current_evidence=extract()
    close_tree({k:v for k,v in current_evidence.items() if k!='source_audit'},
               {k:v for k,v in s['evidence'].items() if k!='source_audit'},'Literature measurements/assumptions')
    for key,value in read(base/'inputs.json').items():close_tree(value,s[key],'Snapshot.'+key)
    jl=runtime();covered=[]
    close_tree(source_shape_diagnostics(jl,s['evidence']),s['source_shape_diagnostics'],'Source interpolation-shape diagnostic')
    for index,c in enumerate(s['cases']+[s['seed_validation']]):
        seeded=index==len(s['cases']);q=c['request'];directory=base/('validation/seed' if seeded else 'process/'+c['id'])
        t=translate(jl,native,source,s['protocol'],q);close_tree(t,c['translation'],'Julia replay')
        close_tree(capacity_decisions(jl,t,s['evidence'],s['protocol']),c['capacity'],'Source capacity replay')
        kw={'seed_multiplier':.5,'tolerance_multiplier':.1,'seed_temperature_offset':-1.} if seeded else {}
        a=inputs(s['hydraulics'],source,t,**kw)
        close_tree(a,read(directory/'inputs.json'),'Actual simulator request')
        actual=read(directory/'summary.json');close_tree(actual,c['dwsim'],'Actual DWSIM output')
        close_tree(validate(actual,a,t['process'],s['hydraulics'],source),c['validation'],'Mass/enthalpy/pressure/reload checks')
        close_tree(topology(directory/'crossflow.dwxmz'),c['saved_connections'],'Actual saved topology')
        if t['sustainable_flux_LMH'] is not None or t['fouling_time_s'] is not None:raise ValueError('Unavailable physical predictions replaced.')
        with (directory/'axial_profile.csv').open(newline='') as f:csvrows=list(csv.DictReader(f))
        a=t['axial'];f=t['feed'];expected=list(zip(a['x_m'],a['pressure_Pa'],a['local_flux_m_s'],a['channel_flow_m3_s'],f['axial_COD_kg_m3'],strict=True))
        if len(csvrows)!=101:raise ValueError('Axial shape changed.')
        for row,values in zip(csvrows,expected,strict=True):
            close_tree([float(row[k]) for k in ('x_m','TMP_Pa','local_flux_m_s','channel_flow_positive_x_m3_s','COD_kg_m3')],list(values),'SI axial export')
        if not seeded:covered.append(q)
    for rejected in s['rejected_scenarios']:
        try:translate(jl,native,source,s['protocol'],rejected['request'])
        except ValueError as exc:
            if str(exc)!=rejected['reason']:raise ValueError('Hydraulic rejection changed.')
        else:raise ValueError('Rejected scenario accepted on replay.')
        covered.append(rejected['request'])
    expected=list(requests(s['protocol']))
    if sorted(map(str,covered))!=sorted(map(str,expected)):raise ValueError('Missing or duplicate design points.')
    close_tree(comparison(s['cases']),s['comparison'],'Energy and material comparison')
    close_tree(selections(s['cases'],s['protocol']),s['selections'],'Qualified discrete selections')
    for filename,name in [('comparison.csv','comparison'),('selections.csv','selections')]:
        with (base/filename).open(newline='') as f:rows=list(csv.DictReader(f))
        if len(rows)!=len(s[name]):raise ValueError('Changed table shape.')
        for row,expected in zip(rows,s[name],strict=True):
            for key,value in expected.items():
                if value is None:
                    if row[key]!='':raise ValueError('Unavailable table value changed.')
                elif isinstance(value,str):close_tree(row[key],value,'Table string')
                else:close_tree(float(row[key]),value,'Table number')
    nominal=next(c for c in s['cases'] if c['request']==s['seed_validation']['request'])
    for k in ('feed_W','circulation_W','cooling_removed_W'):
        if not math.isclose(nominal['dwsim']['powers'][k],s['seed_validation']['dwsim']['powers'][k],rel_tol=1e-8):raise ValueError('Seed invariance failed.')
    if not s['nonconvergence']['nonconvergence_rejected'] or 'RECYCLE_NOT_CONVERGED' not in (base/'validation/nonconvergence/stdout.log').read_text():
        raise ValueError('Failed recycle was not rejected.')
    return {'status':'passed','actual_DWSIM_solved_saved_reloaded_cases':len(s['cases'])+1,
        'rejected_hydraulic_scenarios':len(s['rejected_scenarios']),'capacity_decisions':sum(len(c['capacity']) for c in s['cases']),
        'discrete_policy_selections':len(s['selections']),'artifacts_checked':count,'upstream_rechecked':upstream,'external_solves_by_checker':False}

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('summary',nargs='?',type=Path);a=parser.parse_args()
    path=a.summary or Path(read(ROOT/'results/operating_map/process/latest.json')['summary'])
    import json
    result=verify(path);write_json(path.parent/'verification.json',result)
    print(json.dumps(result,indent=2))
