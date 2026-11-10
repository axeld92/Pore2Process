"""Replay native fields, Julia scale translation, saved DWSIM ports and balances."""
from pore2process_demo.provenance import verified_source_path
import argparse
import csv
import math
from pathlib import Path

from examples.filtration.calibration.check import close_tree
from pore2process_demo.io import sha
from examples.filtration.crossflow_recycle.model import validate, topology
from pore2process_demo.runtime.openlb import verify_vendor
from pore2process_demo.bridge import runtime
from pore2process_demo.paths import ROOT; from pore2process_demo.io import read_json as read
from .native import analyze, verify_suite
from .model import translate, inputs, requests

def verify_native(path):
    path=Path(path);s=read(path);base=path.parent
    if s['status']!='passed':raise ValueError('Native study is incomplete.')
    for name,digest in s['artifacts_sha256'].items():
        f=(base/name).resolve()
        if base.resolve() not in f.parents or sha(f)!=digest:raise ValueError('Native artifact changed: '+name)
    for name,digest in s['source_sha256'].items():
        if sha(verified_source_path(name,digest))!=digest:raise ValueError('OpenLB application changed; rerun numerical verification.')
    cases={name:analyze(base/name) for name in s['cases']}
    close_tree(cases,s['cases'],'native fields')
    from .geometry import polygon
    jl=runtime();vertices=jl.seval('Matrix{Float64}')(polygon())
    for reference in s['references']:
        v=jl.duct_stokes_reference(vertices,reference['resolution'])
        replay={'resolution':reference['resolution'],**{str(k):float(getattr(v,str(k))) for k in jl.fieldnames(jl.typeof(v))}}
        close_tree(replay,reference,'independent duct-reference replay')
    close_tree(verify_suite(cases,s['references'],s['protocol']),s['checks'],'native checks')
    if (base/'operating.json').exists():
        operating=read(base/'operating.json')
        for name,entry in operating['cases'].items():
            actual=analyze(base/name);close_tree(actual,entry,'operating fields')
        if not all(operating['checks'].values()):raise ValueError('Operating-pressure verification failed.')
    return {'status':'passed','native_cases':len(cases),'original_vendor':verify_vendor()}

def verify(path,upstream=True):
    path=Path(path).resolve();base=path.parent;s=read(path)
    if s['status']!='passed':raise ValueError('Incomplete study.')
    files={str(f.relative_to(base)) for f in base.rglob('*') if f.is_file() and f.name not in ('verification.json','status.json') and f!=path}
    if files!=set(s['artifacts_sha256']):raise ValueError('Process artifact coverage changed.')
    for name,digest in s['artifacts_sha256'].items():
        f=(base/name).resolve()
        if base not in f.parents or sha(f)!=digest:raise ValueError('Changed process artifact: '+name)
    for name,value in read(base/'inputs.json').items():close_tree(s[name],value,'input snapshot.'+name)
    native=read(base/'native_summary.json');parent=read(base/'parent_summary.json');source=read(base/'source_summary.json')
    if source['protocol']['numerics'].get('collision_model')!='TRT' or any(v['raw'].get('collision_model')!='TRT' for v in source['native_cases'].values()):
        raise ValueError('Current porous-support evidence must use TRT.')
    for name,key in [('native_summary.json','native_sha256'),('parent_summary.json','parent_sha256'),('source_summary.json','source_sha256')]:
        if sha(base/name)!=s[key]:raise ValueError('Changed ancestor snapshot.')
    if upstream:
        close_tree(verify_native(Path(s['native_path'])),s['native_verification'],'native ancestry')
        from examples.filtration.literature_benchmark.check import verify as verify_parent
        close_tree(verify_parent(Path(s['parent_path'])),s['parent_verification'],'parent ancestry')
        from examples.filtration.graded_support.check import verify as verify_support
        close_tree(verify_support(Path(s['source_path']),write_record=False),s['support_verification'],'TRT support ancestry')
    jl=runtime();expected=list(requests(s['protocol']));covered=[]
    for c in s['cases']+[s['seed_validation']]:
        q=c['request'];t=translate(jl,native,source,s['protocol'],q)
        close_tree(t,c['translation'],'Julia replay')
        is_seed=c is s['seed_validation']
        seed={'seed_multiplier':.99,'tolerance_multiplier':.1,'seed_temperature_offset':2.} if is_seed else {}
        a=inputs(s['hydraulics'],source,t,**seed)
        directory=base/('validation/seed' if is_seed else 'process/'+c['id'])
        actual=read(directory/'summary.json');close_tree(actual,c['dwsim'],'saved DWSIM')
        close_tree(a,read(directory/'inputs.json'),'DWSIM request')
        close_tree(validate(actual,a,t['process'],s['hydraulics'],source),c['validation'],'DWSIM balances')
        close_tree(topology(directory/'crossflow.dwxmz'),c['saved_connections'],'DWSIM topology')
        if not is_seed:
            exported=read(directory/'export/summary.json')
            close_tree(exported['input_snapshot'],q,'exported input snapshot')
            close_tree(exported['case'],c,'exported completed case')
            if exported['case']['translation']['critical_flux_limit_LMH'] is not None:
                raise ValueError('Unavailable sustainable flux must stay JSON null.')
            with (directory/'export/axial_profile.csv').open(newline='') as stream:profile=list(csv.DictReader(stream))
            a=t['axial'];f=t['feed']
            expected_profile=list(zip(a['x_m'],a['pressure_Pa'],a['local_flux_m_s'],a['channel_flow_m3_s'],f['axial_COD_kg_m3'],strict=True))
            if len(profile)!=len(expected_profile):raise ValueError('Exported axial array shape differs.')
            for row,values in zip(profile,expected_profile,strict=True):
                close_tree([float(row[k]) for k in ('x_m','TMP_Pa','local_flux_m_s','channel_flow_positive_x_m3_s','COD_kg_m3')],list(values),'axial CSV')
        if not is_seed:covered.append(q)
        f=t['feed'];scale=t['axial']['module_flow_m3_s']*f['module_inlet_COD_kg_m3']
        if max(abs(f[k]) for k in f if k.endswith('_residual_kg_s'))>1e-10*scale:
            raise ValueError('COD balance failed.')
    for rejected in s['rejected_scenarios']:
        try:translate(jl,native,source,s['protocol'],rejected['request'])
        except ValueError as e:
            if str(e)!=rejected['reason']:raise ValueError('Scenario rejection changed.')
        else:raise ValueError('Scenario no longer rejected.')
        covered.append(rejected['request'])
    if sorted(map(str,expected))!=sorted(map(str,covered)):raise ValueError('Missing or duplicate process scenarios.')
    from .run import comparison
    close_tree(comparison(s['cases']),s['comparison'],'process comparisons')
    with (base/'comparison.csv').open(newline='') as stream:rows=list(csv.DictReader(stream))
    if len(rows)!=len(s['comparison']):raise ValueError('Comparison CSV row count differs.')
    for row,expected in zip(rows,s['comparison'],strict=True):
        if set(row)!=set(expected):raise ValueError('Comparison CSV schema differs.')
        for key,value in expected.items():close_tree(row[key] if isinstance(value,str) else float(row[key]),value,'comparison CSV.'+key)
    nominal=next(c for c in s['cases'] if c['request']==s['seed_validation']['request'])
    for k in ('feed_W','circulation_W','cooling_removed_W'):
        if not math.isclose(nominal['dwsim']['powers'][k],s['seed_validation']['dwsim']['powers'][k],rel_tol=1e-8):
            raise ValueError('Seed invariance failed.')
    if s['nonconvergence']['nonconvergence_rejected'] is not True:raise ValueError('Missing failed recycle check.')
    if 'RECYCLE_NOT_CONVERGED' not in (base/'validation/nonconvergence/stdout.log').read_text():raise ValueError('Missing failed-solve log.')
    return {'status':'passed','DWSIM_saved_reloaded_cases':len(s['cases'])+1,
            'rejected_scenarios':len(s['rejected_scenarios']),'artifacts_checked':len(s['artifacts_sha256']),
            'upstream_rechecked':upstream,'external_solves_by_checker':False}

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('summary',nargs='?',type=Path)
    q=parser.parse_args();from pore2process_demo.io import write_json
    import json
    print(json.dumps(verify(q.summary or Path(read(ROOT/'results/star_channel/process/latest.json')['summary'])),indent=2))
