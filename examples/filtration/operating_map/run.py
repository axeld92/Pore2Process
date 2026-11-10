"""Actual DWSIM operating matrix, with bounded published-capacity comparisons."""
from pore2process_demo.provenance import implementation_files, source_manifest
import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
from datetime import datetime,timezone
from pathlib import Path
import shutil
from uuid import uuid4
from pore2process_demo.io import sha
from pore2process_demo.io import write_json
from pore2process_demo.bridge import runtime,runtime_versions
from pore2process_demo.paths import ROOT; from pore2process_demo.io import read_json as read
from examples.filtration.crossflow_recycle import dwsim
from examples.filtration.crossflow_recycle.model import validate, topology
from examples.filtration.star_channel.model import settings
from examples.filtration.star_channel.run import support_evidence, evidence as ancestors
from .evidence import extract
from .model import requests,key,translate,inputs,capacity_decisions,comparison,source_shape_diagnostics
from .selection import selections


def table(path,rows):
    if not rows:raise ValueError('Empty output table: '+str(path))
    with Path(path).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def solve_prepared(job):
    directory,q,t,a,h,source,decisions=job
    print('Actual DWSIM solve/save/reload: '+key(q),flush=True)
    if (directory/'summary.json').exists():
        actual=read(directory/'summary.json')
        if read(directory/'inputs.json')!=a:raise ValueError('Reused process request changed.')
        print('Replaying matching process output: '+key(q),flush=True)
    else:actual=dwsim.solve(directory,a)
    valid=validate(actual,a,t['process'],h,source)
    scale=t['axial']['module_flow_m3_s']*t['feed']['module_inlet_COD_kg_m3']
    for balance in [t['feed']]+[v['balance'] for v in t['quality_sensitivity']]:
        if max(abs(v) for k,v in balance.items() if k.endswith('_residual_kg_s'))>1e-10*scale:
            raise ValueError('External COD balance failed.')
    if abs(t['axial']['distributed_port_energy_residual_W'])>1e-5*valid['total_shaft_power_W']:
        raise ValueError('Distributed hydraulic energy identity failed.')
    a=t['axial'];f=t['feed']
    with (directory/'axial_profile.csv').open('w',newline='') as stream:
        w=csv.writer(stream);w.writerow(['x_m','TMP_Pa','local_flux_m_s','channel_flow_positive_x_m3_s','COD_kg_m3'])
        w.writerows(zip(a['x_m'],a['pressure_Pa'],a['local_flux_m_s'],a['channel_flow_m3_s'],f['axial_COD_kg_m3'],strict=True))
    result={'id':key(q),'request':q,'translation':t,'capacity':decisions,'dwsim':actual,
            'validation':valid,'saved_connections':topology(directory/'crossflow.dwxmz')}
    write_json(directory/'completed.json',result);return result


def execute(output,resume=None):
    npth=Path(read(ROOT/'results/operating_map/native/latest.json')['summary'])
    spth=support_evidence();_,ppth=ancestors();native,source=read(npth),read(spth)
    from .check import verify_native
    print('Rechecking native and TRT-support evidence...',flush=True)
    vn=verify_native(npth)
    from examples.filtration.graded_support.check import verify as verify_support
    vs=verify_support(spth,write_record=False)
    p=read(ROOT/'examples/filtration/operating_map/protocol.json');e=extract();jl=runtime();h=settings(read(ppth))
    p['scope']+=' Synthetic channel horizontal, gravity neglected and external permeate pressure uniform; source hydrostatic profiles are not transferred.'
    bio=e['biological_feed_extraction']['measurements']
    p['process']['module_inlet_COD_kg_m3']=bio[0]['value']*.001
    p['process']['module_length_m']=.3
    if p['process']['rejection']!=bio[2]['value']:raise ValueError('Nominal retention differs from the retained source extraction.')
    h['protocol_id']='qualified-higher-Re-water-carrier-v1'
    h['recycle']['maximum_iterations']=20000
    h['recycle']['nominal_seed_flow_multiplier']=1.
    h['recycle']['nominal_seed_temperature_offset_K']=2.
    base=Path(resume) if resume else Path(output)/'runs'/(datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid4().hex[:8])
    base=base.resolve();base.mkdir(parents=True,exist_ok=True)
    implementation=list((ROOT/'examples/filtration/operating_map').glob('*.py'))+[ROOT/name for name in (
        'examples/filtration/operating_map/protocol.json','julia/src/OperatingMap.jl','julia/src/StarChannel.jl',
        'julia/src/GradedSupport.jl','julia/src/CrossflowRecycle.jl','julia/src/Pore2ProcessDemo.jl',
        'examples/filtration/crossflow_recycle/model.py','examples/filtration/crossflow_recycle/dwsim_runner/Program.cs','examples/filtration/crossflow_recycle/dwsim.py',
        'uv.lock','julia/Project.toml','julia/Manifest.toml')]
    implementation = implementation_files(implementation)
    state={'protocol':p,'hydraulics':h,'evidence':e,'native_path':str(npth),'native_sha256':sha(npth),
        'support_path':str(spth),'support_sha256':sha(spth),'parent_path':str(ppth),'parent_sha256':sha(ppth),
        'implementation_sha256':source_manifest(implementation),'versions':runtime_versions()}
    if (base/'summary.json').exists():raise ValueError('Accepted run is immutable.')
    if (base/'inputs.json').exists() and read(base/'inputs.json')!=state:raise ValueError('Resume source or configuration changed.')
    write_json(base/'inputs.json',state)
    for file in implementation:
        target=base/'source_snapshot'/file.relative_to(ROOT);target.parent.mkdir(parents=True,exist_ok=True)
        if not target.exists():shutil.copy2(file,target)
    for path,name in [(npth,'native_summary.json'),(spth,'support_summary.json'),(ppth,'parent_summary.json')]:shutil.copy2(path,base/name)
    adapter=dwsim.build(base/'adapter_build');jobs=[];rejected=[]
    from .export import native_field
    if not (base/'native_field').exists():native_field(npth,base/'native_field')
    for q in requests(p):
        try:t=translate(jl,native,source,p,q)
        except ValueError as exc:
            rejected.append({'request':q,'reason':str(exc)});print('Rejected hydraulic scenario: '+key(q)+': '+str(exc),flush=True);continue
        decisions=capacity_decisions(jl,t,e,p)
        jobs.append((base/'process'/key(q),q,t,inputs(h,source,t),h,source,decisions))
    with ThreadPoolExecutor(max_workers=2) as pool:cases=list(pool.map(solve_prepared,jobs))
    nominal=next(c for c in cases if c['request']==dict(inlet_Re=1850.,flux_LMH=8.,fixed_multiplier=1.,support='straight'))
    q=nominal['request'];t=nominal['translation'];a=inputs(h,source,t,seed_multiplier=.5,tolerance_multiplier=.1,seed_temperature_offset=-1.)
    seed=solve_prepared((base/'validation/seed',q,t,a,h,source,nominal['capacity']))
    import math
    for k in ('feed_W','circulation_W','cooling_removed_W'):
        if not math.isclose(seed['dwsim']['powers'][k],nominal['dwsim']['powers'][k],rel_tol=1e-8):raise ValueError('Seed/tolerance invariance failed.')
    bad=dict(a);bad['Maximum_iterations']=1;bad['Seed_module_mass_flow_kg_s']*=.2
    failed=dwsim.reject_nonconvergence(base/'validation/nonconvergence',bad)
    rows=comparison(cases);chosen=selections(cases,p)
    table(base/'comparison.csv',rows);table(base/'selections.csv',chosen)
    write_json(base/'literature_evidence.json',e)
    from .plot import plot
    plot(native,cases,rows,chosen,e,p,base/'design_map.png')
    result={'status':'passed',**state,'cases':cases,'comparison':rows,'selections':chosen,
        'source_shape_diagnostics':source_shape_diagnostics(jl,e),'rejected_scenarios':rejected,'seed_validation':seed,'nonconvergence':failed,
        'native_verification':vn,'support_verification':vs,'adapter':adapter,
        'timestamp_utc':datetime.now(timezone.utc).isoformat(),
        'units':{'pressure':'Pa; TMP gauge relative to permeate, stream pressures absolute',
            'flow':'m3/s positive +x','flux':'m/s internally; LMH in named table columns','length':'m',
            'COD':'kg/m3 COD equivalent; mg/L in named table columns','energy':'kWh/m3 reference filtrate volume'},
        'scope':p['scope']}
    result['artifacts_sha256']={str(f.relative_to(base)):sha(f) for f in base.rglob('*') if f.is_file()}
    write_json(base/'summary.json',result)
    from .check import verify
    checked=verify(base/'summary.json',upstream=False);write_json(base/'verification.json',checked)
    write_json(Path(output)/'latest.json',{'summary':str(base/'summary.json')})
    print('Qualified operating map accepted: '+str(base),flush=True);return result

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,default=ROOT/'results/operating_map/process')
    parser.add_argument('--resume',type=Path);a=parser.parse_args();execute(a.output,a.resume)
