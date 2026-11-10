"""Run and retain the explicit-geometry biological-feed process matrix."""
from pore2process_demo.paths import project_path
from pore2process_demo.provenance import implementation_files, source_manifest
import argparse
import csv
from datetime import datetime, timezone
from pathlib import Path
import shutil
from uuid import uuid4

from examples.filtration.calibration.audit_crossflow_sources import audit
from pore2process_demo.io import sha
from examples.filtration.crossflow_recycle import dwsim
from examples.filtration.crossflow_recycle.model import validate, topology
from pore2process_demo.io import write_json
from pore2process_demo.bridge import runtime, runtime_versions
from pore2process_demo.paths import ROOT; from pore2process_demo.io import read_json as read
from .model import settings, requests, translate, inputs, key

def evidence():
    native_path=Path(read(ROOT/'results/star_channel/native/latest.json')['summary'])
    parent_path=Path(read(ROOT/'results/crossflow_benchmark/latest.json')['summary'])
    return native_path,parent_path

def support_evidence():
    path=Path(read(ROOT/'results/graded_support_trt/latest.json')['summary'])
    source=read(path)
    if source.get('status')!='passed' or source['protocol']['numerics'].get('collision_model')!='TRT' or any(
            c['raw'].get('collision_model')!='TRT' for c in source['native_cases'].values()):
        raise ValueError('The new process case requires revalidated TRT support measurements.')
    return path

def solve_one(directory,request,native,source,p,h,jl=None,**seed):
    jl=runtime() if jl is None else jl
    t=translate(jl,native,source,p,request)
    a=inputs(h,source,t,**seed); actual=dwsim.solve(directory,a)
    checks=validate(actual,a,t['process'],h,source)
    for name,value in t['feed'].items():
        if name.endswith('_residual_kg_s') and abs(value)>1e-15:
            raise ValueError('Distributed COD conservation failed.')
    if abs(t['axial']['distributed_port_energy_residual_W'])>1e-5*checks['total_shaft_power_W']:
        raise ValueError('Distributed hydraulic-port energy identity failed.')
    return {'id':key(request),'request':request,'translation':t,'dwsim':actual,
            'validation':checks,'saved_connections':topology(directory/'crossflow.dwxmz')}

def comparison(cases):
    rows=[]
    for c in cases:
        q=c['request'];a=c['translation']['axial'];f=c['translation']['feed'];v=c['validation']
        baseline=next(vv for vv in cases if vv['request']=={**q,'support':'straight'})
        rows.append({'case_id':c['id'],**q,'inlet_Re':a['inlet_Re'],'inlet_velocity_m_s':a['inlet_velocity_m_s'],
            'mean_TMP_Pa':a['mean_TMP_Pa'],'inlet_TMP_Pa':a['inlet_TMP_Pa'],'outlet_TMP_Pa':a['outlet_TMP_Pa'],
            'channel_drop_Pa':a['channel_drop_Pa'],'fixed_resistance_fraction':a['selective_and_other_fixed_fraction'],
            'maximum_local_flux_LMH':c['translation']['maximum_local_flux_LMH'],
            'fresh_COD_mg_L':f['fresh_COD_kg_m3']*1000,'permeate_COD_mg_L':f['permeate_COD_kg_m3']*1000,
            'retentate_COD_mg_L':f['retentate_COD_kg_m3']*1000,'overall_COD_removal_fraction':f['overall_COD_removal_fraction'],
            'feed_power_W':c['dwsim']['powers']['feed_W'],'circulation_power_W':c['dwsim']['powers']['circulation_W'],
            'total_power_W':v['total_shaft_power_W'],'specific_energy_kWh_m3':v['specific_shaft_energy_kWh_m3_reference_volume'],
            'support_energy_saving_fraction':1-v['total_shaft_power_W']/baseline['validation']['total_shaft_power_W']})
    return rows

def execute(output,require_sources=False):
    npth,ppth=evidence();spth=support_evidence()
    native=read(npth);parent=read(ppth);source=read(spth)
    p=read(ROOT/'examples/filtration/star_channel/protocol.json');h=settings(parent);jl=runtime()
    from .check import verify_native
    native_verification=verify_native(npth)
    from examples.filtration.graded_support.check import verify as verify_support
    support_verification=verify_support(spth,write_record=False)
    from examples.filtration.literature_benchmark.check import verify as verify_parent
    print('Rechecking upstream evidence...',flush=True); parent_verification=verify_parent(ppth)
    source_audit=audit(read(ROOT/'literature/reviews/crossflow/extraction.json'),require_sources)
    bio=read(ROOT/'literature/sources/chiu2006/bio_extraction.json')
    if (p['feed']['COD_mg_L']!=bio['measurements'][0]['value'] or
        p['feed']['reference_rejection']!=bio['measurements'][2]['value']):
        raise ValueError('Protocol differs from the retained biological source extraction.')
    base=Path(output)/'runs'/(datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid4().hex[:8])
    base.mkdir(parents=True)
    implementation=list((ROOT/'examples/filtration/star_channel').glob('*.py'))+[project_path(ROOT/name) for name in (
        'examples/filtration/star_channel/protocol.json','examples/filtration/star_channel/support_protocol.json','julia/src/StarChannel.jl',
        'julia/src/GradedSupport.jl','julia/src/Pore2ProcessDemo.jl','julia/Project.toml','julia/Manifest.toml',
        'uv.lock','src/pore2process_demo/bridge.py','src/pore2process_demo/filtration_ui.py',
        'src/pore2process_demo/ui.py','src/pore2process_demo/visualization.py',
        'examples/filtration/crossflow_recycle/model.py','examples/filtration/crossflow_recycle/dwsim.py','examples/filtration/crossflow_recycle/dwsim_runner/Program.cs','app.py')]
    implementation = implementation_files(implementation)
    snapshot={'protocol':p,'hydraulics':h,'native_path':str(npth),'parent_path':str(ppth),
        'native_sha256':sha(npth),'parent_sha256':sha(ppth),'source_path':str(spth),'source_sha256':sha(spth),
        'source_audit':source_audit,'bio_extraction':bio,'versions':runtime_versions(),
        'implementation_sha256':source_manifest(sorted(implementation))}
    write_json(base/'inputs.json',snapshot)
    shutil.copy2(npth,base/'native_summary.json');shutil.copy2(ppth,base/'parent_summary.json')
    shutil.copy2(spth,base/'source_summary.json')
    adapter=dwsim.build(base/'adapter_build');cases=[];rejections=[]
    for request in requests(p):
        try: translate(jl,native,source,p,request)
        except ValueError as e:
            print('Rejected scenario '+key(request)+': '+str(e),flush=True)
            rejections.append({'request':request,'reason':str(e)});continue
        print('DWSIM solve/save/reload: '+key(request),flush=True)
        case=solve_one(base/'process'/key(request),request,native,source,p,h,jl)
        from .export import export_completed
        export_completed(base/'process'/key(request)/'export',case,npth,spth,p,runtime_versions())
        cases.append(case)
    # A strongly perturbed tear checks genuine iteration at high circulation.
    nominal=next(c for c in cases if c['request']==dict(support='straight',flux_LMH=p['process']['nominal_flux_LMH'],rejection=.65,viscosity_multiplier=1.,fixed_multiplier=1.,inlet_Re=p['process']['inlet_Re']))
    print('DWSIM seed/tolerance check',flush=True)
    validation_case=solve_one(base/'validation/seed',nominal['request'],native,source,p,h,jl,
        seed_multiplier=.99,tolerance_multiplier=.1,seed_temperature_offset=2.)
    for power in ('feed_W','circulation_W','cooling_removed_W'):
        import math
        if not math.isclose(validation_case['dwsim']['powers'][power],nominal['dwsim']['powers'][power],rel_tol=1e-8):
            raise ValueError('Seed/tolerance invariance failed.')
    bad=dict(nominal['dwsim']['inputs']);bad['Maximum_iterations']=1;bad['Seed_module_mass_flow_kg_s']*=.2
    failed=dwsim.reject_nonconvergence(base/'validation/nonconvergence',bad)
    # Prototype bounds are separate from the unavailable sustainable-flux limit.
    try: translate(jl,native,source,p,{**nominal['request'],'flux_LMH':100.})
    except ValueError as e: unsupported_flux_rejection=str(e)
    else: raise ValueError('Out-of-range prototype flux was not rejected.')
    rows=comparison(cases)
    with (base/'comparison.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    from .plot import plot_study
    plot_study(native,cases,rows,base/'study.png')
    result={'status':'passed',**snapshot,'cases':cases,'comparison':rows,'rejected_scenarios':rejections,
        'seed_validation':validation_case,'nonconvergence':failed,'unsupported_flux_rejection':unsupported_flux_rejection,
        'native_verification':native_verification,'support_verification':support_verification,'parent_verification':parent_verification,'adapter':adapter,
        'timestamp_utc':datetime.now(timezone.utc).isoformat()}
    result['artifacts_sha256']={str(f.relative_to(base)):sha(f) for f in base.rglob('*') if f.is_file()}
    write_json(base/'summary.json',result)
    from .check import verify
    checked=verify(base/'summary.json',upstream=False);write_json(base/'verification.json',checked)
    write_json(Path(output)/'latest.json',{'summary':str((base/'summary.json').resolve())})
    print('Verified process study: '+str(base),flush=True);return result

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,default=ROOT/'results/star_channel/process')
    parser.add_argument('--require-sources',action='store_true')
    args=parser.parse_args();execute(args.output,args.require_sources)
