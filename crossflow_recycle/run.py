"""Actual DWSIM recycle scenarios from the preserved OpenLB/Julia data boundary."""
import argparse
import csv
from datetime import datetime,timezone
from pathlib import Path
import platform
import shutil
from time import perf_counter
from uuid import uuid4

from coupling.run import sha
from hydraulic.run import write_json
from pore2process_demo.bridge import runtime,runtime_versions
from process_filtration.check import verify as verify_parent
from process_filtration.run import read,ROOT
from crossflow_recycle import dwsim
from crossflow_recycle.model import translate,adapter_inputs,validate,topology

def check_protocol(p,source):
    if p['protocol_id']!='crossflow-recycle-hydraulics-v1' or source['protocol']['protocol_id']!=p['source_protocol_id']:
        raise RuntimeError('Review changed study/source protocol')
    q=source['protocol']['process']
    for key,value in (('temperature_K',q['temperature_K']),('external_absolute_pressure_Pa',q['inlet_absolute_pressure_Pa']),
                      ('permeate_flux_LMH',q['filtration_flux_LMH']),('membrane_area_m2',q['membrane_area_m2']),
                      ('feed_pump_efficiency_fraction',q['pump_efficiency_fraction'])):
        if p[key]!=value:raise RuntimeError('Do not silently reuse native data at changed anchor conditions')

def jobs(p,source):
    conditions=[('nominal',p['nominal_recovery'],p['nominal_module_to_permeate_ratio'],1.)]
    conditions += [(f'circulation_{c:g}',p['nominal_recovery'],c,1.) for c in p['module_to_permeate_ratios'] if c!=p['nominal_module_to_permeate_ratio']]
    conditions += [(f'recovery_{100*y:g}',y,p['nominal_module_to_permeate_ratio'],1.) for y in p['overall_recoveries'] if y!=p['nominal_recovery']]
    conditions += [(f'channel_{m:g}',p['nominal_recovery'],p['nominal_module_to_permeate_ratio'],m) for m in p['channel_resistance_multipliers'] if m!=1.]
    for label,y,c,m in conditions:
        for support in source['supports']:
            yield {'case_id':label+'/'+support['name'],'condition':label,'support':support['name'],
                'recovery':y,'module_ratio':c,'resistance_multiplier':m,'seed_options':{},'validation_only':False}
    r=p['recycle']
    for seed,offset in zip(r['validation_seed_flow_multipliers'],r['validation_seed_temperature_offsets_K'],strict=True):
        yield {'case_id':f'validation/seed_{seed:g}','condition':'nominal','support':'straight','recovery':p['nominal_recovery'],
            'module_ratio':p['nominal_module_to_permeate_ratio'],'resistance_multiplier':1.,
            'seed_options':{'seed_multiplier':seed,'seed_temperature_offset':offset},'validation_only':True}
    high=max(p['module_to_permeate_ratios'])
    for name,options in (('high_circulation_seed',{'seed_multiplier':.1,'seed_temperature_offset':2.}),
                         ('tight_tolerance',{'tolerance_multiplier':r['validation_tolerance_multiplier']})):
        yield {'case_id':'validation/'+name,'condition':f'circulation_{high:g}','support':'straight',
            'recovery':p['nominal_recovery'],'module_ratio':high,'resistance_multiplier':1.,'seed_options':options,'validation_only':True}

def comparison_rows(cases):
    return [{'case_id':c['request']['case_id'],'support':c['request']['support'],'condition':c['request']['condition'],
        'validation_only':int(c['request']['validation_only']),
        'overall_recovery_fraction':c['julia_process']['overall_recovery_fraction'],
        'module_to_permeate_ratio':c['julia_process']['module_to_permeate_ratio'],
        'channel_resistance_Pa_s_m3':c['julia_process']['channel_resistance_Pa_s_m3'],
        'channel_pressure_drop_Pa':c['julia_process']['channel_pressure_drop_Pa'],
        'TMP_Pa':c['julia_process']['transmembrane_pressure_Pa'],
        'fresh_reference_flow_m3_s':c['julia_process']['fresh_feed_reference_flow_m3_s'],
        'module_reference_flow_m3_s':c['julia_process']['module_reference_flow_m3_s'],
        'recycle_reference_flow_m3_s':c['julia_process']['recycle_reference_flow_m3_s'],
        'actual_feed_power_W':c['dwsim']['powers']['feed_W'],
        'actual_circulation_power_W':c['dwsim']['powers']['circulation_W'],
        'actual_total_shaft_power_W':c['validation']['total_shaft_power_W'],
        'actual_cooling_removed_W':c['dwsim']['powers']['cooling_removed_W'],
        'actual_specific_shaft_energy_kWh_m3_reference_volume':c['validation']['specific_shaft_energy_kWh_m3_reference_volume'],
        'circulation_shaft_fraction':c['validation']['circulation_shaft_fraction'],
        'energy_saving_vs_straight_fraction':c['energy_saving_vs_straight_fraction'],
        'recycle_iterations':c['dwsim']['recycle']['iterations_taken']} for c in cases]

def execute(output,parent_path):
    p=read(ROOT/'crossflow_recycle/protocol.json');parent_path=Path(parent_path).resolve();parent=read(parent_path)
    source_path=ROOT/parent['source_summary_path'];source=read(source_path);check_protocol(p,source)
    t=perf_counter();jl=runtime();cold=perf_counter()-t
    print('Rechecking earlier filtration and immutable OpenLB evidence...',flush=True)
    verified_parent=verify_parent(parent_path)
    run_id=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid4().hex[:8]
    base=Path(output)/'runs'/run_id;base.mkdir(parents=True)
    write_json(base/'status.json',{'status':'running','run_id':run_id})
    try:
        shutil.copy2(parent_path,base/'parent_summary.json');shutil.copy2(source_path,base/'source_summary.json')
        write_json(base/'parent_verification.json',verified_parent)
        adapter=dwsim.build(base/'adapter_build')
        paths=sorted(f for f in (ROOT/'crossflow_recycle').rglob('*') if f.suffix in ('.py','.cs','.csproj','.json'))
        paths += [ROOT/f for f in ('julia/src/CrossflowRecycle.jl','julia/src/Pore2ProcessDemo.jl',
            'julia/src/GradedSupport.jl','julia/src/SeriesCalibration.jl','graded_support/model.py','graded_support/check.py',
            'process_filtration/check.py','process_filtration/run.py','process_filtration/model.py',
            'compatibility/run_dwsim.py','src/pore2process_demo/bridge.py','uv.lock','julia/Project.toml','julia/Manifest.toml')]
        snapshot={'protocol':p,'parent_summary_path':str(parent_path.relative_to(ROOT)),
            'parent_summary_sha256':sha(parent_path),'source_summary_sha256':sha(source_path),
            'source_run_id':source['run_id'],'parent_run_id':parent['run_id'],
            'implementation_sha256':{str(f.relative_to(ROOT)):sha(f) for f in paths},
            'runtime_sha256':{str(f.relative_to(ROOT)):sha(f) for f in (dwsim.DISTRIBUTION/name for name in
                ('DWSIM.Automation.FluentAPI.dll','DWSIM.UnitOperations.dll','DWSIM.Thermodynamics.dll',
                 'DWSIM.Interfaces.dll','DWSIM.FlowsheetSolver.dll','DWSIM.SharedClasses.dll','DWSIM.Automation.dll'))},
            'machine':{'system':platform.system(),'release':platform.release(),'architecture':platform.machine()}}
        write_json(base/'inputs.json',snapshot)
        cases=[]
        for request in jobs(p,source):
            support=next(v for v in source['supports'] if v['name']==request['support'])
            t=perf_counter();membrane,c=translate(jl,source,support,p,request['recovery'],request['module_ratio'],request['resistance_multiplier'])
            call_time=perf_counter()-t;inputs=adapter_inputs(p,source,c,**request['seed_options'])
            print('Solving/reloading actual DWSIM recycle: '+request['case_id'],flush=True)
            directory=base/'process'/request['case_id'];actual=dwsim.solve(directory,inputs)
            v=validate(actual,inputs,c,p,source);connections=topology(directory/'crossflow.dwxmz')
            baseline=next((b for b in cases if b['request']['condition']==request['condition'] and b['request']['support']=='straight' and not b['request']['validation_only']),None)
            saving=0. if baseline is None else 1-v['total_shaft_power_W']/baseline['validation']['total_shaft_power_W']
            cases.append({'request':request,'julia_membrane':membrane,'julia_process':c,'dwsim':actual,
                'validation':v,'saved_connections':connections,'energy_saving_vs_straight_fraction':saving,'julia_translation_s':call_time})
        failed_input=adapter_inputs(p,source,cases[0]['julia_process'],seed_multiplier=.1)
        failed_input['Maximum_iterations']=1
        print('Checking deliberate nonconvergence rejection...',flush=True)
        rejection=dwsim.reject_nonconvergence(base/'validation/nonconvergence',failed_input)
        with (base/'comparison.csv').open('w',newline='') as f:
            rows=comparison_rows(cases);w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
        from crossflow_recycle.plot import plot
        plot(cases,base/'circulation_energy.png')
        shutil.copy2(base/'process/nominal/straight/flowsheet.png',base/'flowsheet.png')
        result={'schema_version':1,'status':'passed','run_id':run_id,'timestamp_utc':datetime.now(timezone.utc).isoformat(),
            **snapshot,'parent_verification':verified_parent,'adapter':adapter,'versions':runtime_versions(),
            'Julia_initialization_s':cold,'cases':cases,'nonconvergence_rejection':rejection,
            'scope':'Clean-water steady-state hydraulic recycle with fixed reference flux, declared channel resistance and temperature control. Shaft work and removed heat are separate. No selectivity/fouling or total-plant electricity prediction.',
            'artifact_sha256':{str(f.relative_to(base)):sha(f) for f in base.rglob('*') if f.is_file() and f.name!='status.json'}}
        write_json(base/'summary.json',result)
        from crossflow_recycle.check import verify
        verification=verify(base/'summary.json',verify_upstream=False)
        write_json(base/'verification.json',verification);write_json(base/'status.json',{'status':'passed','run_id':run_id})
        write_json(Path(output)/'latest.json',{'summary':str((base/'summary.json').resolve())})
        print(f'Passed: {base/"summary.json"}',flush=True)
        return result
    except BaseException as e:
        write_json(base/'status.json',{'status':'failed','run_id':run_id,'error':str(e)});raise

def main():
    a=argparse.ArgumentParser(description=__doc__);a.add_argument('--parent',type=Path)
    a.add_argument('--output',type=Path,default=ROOT/'results/crossflow_recycle');args=a.parse_args()
    parent=args.parent or Path(read(ROOT/'results/process_filtration/latest.json')['summary'])
    execute(args.output,parent)

if __name__=='__main__':main()
