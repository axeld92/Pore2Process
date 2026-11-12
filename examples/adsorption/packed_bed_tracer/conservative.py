"""One tracer-derived D with actual flow and physical liquid inventory fixed."""
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path
import shutil
from uuid import uuid4
import numpy as np
from pore2process_demo.paths import ROOT
from pore2process_demo.io import sha,write_json
from .bridge import call,moments
from .check import hashes

HERE=Path(__file__).parent

def derive(case):
    L=case['transport_length_m'];r=case['raw'];Q=r['flow_m3_s'];V=r['fluid_volume_m3']
    v=Q*L/V;target=case['moments']['mean_s']
    # Broad, declared Pe interval of the independently checked analytical family.
    lo=.05;hi=1e5
    def mean(pe):return call('absorbing_moments',length_m=L,velocity_m_s=v,dispersion_m2_s=v*L/pe)['mean_s']
    if not mean(lo)<target<mean(hi):raise ValueError('Short native mean not representable with fixed physical velocity.')
    for _ in range(70):
        p=(lo*hi)**.5
        if mean(p)<target:lo=p
        else:hi=p
    pe=(lo*hi)**.5
    return dict(effective_dispersion_m2_s=v*L/pe,calibration_case='bed32',holdout_case='length32',
        physical_velocity_m_s=v,physical_fluid_area_m2=V/L,flow_m3_s=Q,fluid_volume_m3=V,
        fitted_quantities=['effective dispersion from short inert-tracer mean only'],
        velocity_fitted=False,inventory_fitted=False,adsorption_data_used=False)

def predict(native,case,derived,*,cells=800,dt_s=.002):
    r=case['raw'];L=case['transport_length_m'];Q=r['flow_m3_s'];V=r['fluid_volume_m3'];area=V/L;v=Q/area
    flux=np.genfromtxt(native/'native'/case['id']/'flux.csv',delimiter=',',names=True)
    p=dict(length_m=L,velocity_m_s=v,dispersion_m2_s=derived['effective_dispersion_m2_s'],
        t_end_s=float(flux['time_s'][-1]),dt_s=dt_s,cells=cells)
    model=call('solve_ade',**p);a=call('absorbing_moments',length_m=L,velocity_m_s=v,dispersion_m2_s=p['dispersion_m2_s'])
    model_fraction=np.interp(flux['time_s'],model['time_s'],np.array(model['outlet_flux_per_area_m_s'])/a['steady_flux_per_area_m_s'])
    observed=flux['outlet_flux_kg_s']/case['moments']['steady_exit_flux_kg_s']
    steady=area*a['steady_flux_per_area_m_s']
    metrics=dict(normalized_profile_RMSE=float(np.sqrt(np.mean((model_fraction-observed)**2))),
        mean_relative_error=abs(a['mean_s']/case['moments']['mean_s']-1),
        variance_relative_error=abs(a['variance_s2']/case['moments']['variance_s2']-1),
        steady_exit_flux_relative_error=abs(steady/case['moments']['steady_exit_flux_kg_s']-1))
    metrics['passed']=metrics['normalized_profile_RMSE']<=.03 and metrics['mean_relative_error']<=.1 and metrics['variance_relative_error']<=.2 and metrics['steady_exit_flux_relative_error']<=.1
    numerical=moments(model['time_s'],model['outlet_flux_per_area_m_s'])
    return dict(inputs=p,physical_fluid_area_m2=area,physical_fluid_volume_m3=V,flow_m3_s=Q,
        numerical_moments=numerical,
        predicted_steady_exit_flux_kg_s=steady,analytical_moments=a,result=model,metrics=metrics)

def refinements(a,b):
    f=np.array(a['result']['outlet_flux_per_area_m_s'])/a['analytical_moments']['steady_flux_per_area_m_s']
    g=np.interp(a['result']['time_s'],b['result']['time_s'],np.array(b['result']['outlet_flux_per_area_m_s'])/b['analytical_moments']['steady_flux_per_area_m_s'])
    error=float(np.max(abs(f-g)))
    dm=abs(b['numerical_moments']['mean_s']/a['numerical_moments']['mean_s']-1)
    dv=abs(b['numerical_moments']['variance_s2']/a['numerical_moments']['variance_s2']-1)
    return dict(maximum_normalized_profile_change=error,mean_relative_change=dm,
        variance_relative_change=dv,passed=error<=.005 and dm<=.02 and dv<=.02)

def run():
    ptr=json.loads((ROOT/'results/adsorption/packed_bed_tracer/latest.json').read_text())
    path=ROOT/ptr['summary'];s=json.loads(path.read_text());native=path.parent
    hashes(ROOT,s['source_sha256']);hashes(native,s['artifacts_sha256'])
    if not s['native_checks_passed'] or not s['numerical_refinement_passed']:raise ValueError('Native prerequisite is unqualified.')
    rid=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid4().hex[:8]
    out=native/'translations'/rid;out.mkdir(parents=True)
    sources={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),HERE/'conservative_protocol.md',
        HERE/'bridge.py',HERE/'check.py',ROOT/'julia/adsorption/TracerContact.jl']}
    for name in sources:
        target=out/'source_snapshot'/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/name,target)
    d=derive(s['cases']['bed32']);predictions={};refinement={}
    for key in ('bed32','length32'):
        case={**s['cases'][key],'id':key};predictions[key]=predict(native,case,d)
        write_json(out/f'{key}.json',predictions[key]);refinement[key]={}
        for label,cells,dt in [('space',1600,.002),('time',800,.001)]:
            p=predict(native,case,d,cells=cells,dt_s=dt);write_json(out/f'{key}_{label}.json',p)
            refinement[key][label]=refinements(predictions[key],p)
    summary=dict(schema='pore2process.conservative-tracer.v1',run_id=rid,native_summary=str(path.relative_to(ROOT)),
        native_summary_sha256=sha(path),native_run_id=s['run_id'],source_sha256=sources,derived=d,
        metrics={key:p['metrics'] for key,p in predictions.items()},refinement=refinement,
        numerical_refinement_passed=all(v['passed'] for r in refinement.values() for v in r.values()),
        comparison_passed=all(p['metrics']['passed'] for p in predictions.values()),
        new_native_runs=0,new_DWSIM_runs=0,physical_inventory_substitution=False)
    summary['artifacts_sha256']={str(p.relative_to(out)):sha(p) for p in out.rglob('*') if p.is_file()}
    write_json(out/'summary.json',summary);write_json(native/'conservative_latest.json',dict(summary=str((out/'summary.json').relative_to(ROOT))))
    print(json.dumps({k:v for k,v in summary.items() if k not in ['source_sha256','artifacts_sha256']},indent=2))

def check():
    p=json.loads((ROOT/'results/adsorption/packed_bed_tracer/latest.json').read_text())
    native=(ROOT/p['summary']).parent
    path=ROOT/json.loads((native/'conservative_latest.json').read_text())['summary'];out=path.parent;s=json.loads(path.read_text())
    hashes(ROOT,s['source_sha256']);hashes(out/'source_snapshot',s['source_sha256']);hashes(out,s['artifacts_sha256'])
    assert sha(ROOT/s['native_summary'])==s['native_summary_sha256']
    n=json.loads((ROOT/s['native_summary']).read_text());hashes((ROOT/s['native_summary']).parent,n['artifacts_sha256'])
    d=derive(n['cases']['bed32']);assert d==s['derived']
    for key in ('bed32','length32'):
        baseline=None
        for label,cells,dt in [('',800,.002),('_space',1600,.002),('_time',800,.001)]:
            saved=json.loads((out/f'{key}{label}.json').read_text())
            actual=predict(native,{**n['cases'][key],'id':key},d,cells=cells,dt_s=dt)
            assert actual==saved
            r=actual['result'];area=actual['physical_fluid_area_m2'];p=actual['inputs']
            assert np.isclose(area*p['length_m'],actual['physical_fluid_volume_m3'],rtol=1e-13)
            assert np.isclose(area*p['velocity_m_s'],actual['flow_m3_s'],rtol=1e-13)
            assert r['normalized_mass_residual']<1e-8
            residue=area*(np.array(r['inventory_per_area_m'])-np.array(r['cumulative_inlet_per_area_m'])+np.array(r['cumulative_outlet_per_area_m']))
            assert np.max(abs(residue))/(area*r['cumulative_inlet_per_area_m'][-1])<1e-8
            if not label:
                baseline=actual;assert actual['metrics']==s['metrics'][key]
            else:assert refinements(baseline,actual)==s['refinement'][key][label[1:]]
    reference={k:s[k] for k in ['schema','run_id','native_run_id','native_summary','derived','metrics','refinement',
        'numerical_refinement_passed','comparison_passed','physical_inventory_substitution','new_native_runs','new_DWSIM_runs']}
    reference['predictions']={key:{k:v for k,v in json.loads((out/f'{key}.json').read_text()).items() if k!='result'} for key in ('bed32','length32')}
    write_json(HERE/'reference/conservative_reference.json',reference)
    report=dict(Julia_predictions_replayed=6,comparison_passed=s['comparison_passed'],
        numerical_refinement_passed=s['numerical_refinement_passed'],flow_and_liquid_inventory_preserved=True,
        source_files_checked=len(s['source_sha256']),artifacts_checked=len(s['artifacts_sha256']))
    write_json(out/'verification.json',report);print(json.dumps(report,indent=2))

def main():
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');args=p.parse_args()
    check() if args.check else run()

if __name__=='__main__':main()
