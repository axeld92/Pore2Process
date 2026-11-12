"""Execute and verify the separate stock flat-wall Poiseuille diagnostic."""
from datetime import datetime,timezone
import argparse
import json
from pathlib import Path
import shutil
import subprocess
from time import perf_counter
from uuid import uuid4

import numpy as np

from pore2process_demo.paths import ROOT
from pore2process_demo.io import sha,write_json
from pore2process_demo.provenance import source_manifest
from pore2process_demo.runtime.openlb import build_case,verify_vendor
from .reference import converged_reference

HERE=Path(__file__).parent
OUTPUT=ROOT/'results/adsorption/transport_root_cause/flat_channel'
PHYSICAL=dict(height_m=.0013891725543141453,mean_speed_m_s=.012937295182166382,
              diffusivity_m2_s=1e-9,dt_s=1.90174557964e-5,steps=11294)


def specifications():
    return [dict(id='full32',n=32,layers=32,amplitude=.25,speed_factor=1.),
            dict(id='thin32',n=32,layers=4,amplitude=.25,speed_factor=1.),
            dict(id='constant32',n=32,layers=4,amplitude=0.,speed_factor=1.),
            dict(id='diffusion32',n=32,layers=4,amplitude=.25,speed_factor=0.),
            dict(id='thin64',n=64,layers=4,amplitude=.25,speed_factor=1.),
            dict(id='thin96',n=96,layers=4,amplitude=.25,speed_factor=1.)]


def analyze(dest):
    dest=Path(dest);raw=json.loads((dest/'raw.json').read_text())
    history=np.genfromtxt(dest/'inventory.csv',delimiter=',',names=True)
    field=np.genfromtxt(dest/'final.csv',delimiter=',',names=True)
    n=raw['n'];H=raw['height_m'];amp=raw['initial_amplitude_kg_m3']
    lo=float(np.min(history['min_concentration_kg_m3']));hi=float(np.max(history['max_concentration_kg_m3']))
    excursion=max(.5-amp-lo,hi-(.5+amp),0.)
    inventory=history['kinetic_inventory_kg'];drift=float(np.max(abs(inventory/inventory[0]-1)))
    physical_drift=float(np.max(abs(history['fluid_inventory_kg']/history['fluid_inventory_kg'][0]-1)))
    expected_y=(field['j']+.5)*H/n;expected_x=(field['i']+.5)*H/n
    np.testing.assert_allclose(field['y_m'],expected_y,atol=1e-16,rtol=1e-13)
    np.testing.assert_allclose(field['x_m'],expected_x,atol=1e-16,rtol=1e-13)
    u=6*raw['mean_speed_m_s']*(field['y_m']/H)*(1-field['y_m']/H)
    np.testing.assert_allclose(field['ux_m_s'],u,atol=1e-14,rtol=1e-12)
    volume=H*H*raw['depth_m']
    mass_relative=abs(float(history['fluid_inventory_kg'][0])/(.5*volume)-1)
    wall_relative=abs(float(history['active_wall_population_kg'][0])/(.5*volume/(3*n))-1)
    checks=dict(finite=bool(all(np.all(np.isfinite(history[key])) for key in history.dtype.names)),
                kinetic_mass=drift<=1e-8,physical_liquid_mass=physical_drift<=1e-8,
                declared_fluid_volume=mass_relative<=1e-11,initial_active_wall_capacity=wall_relative<=1e-11,
                z_invariance=float(np.max(field['z_variation_kg_m3']))<=1e-12,
                strict_bounds=excursion<=1e-8,engineering_bounds=excursion<=.02,
                no_hidden_diffusivity=abs((raw['tau_odd']-.5)*(H/n)**2/(3*raw['dt_s'])/raw['diffusivity_m2_s']-1)<=1e-10)
    result=dict(raw=raw,checks=checks,minimum_concentration_kg_m3=lo,maximum_concentration_kg_m3=hi,
                maximum_bound_excursion_kg_m3=excursion,relative_kinetic_inventory_drift=drift,
                relative_physical_liquid_inventory_drift=physical_drift,
                initial_active_wall_storage_over_physical_liquid=1/(3*n),
                actual_input_snapshot=json.loads((dest/'request.json').read_text()),
                application_transport_qualified=False,measured_column_dispersion_m2_s=None)
    if amp:
        p=dict(height_m=H,mean_speed_m_s=raw['mean_speed_m_s'],diffusivity_m2_s=raw['diffusivity_m2_s'],duration_s=raw['duration_s'])
        exact,reference=converged_reference(field['x_m'],field['y_m'],amplitude_kg_m3=amp,**p)
        difference=field['concentration_kg_m3']-exact
        result.update(reference_convergence=reference,RMS_error_over_amplitude=float(np.sqrt(np.mean(difference*difference)))/amp,
                      maximum_error_over_amplitude=float(np.max(abs(difference)))/amp,
                      first_fluid_row_RMS_error_over_amplitude=float(np.sqrt(np.mean(difference[field['j']==0]**2)))/amp,
                      last_fluid_row_RMS_error_over_amplitude=float(np.sqrt(np.mean(difference[field['j']==n-1]**2)))/amp)
        checks['profile']=result['RMS_error_over_amplitude']<=.01
        if raw['mean_speed_m_s']==0:
            analytical=.5+amp*np.exp(-raw['diffusivity_m2_s']*(2*np.pi/H)**2*raw['duration_s'])*np.sin(2*np.pi*field['x_m']/H)
            result['oracle_analytical_max_error_kg_m3']=float(np.max(abs(exact-analytical)))
            checks['diffusion_analytical_oracle']=result['oracle_analytical_max_error_kg_m3']<=1e-10
    else:
        result['maximum_constant_error_kg_m3']=max(abs(lo-.5),abs(hi-.5))
        checks['constant_preservation']=result['maximum_constant_error_kg_m3']<=1e-9
    result['control_checks_passed']=all(checks.values())
    return result


def equivalence(full,thin):
    a_raw=json.loads((full/'raw.json').read_text());b_raw=json.loads((thin/'raw.json').read_text());comparisons={}
    for name in ('initial.csv','peak.csv','final.csv'):
        a=np.genfromtxt(full/name,delimiter=',',names=True);b=np.genfromtxt(thin/name,delimiter=',',names=True)
        if a.shape!=b.shape or a.dtype.names!=b.dtype.names:raise ValueError('Flat thin/full field shape differs')
        for key in ('i','j','x_m','y_m','ux_m_s'):np.testing.assert_array_equal(a[key],b[key])
        comparisons[name]=float(np.max(abs(a['concentration_kg_m3']-b['concentration_kg_m3'])))
        comparisons[name+'_z_variation']=float(np.max(b['z_variation_kg_m3']))
    a=np.genfromtxt(full/'inventory.csv',delimiter=',',names=True);b=np.genfromtxt(thin/'inventory.csv',delimiter=',',names=True)
    np.testing.assert_array_equal(a['time_s'],b['time_s']);np.testing.assert_array_equal(a['step'],b['step'])
    for name in ('min_concentration_kg_m3','max_concentration_kg_m3'):
        comparisons[name]=float(np.max(abs(a[name]-b[name])))
    for name in ('fluid_inventory_kg','active_wall_population_kg','kinetic_inventory_kg'):
        comparisons[name+'_depth_scaling_relative_error']=abs(float(b[name][0]/a[name][0])/(b_raw['depth_m']/a_raw['depth_m'])-1)
        comparisons[name+'_normalized_history_error']=float(np.max(abs(a[name]/a[name][0]-b[name]/b[name][0])))
    comparisons['peak_step_matches']=a_raw['peak_step']==b_raw['peak_step']
    comparisons['passed']=comparisons['peak_step_matches'] and all(v<=1e-11 for k,v in comparisons.items() if k!='peak_step_matches')
    return comparisons


def execute(base,spec,binary):
    dest=base/'native'/spec['id'];dest.mkdir(parents=True);p=PHYSICAL
    command=[str(binary),str(spec['n']),str(spec['layers']),str(p['height_m']),str(p['diffusivity_m2_s']),
             str(p['dt_s']),str(p['steps']),str(p['mean_speed_m_s']*spec['speed_factor']),str(spec['amplitude'])]
    write_json(dest/'request.json',dict(spec=spec,physical_inputs=p,command=command,binary_sha256=sha(binary)))
    print('Flat channel starts',spec['id'],flush=True);started=perf_counter()
    with (dest/'stdout.log').open('w') as stream:
        result=subprocess.run(command,cwd=dest,stdout=stream,stderr=subprocess.STDOUT,timeout=3600)
    write_json(dest/'execution.json',dict(returncode=result.returncode,wall_s=perf_counter()-started))
    if result.returncode:raise RuntimeError('Flat channel native failed: '+str(dest/'stdout.log'))
    assessment=analyze(dest);write_json(dest/'analysis.json',assessment)
    print('Flat channel done',spec['id'],f'{perf_counter()-started:.1f}s',flush=True)
    return assessment


def run():
    rid=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid4().hex[:8]
    base=OUTPUT/'runs'/rid;base.mkdir(parents=True)
    files=[HERE/'run.py',HERE/'reference.py',HERE/'protocol.md',*list((HERE/'openlb').iterdir()),
           HERE.parent.parent/'channel/reference.py']
    sources=source_manifest(files)
    for name in sources:
        dest=base/'source_snapshot'/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/name,dest)
    builddir=ROOT/'.local/adsorption/transport_root_cause_flat_channel'
    build=build_case(HERE/'openlb',builddir,base/'build');specs=specifications();cases={}
    for index,spec in enumerate(specs):
        cases[spec['id']]=execute(base,spec,builddir/'capillary3d')
        if index==1:
            agreement=equivalence(base/'native/full32',base/'native/thin32')
            write_json(base/'thin_equivalence.json',agreement)
            if not agreement['passed']:raise ValueError('Flat-channel thin reduction failed; no higher grids run')
    for name,digest in sources.items():
        if sha(ROOT/name)!=digest:raise ValueError('Flat channel source changed during run: '+name)
    record=dict(schema='pore2process.flat-wall-poiseuille-diagnostic.v1',run_id=rid,
                source_sha256=sources,build=build,specifications=specs,physical_inputs=PHYSICAL,
                cases=cases,thin_equivalence=agreement,vendor=verify_vendor(),new_native_runs=len(specs),
                new_DWSIM_runs=0,application_transport_qualified=False)
    record['artifacts_sha256']={str(p.relative_to(base)):sha(p) for p in base.rglob('*') if p.is_file()}
    write_json(base/'summary.json',record);write_json(OUTPUT/'latest.json',dict(summary=str((base/'summary.json').relative_to(ROOT))))
    print('Flat channel complete',base,flush=True);return base/'summary.json'


def check():
    p=ROOT/json.loads((OUTPUT/'latest.json').read_text())['summary'];record=json.loads(p.read_text());base=p.parent
    for name,digest in record['source_sha256'].items():
        if sha(ROOT/name)!=digest or sha(base/'source_snapshot'/name)!=digest:raise ValueError('Flat channel source changed: '+name)
    for name,digest in record['artifacts_sha256'].items():
        if sha(base/name)!=digest:raise ValueError('Flat channel artifact changed: '+name)
    for name,expected in record['cases'].items():
        if analyze(base/'native'/name)!=expected:raise ValueError('Flat channel analysis changed: '+name)
    if equivalence(base/'native/full32',base/'native/thin32')!=record['thin_equivalence']:raise ValueError('Flat thin equivalence changed')
    result=dict(run_id=record['run_id'],source_summary=str(p.relative_to(ROOT)),source_summary_sha256=sha(p),
                sources_checked=len(record['source_sha256']),artifacts_checked=len(record['artifacts_sha256']),
                actual_native_runs=record['new_native_runs'],new_native_runs_in_check=0,new_DWSIM_runs_in_check=0,
                thin_equivalence=record['thin_equivalence'],vendor=verify_vendor())
    write_json(OUTPUT/'verification.json',result);print(json.dumps(result,indent=2));return result


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--check',action='store_true');args=parser.parse_args()
    check() if args.check else run()


if __name__=='__main__':main()
