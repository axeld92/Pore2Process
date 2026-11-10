"""First synthetic adsorption loop: real native K -> Julia -> actual DWSIM pump."""
from pore2process_demo.provenance import implementation_files, source_manifest
import csv
from datetime import datetime,timezone
import json
import math
from pathlib import Path
import shutil
from time import perf_counter
from uuid import uuid4
import numpy as np
from pore2process_demo.runtime.openlb import ROOT, verify_vendor
from pore2process_demo.io import write_json
from pore2process_demo.io import sha
from .dwsim import build_adapter,solve_case
from pore2process_demo.bridge import runtime,runtime_versions
from .bridge import module,column
from .evidence import evidence

def read(path):return json.loads(Path(path).read_text())
def table(path,rows):
    with path.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)

def execute():
    native=ROOT/'results/adsorption/bed/summary.json';transport=ROOT/'results/adsorption/transport/summary.json'
    n=read(native);tr=read(transport)
    if n['status']!='passed' or tr['status']!='passed':raise ValueError('Complete both native verification stages first.')
    from .bed import analyze
    for name,c in n['cases'].items():
        actual=analyze(native.parent/name)
        if actual!=c:raise ValueError('Native field differs from accepted analysis.')
    e=evidence();start=perf_counter();jl=runtime();ads=module();initialization=perf_counter()-start
    reported={v['name']:v['SI_value'] for v in e['measurements_and_fitted_parameters']}
    base=ROOT/'results/adsorption/process/runs'/(datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid4().hex[:8])
    base.mkdir(parents=True)
    sources=list((ROOT/'examples/adsorption').glob('*.py'))+list((ROOT/'examples/adsorption/openlb').glob('*'))+list((ROOT/'examples/adsorption/openlb_bed').glob('*'))+list((ROOT/'examples/adsorption/dwsim_runner').glob('*'))+[
        ROOT/'julia/adsorption/AdsorptionColumn.jl',ROOT/'uv.lock',ROOT/'julia/Project.toml',ROOT/'julia/Manifest.toml',
        ROOT/'src/pore2process_demo/runtime/dwsim_runner/Program.cs',ROOT/'src/pore2process_demo/runtime/dwsim.py',ROOT/'src/pore2process_demo/bridge.py']
    sources = implementation_files(sources)
    state={'evidence':e,'native_sha256':sha(native),'transport_sha256':sha(transport),
        'source_sha256':source_manifest(sources),'versions':runtime_versions(),
        'Julia_runtime_and_module_initialization_s':initialization}
    write_json(base/'inputs.json',state)
    for f in sources:
        if f.is_file():
            target=base/'source_snapshot'/f.relative_to(ROOT);target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(f,target)
    shutil.copy2(native,base/'native_summary.json');shutil.copy2(transport,base/'transport_summary.json')
    shutil.copytree(native.parent,base/'native_evidence')
    shutil.copytree(transport.parent,base/'transport_evidence')
    from .export import export_native
    export_native(base/'native_evidence'/n['accepted_case'],base/'native_bed_reference.vtp')
    build=build_adapter(base/'adapter_build')
    anchor_inputs={'Temperature_K':298.15,'Inlet_pressure_Pa':1e5,'Outlet_pressure_Pa':2e5,'Mass_flow_kg_s':.1,'Pump_efficiency_fraction':.75,
        'Scope':'Actual water-carrier property anchor for the new synthetic adsorption case; not an adsorption solve.'}
    anchor=solve_case(base/'water_anchor',anchor_inputs)
    mu=anchor['inlet_viscosity_Pa_s'];rho=anchor['inlet_density_kg_m3']
    # Preserve reference particle Re under the actual carrier viscosity/density.
    u=n['cases'][n['accepted_case']]['superficial_velocity_m_s']*(mu/rho)/(1e-6)
    parameters={'length_m':1.,'area_m2':.01,'porosity':.6,'particle_density_kg_m3':1000.,
        'superficial_velocity_m_s':u,'inlet_concentration_kg_m3':.05,'capacity_kg_kg':reported['Langmuir_capacity'],
        'affinity_m3_kg':reported['Langmuir_affinity'],'axial_dispersion_m2_s':0.,'breakthrough_fraction':.05,
        'initial_concentration_kg_m3':0.,'initial_loading_kg_kg':0.,'linear_partition_m3_kg':None,'equilibrium':False}
    equilibrium=float(ads.equilibrium_loading(.05,reported['Langmuir_capacity'],reported['Langmuir_affinity']))
    capacity_time=(.6+400*equilibrium/.05)/u
    cases=[]
    for dp in (.000575,.00115,.0023):
        print('Julia adsorption and real DWSIM: dp='+str(dp),flush=True)
        native_case={.000575:'smallest',.00115:'scaled',.0023:n['accepted_case']}[dp]
        K=n['cases'][native_case]['permeability_m2']
        rate=float(ads.ldf_rate(reported['effective_intraparticle_diffusivity_at_50mg_L'],dp))
        p={**parameters,'ldf_s_inv':rate,'t_end_s':2*capacity_time,'cells':640,'dt_s':capacity_time/4000}
        value=column(**p)
        if value['breakthrough_time_s'] is None or value['normalized_mass_residual']>1e-8:raise ValueError('Unresolved breakthrough or mass imbalance.')
        load=ads.pressure_load(K,mu,p['length_m'],u,p['area_m2'],.75)
        pump_inputs={'Temperature_K':298.15,'Inlet_pressure_Pa':1e5,'Outlet_pressure_Pa':1e5+float(load.pressure_drop_Pa),
            'Mass_flow_kg_s':rho*float(load.flow_m3_s),'Pump_efficiency_fraction':.75,
            'Scope':'Actual pump load from native FCC permeability, geometrically scaled within low-Re scope. Julia alone calculates dilute phenol uptake. No DWSIM adsorption unit or experimental material calibration claimed.'}
        directory=base/f'dp{dp*1e6:g}um';directory.mkdir();write_json(directory/'inputs.json',{'column':p,'particle_diameter_m':dp,'permeability_m2':K,'pump':pump_inputs})
        pump=solve_case(directory/'dwsim',pump_inputs)
        if pump['pump_power_relative_error']>1e-8 or abs(pump['mass_balance_error_kg_s'])>1e-10*pump_inputs['Mass_flow_kg_s']:
            raise ValueError('Actual pump work or mass balance failed.')
        if max(pump['inlet_vapor_fraction'],pump['outlet_vapor_fraction'])>1e-10:raise ValueError('Water carrier is not liquid.')
        if abs(pump['energy_balance_error_W'])>1e-3*pump['shaft_power_W']:raise ValueError('Actual pump enthalpy balance failed.')
        sensitivity=[]
        for factor in (.1,1.):
            alt={**p,'axial_dispersion_m2_s':factor*u/p['porosity']*dp}
            sensitivity.append({'chosen_dispersion_factor':factor,'inputs':alt,'column':column(**alt),
                'scope':'Chosen Dax = factor * interstitial velocity * particle diameter; not a measured dispersion closure or experimental uncertainty.'})
        result={'inputs':p,'particle_diameter_m':dp,'native_case':native_case,'permeability_m2':K,'column':value,'dwsim':pump,
            'dispersion_sensitivity':sensitivity,
            'pump_energy_to_breakthrough_J':pump['shaft_power_W']*value['breakthrough_time_s'],
            'specific_pump_energy_kWh_kg_adsorbed':pump['shaft_power_W']*value['breakthrough_time_s']/3.6e6/value['adsorbed_solute_at_breakthrough_kg'],
            'specific_pump_energy_kWh_m3':pump['shaft_power_W']/(u*.01)/3.6e6}
        write_json(directory/'summary.json',result)
        table(directory/'breakthrough.csv',[dict(time_s=t,outlet_concentration_kg_m3=c,cumulative_inlet_kg=mi,cumulative_outlet_kg=mo,inventory_kg=inv,mass_residual_kg=err) for t,c,mi,mo,inv,err in zip(value['time_s'],value['outlet_concentration_kg_m3'],value['cumulative_inlet_kg'],value['cumulative_outlet_kg'],value['inventory_kg'],value['mass_residual_kg'],strict=True)])
        table(directory/'final_profile.csv',[dict(x_m=x,liquid_concentration_kg_m3=c,solid_loading_kg_kg=q) for x,c,q in zip(value['x_m'],value['final_concentration_kg_m3'],value['final_loading_kg_kg'],strict=True)])
        cases.append(result)
    # Distinct spatial and temporal checks on the slowest-transfer case.
    p=cases[-1]['inputs'];fine=column(**{**p,'cells':1280});time=column(**{**p,'dt_s':p['dt_s']/2})
    refinement={'space_breakthrough_relative':abs(fine['breakthrough_time_s']/cases[-1]['column']['breakthrough_time_s']-1),
        'time_breakthrough_relative':abs(time['breakthrough_time_s']/cases[-1]['column']['breakthrough_time_s']-1),
        'fine_grid':fine,'half_time_step':time}
    if refinement['space_breakthrough_relative']>.03 or refinement['time_breakthrough_relative']>.01:raise ValueError('Column refinement failed.')
    rows=[dict(particle_diameter_m=c['particle_diameter_m'],permeability_m2=c['permeability_m2'],ldf_s_inv=c['inputs']['ldf_s_inv'],
        adsorbent_mass_kg=c['column']['adsorbent_mass_kg'],breakthrough_time_s=c['column']['breakthrough_time_s'],
        treated_volume_m3=c['column']['treated_volume_m3'],bed_pressure_drop_Pa=c['dwsim']['outlet_pressure_Pa']-c['dwsim']['inlet_pressure_Pa'],
        solid_utilization_at_breakthrough=c['column']['solid_utilization_at_breakthrough'],
        adsorbed_solute_at_breakthrough_kg=c['column']['adsorbed_solute_at_breakthrough_kg'],
        actual_pump_power_W=c['dwsim']['shaft_power_W'],specific_pump_energy_kWh_m3=c['specific_pump_energy_kWh_m3']) for c in cases]
    table(base/'comparison.csv',rows)
    table(base/'dispersion_sensitivity.csv',[dict(particle_diameter_m=c['particle_diameter_m'],chosen_factor=v['chosen_dispersion_factor'],
        axial_dispersion_m2_s=v['inputs']['axial_dispersion_m2_s'],breakthrough_time_s=v['column']['breakthrough_time_s'],
        treated_volume_m3=v['column']['treated_volume_m3'],solid_utilization_at_breakthrough=v['column']['solid_utilization_at_breakthrough'])
        for c in cases for v in c['dispersion_sensitivity']])
    from .plot import plot
    plot(cases,n,tr,base/'adsorption.png')
    s={'status':'passed_synthetic_example','experimental_validation':'not performed; source conflicts and transfer assumptions retained',**state,
        'cases':cases,'comparison':rows,'refinement':refinement,'property_anchor':anchor,'adapter':build,'original_vendor':verify_vendor(),
        'scope':'Isothermal single-solute clean-bed LDF prototype. Native bed hydraulics and free-fluid scalar transport are verified separately. Bed-scale dispersion and film transfer are not resolved native closures; no reaction, regeneration, heat model or economic optimum.'}
    s['artifacts_sha256']={str(f.relative_to(base)):sha(f) for f in base.rglob('*') if f.is_file()}
    write_json(base/'summary.json',s)
    from .check import verify
    checked=verify(base/'summary.json');write_json(base/'verification.json',checked)
    write_json(ROOT/'results/adsorption/process/latest.json',{'summary':str(base/'summary.json')})
    print(json.dumps({'run':str(base),'comparison':rows,'verification':checked},indent=2))
if __name__=='__main__':execute()
