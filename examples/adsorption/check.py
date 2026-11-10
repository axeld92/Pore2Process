"""Replay SI/native/Julia/process contracts; no new external solves."""
from pore2process_demo.provenance import verified_source_path
import csv
import json
import math
from pathlib import Path
import zipfile
from xml.etree import ElementTree as ET
from examples.filtration.calibration.check import close_tree
from pore2process_demo.runtime.openlb import ROOT, verify_vendor
from pore2process_demo.io import write_json
from pore2process_demo.io import sha
from .bridge import column,module
from .evidence import evidence

def read(p):return json.loads(Path(p).read_text())

def topology(path):
    with zipfile.ZipFile(path) as z:
        names=z.namelist()
        if len(names)!=1:raise ValueError('Unexpected pump archive.')
        root=ET.fromstring(z.read(names[0]))
    objs=root.findall('./GraphicObjects/GraphicObject')
    if not objs:objs=list(root.find('GraphicObjects'))
    by_id={o.findtext('Name'):o.findtext('Tag') for o in objs}
    found=set()
    for o in objs:
        for index,c in enumerate(o.findall('./OutputConnectors/Connector')):
            if c.get('IsAttached')=='true':found.add((o.findtext('Tag'),index,by_id[c.get('AttachedToObjID')],int(c.get('AttachedToConnIndex'))))
    expected={('Water-in',0,'P-1',0),('Pump-power',0,'P-1',1),('P-1',0,'Water-out',0)}
    if found!=expected:raise ValueError('Saved pump connector topology differs.')
    return sorted(found)

def artifacts(path,s):
    base=path.parent.resolve()
    actual={str(f.relative_to(base)) for f in base.rglob('*') if f.is_file() and f.resolve()!=path.resolve() and f.name!='verification.json'}
    if actual!=set(s['artifacts_sha256']):raise ValueError('Artifact coverage changed.')
    for name,digest in s['artifacts_sha256'].items():
        f=(base/name).resolve()
        if base not in f.parents or sha(f)!=digest:raise ValueError('Changed artifact: '+name)

def verify(path):
    path=Path(path).resolve();base=path.parent;s=read(path)
    if s['status']!='passed_synthetic_example':raise ValueError('Incomplete adsorption example.')
    artifacts(path,s)
    for name,digest in s['source_sha256'].items():
        if sha(verified_source_path(name,digest))!=digest:raise ValueError('Implementation changed: '+name)
        if sha(base/'source_snapshot'/name)!=digest:raise ValueError('Source snapshot changed.')
    for key,value in read(base/'inputs.json').items():close_tree(s[key],value,'Run input snapshot.'+key)
    close_tree({k:v for k,v in evidence().items() if k!='sources'}, {k:v for k,v in s['evidence'].items() if k!='sources'},'Retained measured facts and assumptions')
    n=read(base/'native_summary.json');tr=read(base/'transport_summary.json')
    if sha(base/'native_summary.json')!=s['native_sha256'] or sha(base/'transport_summary.json')!=s['transport_sha256']:raise ValueError('Native summary changed.')
    from .bed import analyze
    for name,value in n['cases'].items():
        close_tree(analyze(base/'native_evidence'/name),value,'Native bed field')
    from .transport import analyze as scalar
    for name,value in tr['cases'].items():close_tree(scalar(base/'transport_evidence'/name),value,'Native scalar field')
    if not all(n['checks'].values()) or not all(tr['checks'].values()):raise ValueError('Unaccepted native check.')
    module_ = module();topologies=[topology(base/'water_anchor/water_pump.dwxmz')]
    anchor=read(base/'water_anchor/summary.json');close_tree(anchor,s['property_anchor'],'Actual water anchor')
    close_tree(anchor['inputs'],read(base/'water_anchor/inputs.json'),'Actual water-anchor request')
    if not anchor['forced_liquid'] or not anchor['save_reload_solve_checked'] or anchor['pump_power_relative_error']>1e-8 or abs(anchor['energy_balance_error_W'])>1e-3*anchor['shaft_power_W']:raise ValueError('Water property-anchor solve failed.')
    replayed=0
    for c in s['cases']:
        directory=base/f"dp{c['particle_diameter_m']*1e6:g}um"
        replay=column(**c['inputs']);replay.pop('total_s');original={k:v for k,v in c['column'].items() if k!='total_s'}
        close_tree(replay,original,'Julia conservative column replay')
        replayed+=1
        if c['column']['normalized_mass_residual']>1e-8:raise ValueError('Column mass balance failed.')
        close_tree(c['permeability_m2'],n['cases'][c['native_case']]['permeability_m2'],'Actual measured native K')
        a=c['dwsim']['inputs'];actual=read(directory/'dwsim/summary.json');close_tree(actual,c['dwsim'],'Actual DWSIM output')
        close_tree(a,read(directory/'dwsim/inputs.json'),'Actual pump request')
        p=c['inputs'];load=module_.pressure_load(c['permeability_m2'],s['property_anchor']['inlet_viscosity_Pa_s'],p['length_m'],p['superficial_velocity_m_s'],p['area_m2'],a['Pump_efficiency_fraction'])
        close_tree(a['Outlet_pressure_Pa']-a['Inlet_pressure_Pa'],float(load.pressure_drop_Pa),'Native→Julia→DWSIM pressure')
        expected=a['Mass_flow_kg_s']*(a['Outlet_pressure_Pa']-a['Inlet_pressure_Pa'])/(actual['inlet_density_kg_m3']*a['Pump_efficiency_fraction'])
        close_tree(expected,actual['shaft_power_W'],'Actual pump work')
        if not actual['save_reload_solve_checked'] or actual['pump_power_relative_error']>1e-8 or abs(actual['energy_balance_error_W'])>1e-3*actual['shaft_power_W']:raise ValueError('Actual DWSIM pump checks failed.')
        if not actual['forced_liquid'] or min(actual['inlet_pressure_Pa'],actual['outlet_pressure_Pa'])<=2*actual['saturation_pressure_Pa']:raise ValueError('Liquid phase assumption failed.')
        close_tree(actual['shaft_power_W']*replay['breakthrough_time_s'],c['pump_energy_to_breakthrough_J'],'Cycle pumping energy')
        close_tree(actual['shaft_power_W']/(p['superficial_velocity_m_s']*p['area_m2'])/3.6e6,c['specific_pump_energy_kWh_m3'],'kWh per reference treated m3')
        close_tree(c['pump_energy_to_breakthrough_J']/3.6e6/replay['adsorbed_solute_at_breakthrough_kg'],c['specific_pump_energy_kWh_kg_adsorbed'],'kWh per kg adsorbed')
        topologies.append(topology(directory/'dwsim/water_pump.dwxmz'))
        with (directory/'breakthrough.csv').open(newline='') as f:rows=list(csv.DictReader(f))
        if len(rows)!=len(replay['time_s']):raise ValueError('Exported time grid changed.')
        for key,original_key in [('time_s','time_s'),('outlet_concentration_kg_m3','outlet_concentration_kg_m3'),('mass_residual_kg','mass_residual_kg')]:
            close_tree([float(r[key]) for r in rows],replay[original_key],'SI export.'+key)
        with (directory/'final_profile.csv').open(newline='') as f:rows=list(csv.DictReader(f))
        for key,original_key in [('x_m','x_m'),('liquid_concentration_kg_m3','final_concentration_kg_m3'),('solid_loading_kg_kg','final_loading_kg_kg')]:
            close_tree([float(r[key]) for r in rows],replay[original_key],'Final SI profile.'+key)
        if abs(actual['mass_balance_error_kg_s'])>1e-10*a['Mass_flow_kg_s'] or max(actual['inlet_vapor_fraction'],actual['outlet_vapor_fraction'])>1e-10:raise ValueError('Carrier mass/phase check failed.')
        if abs(actual['inlet_viscosity_Pa_s']/s['property_anchor']['inlet_viscosity_Pa_s']-1)>.001 or abs(actual['outlet_viscosity_Pa_s']/actual['inlet_viscosity_Pa_s']-1)>.001:raise ValueError('Isothermal hydraulic viscosity screen failed.')
        for alternative in c['dispersion_sensitivity']:
            a=alternative['inputs'];close_tree(a['axial_dispersion_m2_s'],alternative['chosen_dispersion_factor']*a['superficial_velocity_m_s']/a['porosity']*c['particle_diameter_m'],'Declared dispersion scenario')
            r=column(**a);r.pop('total_s')
            close_tree(r,{k:v for k,v in alternative['column'].items() if k!='total_s'},'Dispersion sensitivity replay')
            if r['normalized_mass_residual']>1e-8:raise ValueError('Dispersion mass balance failed.')
            replayed+=1
    with (base/'comparison.csv').open(newline='') as f:rows=list(csv.DictReader(f))
    if len(rows)!=len(s['comparison']):raise ValueError('Comparison shape changed.')
    for row,original in zip(rows,s['comparison'],strict=True):
        close_tree({k:float(v) for k,v in row.items()},original,'Comparison SI values')
    for key,changes in [('fine_grid',{'cells':1280}),('half_time_step',{'dt_s':s['cases'][-1]['inputs']['dt_s']/2})]:
        replay=column(**{**s['cases'][-1]['inputs'],**changes});replay.pop('total_s')
        close_tree(replay,{k:v for k,v in s['refinement'][key].items() if k!='total_s'},'Refinement replay')
    if s['refinement']['space_breakthrough_relative']>.03 or s['refinement']['time_breakthrough_relative']>.01:raise ValueError('Column convergence failed.')
    return {'status':'passed','native_bed_cases':len(n['cases']),'native_transport_cases':len(tr['cases']),
        'actual_DWSIM_cases':len(s['cases'])+1,'artifacts_checked':len(s['artifacts_sha256']),
        'Julia_columns_replayed':replayed+2,'saved_pump_topologies':len(topologies),
        'experimental_column_validation_passed':False,'external_solves_by_checker':False,'original_vendor':verify_vendor()}

if __name__=='__main__':
    pointer=ROOT/'results/adsorption/process/latest.json';p=Path(read(pointer)['summary'])
    result=verify(p);write_json(p.parent/'verification.json',result);print(json.dumps(result,indent=2))
