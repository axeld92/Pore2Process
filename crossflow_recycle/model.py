"""SI contracts, real Julia closure and independent DWSIM balance checks."""
import math
from pathlib import Path
import xml.etree.ElementTree as ET
import zipfile

from calibration.check import close_tree
from graded_support.model import load_case

STREAMS={'Fresh','Fresh-high','Mixed','Tear-in','Circulation-in','Circulation-out','Membrane-feed',
    'Permeate-high','Permeate','Module-retentate','Purge-high','Purge','Return'}
BALANCES={
    'feed_pump':(['Fresh'],['Fresh-high']),
    'mixer':(['Fresh-high','Return'],['Mixed']),
    'cooler':(['Mixed'],['Tear-in']),
    'recycle':(['Tear-in'],['Circulation-in']),
    'circulation_pump':(['Circulation-in'],['Circulation-out']),
    'channel':(['Circulation-out'],['Membrane-feed']),
    'membrane_splitter':(['Membrane-feed'],['Permeate-high','Module-retentate']),
    'membrane_valve':(['Permeate-high'],['Permeate']),
    'purge_splitter':(['Module-retentate'],['Purge-high','Return']),
    'purge_valve':(['Purge-high'],['Purge']),
    'overall':(['Fresh'],['Permeate','Purge']),
}
CONNECTIONS={
    ('Fresh','0','P-feed','0'),('P-feed','0','Fresh-high','0'),('Fresh-high','0','M-1','0'),
    ('Return','0','M-1','1'),('M-1','0','Mixed','0'),('Mixed','0','C-1','0'),
    ('C-1','0','Tear-in','0'),('C-1','energy','Cooling-duty','0'),('Tear-in','0','REC-1','0'),
    ('REC-1','0','Circulation-in','0'),('Circulation-in','0','P-circulation','0'),
    ('P-circulation','0','Circulation-out','0'),('Circulation-out','0','V-channel','0'),
    ('V-channel','0','Membrane-feed','0'),('Membrane-feed','0','S-membrane','0'),
    ('S-membrane','0','Permeate-high','0'),('S-membrane','1','Module-retentate','0'),
    ('Permeate-high','0','V-membrane','0'),('V-membrane','0','Permeate','0'),
    ('Module-retentate','0','S-purge','0'),('S-purge','0','Purge-high','0'),('S-purge','1','Return','0'),
    ('Purge-high','0','V-purge','0'),('V-purge','0','Purge','0'),
    ('Feed-power','0','P-feed','1'),('Circulation-power','0','P-circulation','1'),
}

def translate(jl,source,support,p,y,C,resistance_multiplier=1.):
    fixed=source['literature_reference']['fit']['fixed_impedance_Pa_s_m']
    membrane=load_case(jl,support['body_resistance_Pa_s_m3'],source['supports'][0]['body_resistance_Pa_s_m3'],source['protocol'],fixed)
    v=jl.crossflow_recycle_load(membrane['pressure_drop_Pa'],membrane['flow_m3_s'],y,C,
        p['nominal_channel_resistance_Pa_s_m3']*resistance_multiplier,
        p['feed_pump_efficiency_fraction'],p['circulation_pump_efficiency_fraction'])
    return membrane,{str(k):float(getattr(v,str(k))) for k in jl.fieldnames(jl.typeof(v))}

def adapter_inputs(p,source,c,*,seed_multiplier=None,seed_temperature_offset=None,tolerance_multiplier=1.):
    r=p['recycle'];rho=source['fluid_anchor']['inlet_density_kg_m3']
    seed_multiplier=r['nominal_seed_flow_multiplier'] if seed_multiplier is None else seed_multiplier
    offset=r['nominal_seed_temperature_offset_K'] if seed_temperature_offset is None else seed_temperature_offset
    return {'Temperature_K':p['temperature_K'],'External_pressure_Pa':p['external_absolute_pressure_Pa'],
        'TMP_Pa':c['transmembrane_pressure_Pa'],'Channel_drop_Pa':c['channel_pressure_drop_Pa'],
        'Fresh_mass_flow_kg_s':rho*c['fresh_feed_reference_flow_m3_s'],
        'Module_to_permeate_ratio':c['module_to_permeate_ratio'],'Recovery_fraction':c['overall_recovery_fraction'],
        'Feed_efficiency_fraction':p['feed_pump_efficiency_fraction'],
        'Circulation_efficiency_fraction':p['circulation_pump_efficiency_fraction'],
        'Seed_module_mass_flow_kg_s':rho*c['module_reference_flow_m3_s']*seed_multiplier,
        'Seed_temperature_K':p['temperature_K']+offset,
        'Mass_tolerance_kg_s':r['mass_tolerance_kg_s']*tolerance_multiplier,
        'Temperature_tolerance_K':r['temperature_tolerance_K']*tolerance_multiplier,
        'Pressure_tolerance_Pa':r['pressure_tolerance_Pa']*tolerance_multiplier,
        'Maximum_iterations':r['maximum_iterations'],
        'Scope':'Pure-water hydraulic recycle: stock mixer, cooler, two pumps, splitters, pressure-loss valves and recycle; Julia supplies TMP and reference channel law. Cooling heat is separate from shaft energy.'}

def topology(path):
    with zipfile.ZipFile(path) as z:
        names=z.namelist()
        if len(names)!=1 or not names[0].endswith('.xml'):raise RuntimeError('Changed DWSIM archive contract')
        root=ET.fromstring(z.read(names[0]))
    graphics=root.findall('./GraphicObjects/GraphicObject');tags={g.findtext('Name'):g.findtext('Tag') for g in graphics}
    actual=set()
    for g in graphics:
        for group in ('OutputConnectors','EnergyConnector'):
            for i,c in enumerate(g.findall('./'+group+'/Connector')):
                if c.get('IsAttached')=='true' and c.get('AttachedToObjID'):
                    actual.add((g.findtext('Tag'),str(i) if group=='OutputConnectors' else 'energy',
                        tags[c.get('AttachedToObjID')],c.get('AttachedToConnIndex')))
    if actual!=CONNECTIONS:raise RuntimeError('Saved recycle topology differs from contract')
    # Confirm that the saved numerical tear setting survives the archive.
    sim=root.findall('./SimulationObjects/SimulationObject')
    material=[o for o in sim if o.findtext('Name') in {g.findtext('Name') for g in graphics if g.findtext('Tag') in ('Tear-in','Circulation-in')}]
    if len(material)!=2 or any(o.findtext('OverrideSingleCompoundFlashBehavior')!='true' for o in material):
        raise RuntimeError('Saved tear TP setting differs')
    return sorted([list(v) for v in actual])

def balances(s,powers):
    def total(names,key):return sum(s[n]['mass_flow_kg_s']*(s[n]['enthalpy_J_kg'] if key=='enthalpy' else 1.) for n in names)
    mass={};energy={}
    for name,(ins,outs) in BALANCES.items():
        mass[name]=total(outs,'mass')-total(ins,'mass')
        work=powers['feed_W'] if name=='feed_pump' else powers['circulation_W'] if name=='circulation_pump' else powers['feed_W']+powers['circulation_W'] if name=='overall' else 0.
        cooling=powers['cooling_removed_W'] if name in ('cooler','overall') else 0.
        energy[name]=total(outs,'enthalpy')-total(ins,'enthalpy')-work+cooling
    return mass,energy

def validate(actual,inputs,closure,p,source):
    close_tree(actual['inputs'],inputs,'DWSIM inputs')
    if actual['tool']!='DWSIM' or actual['property_package']!=source['fluid_anchor']['property_package'] or actual['compounds']!=['Water']:
        raise RuntimeError('Changed simulator/property contract')
    s=actual['streams'];lim=p['acceptance'];power=actual['powers']
    if set(s)!=STREAMS or set(actual['reloaded_streams'])!=STREAMS:raise RuntimeError('Missing/extra material stream')
    mf=inputs['Fresh_mass_flow_kg_s'];mp=mf*inputs['Recovery_fraction'];mm=mp*inputs['Module_to_permeate_ratio']
    expected={n:mm for n in STREAMS}
    expected.update({n:mf for n in ('Fresh','Fresh-high')});expected.update({n:mp for n in ('Permeate','Permeate-high')})
    expected['Module-retentate']=mm-mp
    expected.update({n:mf-mp for n in ('Purge','Purge-high')});expected['Return']=mm-mf
    rho=source['fluid_anchor']['inlet_density_kg_m3'];mu=source['fluid_anchor']['inlet_viscosity_Pa_s']
    close_tree(s['Fresh']['density_kg_m3'],rho,'anchor density');close_tree(s['Fresh']['viscosity_Pa_s'],mu,'anchor viscosity')
    max_mu=max_rho=0.
    for n,state in s.items():
        if not all(math.isfinite(v) for v in state.values()) or min(state['mass_flow_kg_s'],state['density_kg_m3'],state['viscosity_Pa_s'],state['temperature_K'])<=0:
            raise RuntimeError('Invalid stream state')
        if state['vapor_fraction']!=0 or state['water_mass_fraction']!=1:raise RuntimeError('Require pure liquid water')
        target_p=inputs['External_pressure_Pa']+(0. if n in ('Fresh','Permeate','Purge') else inputs['TMP_Pa'])
        if n=='Circulation-out':target_p+=inputs['Channel_drop_Pa']
        if abs(state['pressure_Pa']-target_p)>lim['pressure_absolute_Pa']:raise RuntimeError('Pressure profile differs')
        if not math.isclose(state['mass_flow_kg_s'],expected[n],rel_tol=lim['stream_mass_relative']):raise RuntimeError('Flow/recovery differs')
        close_tree(state['actual_volume_flow_m3_s'],state['mass_flow_kg_s']/state['density_kg_m3'],'actual volume')
        after=actual['reloaded_streams'][n]
        if not all(math.isfinite(v) for v in after.values()):raise RuntimeError('Invalid reloaded state')
        for key,tol in (('pressure_Pa',lim['pressure_absolute_Pa']),('mass_flow_kg_s',lim['reload_relative_mass']*expected[n]),
                        ('temperature_K',lim['reload_temperature_K']),('enthalpy_J_kg',lim['reload_enthalpy_J_kg'])):
            if abs(after[key]-state[key])>tol:raise RuntimeError('Reloaded state changed')
        max_mu=max(max_mu,abs(state['viscosity_Pa_s']/mu-1));max_rho=max(max_rho,abs(state['density_kg_m3']/rho-1))
    for n in ('Fresh','Tear-in','Circulation-in'):
        if abs(s[n]['temperature_K']-p['temperature_K'])>1e-8:raise RuntimeError('Temperature control differs')
    for label in ('recycle','reloaded_recycle'):
        r=actual[label]
        if r['converged'] is not True or not (1<=r['iterations_taken']<=inputs['Maximum_iterations']):raise RuntimeError('Recycle did not converge')
        if r['acceleration']!=p['recycle']['acceleration'] or r['legacy_mode']!=p['recycle']['legacy_mode']:raise RuntimeError('Changed recycle method')
        for key,requested in (('maximum_iterations','Maximum_iterations'),('temperature_tolerance_K','Temperature_tolerance_K'),
                              ('pressure_tolerance_Pa','Pressure_tolerance_Pa'),('mass_tolerance_kg_s','Mass_tolerance_kg_s')):
            close_tree(r[key],inputs[requested],'recycle settings')
    tear={k:abs(s['Tear-in'][k]-s['Circulation-in'][k]) for k in ('temperature_K','pressure_Pa','mass_flow_kg_s')}
    total=power['feed_W']+power['circulation_W'];mass,energy=balances(s,power)
    if not all(math.isfinite(v) and v>0 for v in power.values()):raise RuntimeError('Invalid pump/cooling duty')
    corrected={}
    for pump,inlet,drop,eta in (('feed','Fresh',inputs['TMP_Pa'],inputs['Feed_efficiency_fraction']),
                               ('circulation','Circulation-in',inputs['Channel_drop_Pa'],inputs['Circulation_efficiency_fraction'])):
        st=s[inlet];oracle=st['mass_flow_kg_s']/st['density_kg_m3']*drop/eta
        corrected[pump+'_W']=oracle
        if not math.isclose(power[pump+'_W'],oracle,rel_tol=lim['pump_property_corrected_relative_error']):raise RuntimeError('Property-corrected pump oracle failed')
        if not math.isclose(power[pump+'_energy_stream_W'],power[pump+'_W'],rel_tol=1e-11):raise RuntimeError('Pump energy stream differs')
    close_tree(power['cooling_removed_W'],power['cooling_energy_stream_W'],'cooler energy stream')
    for name,key in (('Feed-power','feed_W'),('Circulation-power','circulation_W'),('Cooling-duty','cooling_removed_W')):
        if not math.isclose(actual['reloaded_energy_streams_W'][name],power[key],rel_tol=lim['seed_invariance_relative_power']):raise RuntimeError('Reloaded duty differs')
    checks={
        'all_unit_and_overall_mass':max(map(abs,mass.values()))<=lim['unit_mass_relative']*mf,
        'all_unit_and_overall_enthalpy':max(map(abs,energy.values()))<=lim['unit_and_overall_energy_relative_to_total_shaft']*total,
        'tear_mass':tear['mass_flow_kg_s']<=inputs['Mass_tolerance_kg_s']*1.01+1e-15,
        'tear_temperature':tear['temperature_K']<=inputs['Temperature_tolerance_K']*1.01+1e-12,
        'tear_pressure':tear['pressure_Pa']<=inputs['Pressure_tolerance_Pa']*1.01+1e-9,
        'reference_density_power_screen':abs(total/closure['total_shaft_power_W']-1)<=lim['reference_density_power_relative_error'],
        'constant_viscosity_screen':max_mu<=lim['relative_viscosity_change'],
        'constant_density_screen':max_rho<=lim['relative_density_change'],
        'tear_TP_configuration':actual['tear_single_compound_flash_override'] is True,
        'save_reload_solve':actual['save_reload_solve_checked'] is True,
    }
    expected_units={**{n:'DWSIM.Thermodynamics.Streams.MaterialStream' for n in STREAMS},
        **{n:'DWSIM.UnitOperations.Streams.EnergyStream' for n in ('Feed-power','Circulation-power','Cooling-duty')},
        **{n:'DWSIM.UnitOperations.UnitOperations.Pump' for n in ('P-feed','P-circulation')},
        **{n:'DWSIM.UnitOperations.UnitOperations.Valve' for n in ('V-channel','V-membrane','V-purge')},
        **{n:'DWSIM.UnitOperations.UnitOperations.Splitter' for n in ('S-membrane','S-purge')},
        'M-1':'DWSIM.UnitOperations.UnitOperations.Mixer','C-1':'DWSIM.UnitOperations.UnitOperations.Cooler',
        'REC-1':'DWSIM.UnitOperations.SpecialOps.Recycle'}
    close_tree(actual['units'],expected_units,'stock unit types')
    if inputs.get('Force_liquid',False):
        expected_phases={name:'Liquid' for name in STREAMS}
        close_tree(actual['forced_phases'],expected_phases,'forced liquid flash configuration')
        close_tree(actual['reloaded_forced_phases'],expected_phases,'reloaded liquid flash configuration')
        for name,stream in s.items():
            psat=actual['saturation_pressure_Pa'][name]
            if not math.isfinite(psat) or psat<=0 or stream['temperature_K']<=273.15 or stream['pressure_Pa']<=2*psat:
                raise RuntimeError('Forced-liquid phase assumption is outside its saturation/freezing guard.')
    if not all(checks.values()):raise RuntimeError(f'Recycle acceptance failed: {checks}')
    return {'checks':checks,'mass_residuals_kg_s':mass,'enthalpy_residuals_W':energy,'tear_residuals':tear,
        'max_mass_relative_residual':max(map(abs,mass.values()))/mf,
        'max_enthalpy_relative_residual':max(map(abs,energy.values()))/total,
        'max_relative_viscosity_change':max_mu,'max_relative_density_change':max_rho,
        'property_corrected_reference_powers':corrected,'total_shaft_power_W':total,
        'specific_shaft_energy_kWh_m3_reference_volume':total/(s['Permeate']['mass_flow_kg_s']/rho)/3.6e6,
        'circulation_shaft_fraction':power['circulation_W']/total,
        'reference_density_power_relative_difference':total/closure['total_shaft_power_W']-1}
