"""Call Julia closures and check the actual DWSIM streams independently."""
import math
from pathlib import Path
import xml.etree.ElementTree as ET
import zipfile

from calibration.check import close_tree
from graded_support.model import load_case

STREAMS={"Feed","Pump-discharge","Membrane-feed","Permeate-high","Retentate-high","Permeate","Retentate"}
BALANCES={
    "pump":(["Feed"],["Pump-discharge"]),
    "channel":(["Pump-discharge"],["Membrane-feed"]),
    "splitter":(["Membrane-feed"],["Permeate-high","Retentate-high"]),
    "membrane":(["Permeate-high"],["Permeate"]),
    "retentate_valve":(["Retentate-high"],["Retentate"]),
    "overall":(["Feed"],["Permeate","Retentate"]),
}
CONNECTIONS={
    ("Feed",0,"P-1",0),("P-1",0,"Pump-discharge",0),("Pump-power",0,"P-1",1),
    ("Pump-discharge",0,"V-channel",0),("V-channel",0,"Membrane-feed",0),
    ("Membrane-feed",0,"S-1",0),("S-1",0,"Permeate-high",0),("S-1",1,"Retentate-high",0),
    ("Permeate-high",0,"V-membrane",0),("V-membrane",0,"Permeate",0),
    ("Retentate-high",0,"V-retentate",0),("V-retentate",0,"Retentate",0),
}

def translate(jl,source,support,recovery,channel):
    p=source["protocol"]
    Z=source["literature_reference"]["fit"]["fixed_impedance_Pa_s_m"]
    load=load_case(jl,support["body_resistance_Pa_s_m3"],source["supports"][0]["body_resistance_Pa_s_m3"],p,Z)
    v=jl.once_through_filtration(load["pressure_drop_Pa"],load["flow_m3_s"],recovery,
                                 p["process"]["pump_efficiency_fraction"],channel)
    return load,{str(k):float(getattr(v,str(k))) for k in jl.fieldnames(jl.typeof(v))}

def adapter_inputs(protocol,source,closure):
    return {"Temperature_K":protocol["temperature_K"],
        "External_pressure_Pa":protocol["external_absolute_pressure_Pa"],
        "TMP_Pa":closure["transmembrane_pressure_Pa"],
        "Channel_drop_Pa":closure["channel_pressure_drop_Pa"],
        "Feed_mass_flow_kg_s":source["fluid_anchor"]["inlet_density_kg_m3"]*closure["feed_reference_flow_m3_s"],
        "Recovery_fraction":closure["recovery_fraction"],
        "Pump_efficiency_fraction":protocol["pump_efficiency_fraction"],
        "Scope":"Once-through pure-water hydraulic surrogate: stock pump/splitter/pressure-loss valves; Julia supplies TMP. No separation selectivity, recycle or pressure recovery."}

def topology(path):
    """Inspect actual saved DWSIM connector IDs, including splitter outlet order."""
    with zipfile.ZipFile(path) as z:
        names=z.namelist()
        if len(names)!=1 or not names[0].endswith('.xml'): raise RuntimeError("Changed DWSIM archive contract")
        root=ET.fromstring(z.read(names[0]))
    objects=root.findall('./GraphicObjects/GraphicObject')
    tags={v.findtext('Name'):v.findtext('Tag') for v in objects}
    connections=set()
    for obj in objects:
        for i,c in enumerate(obj.findall('./OutputConnectors/Connector')):
            if c.get('IsAttached')=='true':
                connections.add((obj.findtext('Tag'),i,tags[c.get('AttachedToObjID')],int(c.get('AttachedToConnIndex'))))
    if connections!=CONNECTIONS: raise RuntimeError("Saved DWSIM topology differs from hydraulic contract")
    return sorted([list(v) for v in connections])

def validate(actual,inputs,closure,protocol,source):
    """Recompute balances from stream states, not stored adapter residual flags."""
    close_tree(actual['inputs'],inputs,'DWSIM snapshot')
    s=actual['streams'];lim=protocol['acceptance'];power=actual['shaft_power_W']
    if set(s)!=STREAMS or set(actual['reloaded_streams'])!=STREAMS: raise RuntimeError("Missing/extra streams")
    if actual['tool']!='DWSIM' or actual['property_package']!=source['fluid_anchor']['property_package'] or actual['compounds']!=['Water']:
        raise RuntimeError("Changed simulator/property contract")
    rho=s['Feed']['density_kg_m3'];mu=s['Feed']['viscosity_Pa_s'];mf=inputs['Feed_mass_flow_kg_s']
    if abs(s['Feed']['temperature_K']-inputs['Temperature_K'])>1e-8:
        raise RuntimeError('Fresh-feed temperature differs from request')
    if not all(math.isclose(a,b,rel_tol=1e-12) for a,b in (
        (rho,source['fluid_anchor']['inlet_density_kg_m3']),(mu,source['fluid_anchor']['inlet_viscosity_Pa_s']))):
        raise RuntimeError("Fluid anchor changed; do not reuse OpenLB properties silently")
    expected_p={n:inputs['External_pressure_Pa'] for n in STREAMS}
    expected_p.update({n:inputs['External_pressure_Pa']+inputs['TMP_Pa'] for n in ('Membrane-feed','Permeate-high','Retentate-high')})
    expected_p['Pump-discharge']+=inputs['TMP_Pa']+inputs['Channel_drop_Pa']
    expected_m={n:mf for n in STREAMS}
    expected_m.update({n:mf*inputs['Recovery_fraction'] for n in ('Permeate-high','Permeate')})
    expected_m.update({n:mf*(1-inputs['Recovery_fraction']) for n in ('Retentate-high','Retentate')})
    max_mu=max_rho=0.
    for name,state in s.items():
        if not all(math.isfinite(v) for v in state.values()): raise RuntimeError("Non-finite stream")
        if min(state['temperature_K'],state['density_kg_m3'],state['viscosity_Pa_s'],state['mass_flow_kg_s'])<=0:
            raise RuntimeError("Nonpositive water state")
        if state['vapor_fraction']!=0 or state['water_mass_fraction']!=1: raise RuntimeError("Case must stay pure liquid water")
        if abs(state['pressure_Pa']-expected_p[name])>lim['pressure_absolute_Pa']: raise RuntimeError("Stream pressure differs")
        if abs(state['mass_flow_kg_s']-expected_m[name])>lim['mass_balance_relative']*mf: raise RuntimeError("Stream flow/recovery differs")
        close_tree(state['actual_volume_flow_m3_s'],state['mass_flow_kg_s']/state['density_kg_m3'],'actual volume')
        after=actual['reloaded_streams'][name]
        for field,tol in (('pressure_Pa',lim['pressure_absolute_Pa']),('temperature_K',lim['reload_temperature_absolute_K']),
                          ('enthalpy_J_kg',lim['reload_enthalpy_absolute_J_kg']),('mass_flow_kg_s',lim['mass_balance_relative']*mf)):
            if not math.isfinite(after[field]) or abs(after[field]-state[field])>tol: raise RuntimeError("Saved/reloaded stream differs")
        max_mu=max(max_mu,abs(state['viscosity_Pa_s']/mu-1))
        max_rho=max(max_rho,abs(state['density_kg_m3']/rho-1))
    mass_errors={};energy_errors={}
    for unit,(ins,outs) in BALANCES.items():
        mass_errors[unit]=sum(s[n]['mass_flow_kg_s'] for n in outs)-sum(s[n]['mass_flow_kg_s'] for n in ins)
        energy_errors[unit]=sum(s[n]['mass_flow_kg_s']*s[n]['enthalpy_J_kg'] for n in outs)-sum(s[n]['mass_flow_kg_s']*s[n]['enthalpy_J_kg'] for n in ins)-(power if unit in ('pump','overall') else 0.)
    close_tree(actual['mass_residuals_kg_s'],mass_errors,'mass residuals')
    close_tree(actual['energy_residuals_W'],energy_errors,'energy residuals')
    # Independent incompressible oracle: full feed is pumped, not just filtrate.
    oracle=mf/rho*(inputs['TMP_Pa']+inputs['Channel_drop_Pa'])/inputs['Pump_efficiency_fraction']
    expected_objects={**{n:'DWSIM.Thermodynamics.Streams.MaterialStream' for n in STREAMS},
        'Pump-power':'DWSIM.UnitOperations.Streams.EnergyStream',
        'P-1':'DWSIM.UnitOperations.UnitOperations.Pump','S-1':'DWSIM.UnitOperations.UnitOperations.Splitter',
        **{n:'DWSIM.UnitOperations.UnitOperations.Valve' for n in ('V-channel','V-membrane','V-retentate')}}
    close_tree(actual['objects'],expected_objects,'stock DWSIM units')
    close_tree(closure['shaft_power_W'],oracle,'Julia process power')
    specific=power/(s['Permeate']['mass_flow_kg_s']/rho)/3.6e6
    checks={
        'all_unit_and_overall_mass_balances':max(map(abs,mass_errors.values()))<=lim['mass_balance_relative']*mf,
        'all_unit_and_overall_energy_balances':max(map(abs,energy_errors.values()))<=lim['energy_balance_relative_to_shaft']*power,
        'pump_independent_oracle':math.isclose(power,oracle,rel_tol=lim['pump_power_relative']),
        'energy_stream':math.isclose(actual['energy_stream_power_W'],power,rel_tol=lim['pump_power_relative']),
        'reload_energy_stream':math.isclose(actual['reloaded_energy_stream_power_W'],power,rel_tol=lim['pump_power_relative']),
        'small_constant_viscosity_approximation':max_mu<=lim['relative_viscosity_change'],
        'small_constant_density_approximation':max_rho<=lim['relative_density_change'],
        'saved_reloaded_solved':actual['save_reload_solve_checked'] is True,
    }
    if not all(checks.values()): raise RuntimeError(f"Filtration acceptance failed: {checks}")
    return {'checks':checks,'max_mass_relative_residual':max(map(abs,mass_errors.values()))/mf,
        'max_energy_relative_residual':max(map(abs,energy_errors.values()))/power,
        'max_relative_viscosity_change':max_mu,'max_relative_density_change':max_rho,
        'actual_specific_shaft_energy_kWh_m3_reference_volume':specific,
        'actual_permeate_reference_flow_m3_s':s['Permeate']['mass_flow_kg_s']/rho,
        'actual_recovery_mass_fraction':s['Permeate']['mass_flow_kg_s']/s['Feed']['mass_flow_kg_s']}
