"""SI contracts and scale translation through the shared Julia package."""
import copy
import math
import numpy as np

from crossflow_recycle.model import adapter_inputs
from graded_support.model import load_case

def finite_number(value):
    try:
        return not isinstance(value,bool) and math.isfinite(value)
    except (TypeError,ValueError):
        return False

def record(jl,value):
    result={}
    for field in jl.fieldnames(jl.typeof(value)):
        key=str(field); v=getattr(value,key)
        result[key]=np.asarray(v,dtype=float).tolist() if key in (
            'x_m','pressure_Pa','local_flux_m_s','channel_flow_m3_s','axial_COD_kg_m3') else float(v)
    return result

def settings(parent):
    p={k:copy.deepcopy(parent['hydraulics'][k]) for k in (
        'temperature_K','external_absolute_pressure_Pa','feed_pump_efficiency_fraction',
        'circulation_pump_efficiency_fraction','recycle','acceptance')}
    p['protocol_id']='explicit-star-water-carrier-equipment-v1'
    p['recycle']['nominal_seed_flow_multiplier']=.5
    p['recycle']['maximum_iterations']=5000
    p['recycle']['mass_tolerance_kg_s']=1e-14
    p['force_liquid_water']=True
    return p

def requests(p):
    nominal=p['process']['nominal_flux_LMH']
    for flux in p['process']['fluxes_LMH']:
        for support in ('straight','widened_5um','widened_6um'):
            yield dict(support=support,flux_LMH=flux,rejection=.65,viscosity_multiplier=1.,fixed_multiplier=1.,inlet_Re=p['process']['inlet_Re'])
    for support in ('straight','widened_5um','widened_6um'):
        for rejection in (.5,.8):
            yield dict(support=support,flux_LMH=nominal,rejection=rejection,viscosity_multiplier=1.,fixed_multiplier=1.,inlet_Re=p['process']['inlet_Re'])
        yield dict(support=support,flux_LMH=nominal,rejection=.65,viscosity_multiplier=1.2,fixed_multiplier=1.,inlet_Re=p['process']['inlet_Re'])
        for fixed in (.5,2.):
            yield dict(support=support,flux_LMH=p['process']['fluxes_LMH'][0],rejection=.65,viscosity_multiplier=1.,fixed_multiplier=fixed,inlet_Re=p['process']['inlet_Re'])
        for re in (10.,30.):
            yield dict(support=support,flux_LMH=nominal,rejection=.65,viscosity_multiplier=1.,fixed_multiplier=1.,inlet_Re=re)

def translate(jl,native,source,p,request):
    q=request; sp=p['process']; f=p['feed']
    for key in ('flux_LMH','viscosity_multiplier','fixed_multiplier'):
        if not finite_number(q[key]) or q[key]<=0: raise ValueError(key+' must be finite and positive.')
    if not 1.<=q['viscosity_multiplier']<=1.2 or not .5<=q['fixed_multiplier']<=2.:
        raise ValueError('Use viscosity multiplier 1..1.2 and fixed-impedance multiplier 0.5..2 for this verified showcase.')
    if not finite_number(q['rejection']) or not 0<=q['rejection']<=1:
        raise ValueError('Rejection must be in [0,1].')
    if not min(sp['fluxes_LMH'])<=q['flux_LMH']<=max(sp['fluxes_LMH']):
        raise ValueError('This verified showcase supports mean flux 8..12 LMH; sustainable flux is unavailable for the new geometry and speed.')
    re=q.get('inlet_Re',sp['inlet_Re'])
    if not finite_number(re) or not sp['verified_Re_range'][0]<=re<=sp['verified_Re_range'][1]:
        raise ValueError('Inlet Reynolds number must stay within the verified 10..30 range.')
    support=next((v for v in source['supports'] if v['name']==q['support']),None)
    if support is None: raise ValueError('Select a verified support geometry.')
    reference=load_case(jl,support['body_resistance_Pa_s_m3'],source['supports'][0]['body_resistance_Pa_s_m3'],
        source['protocol'],source['literature_reference']['fit']['fixed_impedance_Pa_s_m'])
    g=native['geometry']; anchor=source['fluid_anchor']
    load=jl.star_membrane_load(native['cases'][native.get('accepted_case','grid48')]['geometry_coefficient_m4'],g['open_area_m2'],
        g['wetted_perimeter_m'],p['geometry']['module_length_m'],
        anchor['inlet_viscosity_Pa_s']*q['viscosity_multiplier'],anchor['inlet_density_kg_m3'],
        re,q['flux_LMH']/3.6e6,sp['area_m2'],
        reference['support_impedance_Pa_s_m']*q['viscosity_multiplier'],
        reference['non_support_impedance_Pa_s_m']*q['fixed_multiplier']*q['viscosity_multiplier'])
    a=record(jl,load)
    if a['slender_channel_inertia_screen']>.05 or a['wall_suction_velocity_ratio']>.002:
        raise ValueError('Slowly varying channel approximation leaves its declared inertia/suction screens.')
    process=jl.crossflow_recycle_load(a['outlet_TMP_Pa'],a['permeate_flow_m3_s'],sp['recovery'],
        a['module_to_permeate_ratio'],a['channel_resistance_Pa_s_m3'],sp['efficiency'],sp['efficiency'])
    solute=jl.star_feed_balance(f['COD_mg_L']*.001,q['rejection'],load,sp['recovery'])
    return {'axial':a,'process':record(jl,process),'feed':record(jl,solute),
            'reference_membrane':reference,'critical_flux_limit_LMH':None,
            'maximum_local_flux_LMH':max(a['local_flux_m_s'])*3.6e6,
            'critical_flux_status':'unavailable for the synthetic geometry at the verified low Re',
            'source_critical_flux_reference_LMH':f['critical_flux_point']['flux_LMH'],
            'sustainable_flux_LMH':None,'fouling_time_s':None,'measured_feed_viscosity_Pa_s':None,
            'retention_validated_at_operating_point':False,
            'local_flux_admissible':None,'feed_scope':'Clean initial hydraulic case; sustainable operation and physical fouling time are not predicted.'}

def inputs(h,source,t,**kwargs):
    a=adapter_inputs(h,source,t['process'],**kwargs)
    a['Force_liquid']=h['force_liquid_water']
    a['Scope']='Distributed low-Re synthetic star hydraulic ports from OpenLB/Julia; actual DWSIM water-carrier recycle. Julia external illustrative COD retention; sustainable flux unavailable. Stock lumped channel/membrane allocation preserves port flows, pressure and total hydraulic work; internal distributed losses are reported separately.'
    return a

def key(request):
    return f"{request['support']}-J{request['flux_LMH']:g}-R{request['rejection']:g}-mu{request['viscosity_multiplier']:g}-Z{request['fixed_multiplier']:g}-Re{request.get('inlet_Re',20):g}"
