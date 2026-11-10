"""Native-measurement to Julia to DWSIM contracts; no Python flow solver."""
import math
from itertools import product
from examples.filtration.graded_support.model import load_case
from examples.filtration.star_channel.model import record, finite_number, inputs as port_inputs


def requests(p):
    for re,flux,fixed,support in product(p['process']['Re_targets'],p['process']['fluxes_LMH'],
            p['process']['fixed_multipliers'],('straight','widened_5um','widened_6um')):
        yield dict(inlet_Re=float(re),flux_LMH=float(flux),fixed_multiplier=float(fixed),support=support)
    for re,support in product(p['process']['Re_targets'],('straight','widened_5um','widened_6um')):
        yield dict(inlet_Re=float(re),flux_LMH=12.,fixed_multiplier=1.,support=support,viscosity_multiplier=1.2)


def key(q):
    return f"{q['support']}-J{q['flux_LMH']:g}-Z{q['fixed_multiplier']:g}-Re{q['inlet_Re']:g}-mu{q.get('viscosity_multiplier',1.):g}"


def translate(jl,native,source,p,q):
    for k in ('inlet_Re','flux_LMH','fixed_multiplier'):
        if not finite_number(q.get(k)) or q[k]<=0:raise ValueError(k+' must be finite and positive.')
    if not min(p['process']['fluxes_LMH'])<=q['flux_LMH']<=max(p['process']['fluxes_LMH']):
        raise ValueError('Mean flux is outside the declared operating-map range.')
    if not .5<=q['fixed_multiplier']<=2.:raise ValueError('Fixed-impedance multiplier must be 0.5..2.')
    lo,hi=native['verified_Re_range']
    if not lo<=q['inlet_Re']<=hi:raise ValueError('Reynolds number outside the native verified envelope.')
    viscosity=q.get('viscosity_multiplier',1.)
    if not finite_number(viscosity) or not 1.<=viscosity<=1.2:raise ValueError('Hydraulic-viscosity multiplier must be 1..1.2.')
    support=next((s for s in source['supports'] if s['name']==q.get('support')),None)
    if support is None:raise ValueError('Unknown verified TRT support.')
    ref=load_case(jl,support['body_resistance_Pa_s_m3'],source['supports'][0]['body_resistance_Pa_s_m3'],
        source['protocol'],source['literature_reference']['fit']['fixed_impedance_Pa_s_m'])
    g=native['geometry'];anchor=source['fluid_anchor'];sp=p['process']
    load=jl.star_membrane_load(native['geometry_coefficient_m4'],g['open_area_m2'],g['wetted_perimeter_m'],sp.get('module_length_m',.3),
        anchor['inlet_viscosity_Pa_s']*viscosity,anchor['inlet_density_kg_m3'],q['inlet_Re'],q['flux_LMH']/3.6e6,sp['area_m2'],
        ref['support_impedance_Pa_s_m']*viscosity,ref['non_support_impedance_Pa_s_m']*q['fixed_multiplier']*viscosity)
    a=record(jl,load)
    measured=native['cases'][native['accepted_case']]
    beta_bound=measured['raw']['max_speed_m_s']/abs(measured['mean_velocity_m_s'])
    # Momentum-flux correction beta <= umax/umean for nonnegative axial flow.
    # With q'=-P*J, the neglected axial acceleration / friction is bounded by
    # 2*beta*rho*C*P*J/(mu*A^2) under the local-profile approximation.
    screen=float(jl.axial_momentum_screen(a['geometry_coefficient_m4'],g['open_area_m2'],g['wetted_perimeter_m'],a['viscosity_Pa_s'],a['density_kg_m3'],max(a['local_flux_m_s']),beta_bound))
    a['momentum_shape_factor_upper_bound']=beta_bound
    a['neglected_axial_momentum_to_friction_upper_screen']=screen
    if screen>.02:raise ValueError('Maximum-local-flux momentum/friction screen exceeds 2%.')
    if a['slender_channel_inertia_screen']>.05 or a['wall_suction_velocity_ratio']>.002:
        raise ValueError('Slowly varying laminar-channel inertia/suction screen failed.')
    process=jl.crossflow_recycle_load(a['outlet_TMP_Pa'],a['permeate_flow_m3_s'],sp['recovery'],
        a['module_to_permeate_ratio'],a['channel_resistance_Pa_s_m3'],sp['efficiency'],sp['efficiency'])
    cm=sp.get('module_inlet_COD_kg_m3',.0687)
    feed=jl.star_feed_balance(cm,sp['rejection'],load,sp['recovery'])
    return {'axial':a,'process':record(jl,process),'feed':record(jl,feed),'reference_membrane':ref,
        'quality_sensitivity':[{'assumed_local_rejection':r,'balance':record(jl,jl.star_feed_balance(cm,r,load,sp['recovery'])),
            'scope':'External passive COD balance on the same hydraulic result; chosen retention scenario, not another DWSIM species simulation.'} for r in (.5,.8)],
        'sustainable_flux_LMH':None,'fouling_time_s':None,'measured_feed_viscosity_Pa_s':None,
        'retention_validated_at_operating_point':False}


def inputs(h,source,t,**kw):
    a=port_inputs(h,source,t,**kw)
    a['Scope']='Higher-Re synthetic fully developed laminar-branch OpenLB/Julia hydraulic ports; actual DWSIM water-carrier recycle. Critical-flux transfer screened separately; no physical fouling time or sustainable-operation prediction.'
    return a


def capacity_decisions(jl,t,evidence,p):
    points=evidence['points'];vec=jl.seval('Vector{Float64}')
    arrays=[vec([v[k] for v in points]) for k in ('source_Re',)]
    arrays += [vec([v['critical_flux_LMH']/3.6e6 for v in points])]
    arrays += [vec([v['reading_interval_LMH'][i]/3.6e6 for v in points]) for i in (0,1)]
    a=t['axial'];peak=max(a['local_flux_m_s']);results=[]
    for mapping,fraction,factor in product(p['capacity']['mappings'],p['capacity']['fractions'],p['capacity']['transfer_factors']):
        query=a['inlet_Re'] if mapping=='equal_Re' else evidence['reference_velocity_pair']['Re']*a['inlet_velocity_m_s']/evidence['reference_velocity_pair']['velocity_m_s']
        row={'mapping':mapping,'capacity_fraction':fraction,'transfer_factor':factor,'source_query_Re':query,
             'maximum_local_flux_LMH':peak*3.6e6,'sustainable_operation_validated':False}
        if not points[0]['source_Re']<=query<=points[-1]['source_Re']:
            results.append({**row,'status':'unavailable','reason':'Outside measured source range; no extrapolation.',
                            'nominal_admissible':None,'linear_reading_admissible':None,'monotonic_classification':'unavailable'})
            continue
        v=jl.critical_flux_capacity(*arrays,query,peak,fraction,factor)
        out={}
        for name in jl.fieldnames(jl.typeof(v)):
            k=str(name);val=getattr(v,k)
            out[k]=str(val) if k=='monotonic_classification' else bool(val) if k.endswith('_admissible') else int(val) if k.endswith('_index') else float(val)
        results.append({**row,**out,'status':'evaluated'})
    return results


def comparison(cases):
    rows=[]
    for c in cases:
        q=c['request'];a=c['translation']['axial'];v=c['validation'];f=c['translation']['feed']
        base=next((b for b in cases if b['request']=={**q,'support':'straight'}),None)
        energy=v['specific_shaft_energy_kWh_m3_reference_volume']
        rows.append({'case_id':c['id'],'viscosity_multiplier':q.get('viscosity_multiplier',1.),**q,'specific_shaft_energy_kWh_m3':energy,
            'support_energy_saving_fraction':None if base is None else 1-energy/base['validation']['specific_shaft_energy_kWh_m3_reference_volume'],
            'fixed_resistance_fraction':a['selective_and_other_fixed_fraction'],
            'mean_TMP_Pa':a['mean_TMP_Pa'],'maximum_local_flux_LMH':max(a['local_flux_m_s'])*3.6e6,
            'inlet_TMP_Pa':a['inlet_TMP_Pa'],'outlet_TMP_Pa':a['outlet_TMP_Pa'],
            'channel_drop_Pa':a['channel_drop_Pa'],'inlet_velocity_m_s':a['inlet_velocity_m_s'],
            'feed_power_W':c['dwsim']['powers']['feed_W'],'circulation_power_W':c['dwsim']['powers']['circulation_W'],
            'circulation_power_fraction':c['dwsim']['powers']['circulation_W']/v['total_shaft_power_W'],
            'permeate_COD_mg_L':f['permeate_COD_kg_m3']*1000,
            'module_permeate_ratio':a['module_to_permeate_ratio'],
            'viscosity_Pa_s':a['viscosity_Pa_s'],'geometry_coefficient_m4':a['geometry_coefficient_m4']})
    return rows


def source_shape_diagnostics(jl,evidence):
    """Leave-one-interior-point-out shape check, not independent validation."""
    points=evidence['points'];vec=jl.seval('Vector{Float64}');rows=[]
    for index in range(1,len(points)-1):
        held=points[index];training=[v for i,v in enumerate(points) if i!=index]
        x=vec([v['source_Re'] for v in training]);y=vec([v['critical_flux_LMH']/3.6e6 for v in training])
        v=jl.critical_flux_capacity(x,y,y,y,held['source_Re'],held['critical_flux_LMH']/3.6e6,1.,1.)
        predicted=float(v.nominal_capacity_m_s)*3.6e6
        rows.append({'held_out_source_Re':held['source_Re'],'reported_or_digitized_flux_LMH':held['critical_flux_LMH'],
            'interpolated_without_point_LMH':predicted,'relative_error':predicted/held['critical_flux_LMH']-1,
            'scope':'Internal interpolation-shape diagnostic from the same four-point series; not independent validation or experimental uncertainty.'})
    return rows
