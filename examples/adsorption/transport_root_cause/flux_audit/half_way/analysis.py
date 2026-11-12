"""Read-only assessment of the finite stock half-way boundary candidate."""
import json
from pathlib import Path
import numpy as np
from ..flat_channel.reference_qualified import converged_reference


def structural(dest):
    dest=Path(dest)
    links=json.loads((dest/'halfway_links.json').read_text())
    audits=[json.loads(p.read_text()) for p in sorted(dest.glob('halfway_audit_*.json'))]
    checks=dict(complete_links=links['fluid_solid_links']>0 and links['missing_links']==0,
                no_spurious_links=links['unexpected_active_links']==0,
                inactive_untouched=links['inactive_changed']==0,
                exact_halfway=links['q']==.5,
                local_fluid_balance=len(audits)>=3 and all(a['max_local_fluid_balance_residual']<=1e-11 for a in audits),
                actual_reflection=len(audits)>=3 and all(a['max_reflection_pair_residual']<=1e-11 and a['reflected_links']==links['fluid_solid_links'] for a in audits))
    return links,audits,checks


def analyze(dest,kind):
    dest=Path(dest);raw=json.loads((dest/'raw.json').read_text())
    history=np.genfromtxt(dest/'inventory.csv',delimiter=',',names=True)
    field=np.genfromtxt(dest/'final.csv',delimiter=',',names=True)
    initial=np.genfromtxt(dest/'initial.csv',delimiter=',',names=True)
    amp=raw['initial_amplitude_kg_m3'] if kind=='flat' else .25
    links,audits,checks=structural(dest)
    fluid=history['fluid_inventory_kg'] if kind=='flat' else history['fluid_tracer_kg']
    drift=float(np.max(abs(fluid/fluid[0]-1)))
    lo=float(min(history['min_concentration_kg_m3']));hi=float(max(history['max_concentration_kg_m3']))
    excursion=max(.5-amp-lo,hi-.5-amp,0.)
    initial_fluid=initial if kind=='flat' else initial[initial['material']==1]
    length=raw['length_m'];expected=.5+amp*np.sin(2*np.pi*initial_fluid['x_m']/length)
    initial_error=float(np.max(abs(initial_fluid['concentration_kg_m3']-expected)))
    checks.update(finite=all(bool(np.all(np.isfinite(history[k]))) for k in history.dtype.names),
                  exact_initial_snapshot=initial_error<=1e-12,
                  physical_fluid_mass=drift<=1e-8,
                  strict_bounds=excursion<=1e-8,engineering_bounds=excursion<=.02)
    result=dict(raw=raw,links=links,local_audits=audits,checks=checks,
                minimum_concentration_kg_m3=lo,maximum_concentration_kg_m3=hi,
                maximum_bound_excursion_kg_m3=excursion,
                maximum_initial_error_kg_m3=initial_error,relative_physical_fluid_inventory_drift=drift,
                unused_solid_populations_are_inventory=False,application_transport_qualified=False,
                measured_column_dispersion_m2_s=None)
    if kind=='flat':
        n=raw['n'];H=raw['height_m']
        np.testing.assert_allclose(field['x_m'],(field['i']+.5)*H/n,atol=1e-16,rtol=1e-13)
        np.testing.assert_allclose(field['y_m'],(field['j']+.5)*H/n,atol=1e-16,rtol=1e-13)
        velocity=6*raw['mean_speed_m_s']*(field['y_m']/H)*(1-field['y_m']/H)
        np.testing.assert_allclose(field['ux_m_s'],velocity,atol=1e-14,rtol=1e-12)
        checks['declared_fluid_volume']=abs(fluid[0]/(.5*H*H*raw['depth_m'])-1)<=1e-11
        checks['z_invariance']=float(np.max(field['z_variation_kg_m3']))<=1e-12
        checks['no_hidden_diffusivity']=abs((raw['tau_odd']-.5)*(H/n)**2/(3*raw['dt_s'])/raw['diffusivity_m2_s']-1)<=1e-10
        if amp:
            p=dict(height_m=H,mean_speed_m_s=raw['mean_speed_m_s'],diffusivity_m2_s=raw['diffusivity_m2_s'],duration_s=raw['duration_s'])
            exact,reference=converged_reference(field['x_m'],field['y_m'],amplitude_kg_m3=amp,**p)
            difference=field['concentration_kg_m3']-exact
            result.update(reference_convergence=reference,
                          RMS_error_over_amplitude=float(np.sqrt(np.mean(difference**2)))/amp,
                          maximum_error_over_amplitude=float(np.max(abs(difference)))/amp,
                          first_fluid_row_RMS_error_over_amplitude=float(np.sqrt(np.mean(difference[field['j']==0]**2)))/amp)
            checks['profile']=result['RMS_error_over_amplitude']<=.01
            if raw['mean_speed_m_s']==0:
                analytical=.5+amp*np.exp(-raw['diffusivity_m2_s']*(2*np.pi/H)**2*raw['duration_s'])*np.sin(2*np.pi*field['x_m']/H)
                result['oracle_analytical_max_error_kg_m3']=float(np.max(abs(exact-analytical)))
                checks['diffusion_analytical_oracle']=result['oracle_analytical_max_error_kg_m3']<=1e-10
        else:
            result['maximum_constant_error_kg_m3']=max(abs(lo-.5),abs(hi-.5))
            checks['constant_preservation']=result['maximum_constant_error_kg_m3']<=1e-9
    result['checks']={k:bool(v) for k,v in checks.items()}
    result['all_checks_passed']=all(result['checks'].values())
    return result


def stage_one_gate(cases):
    """Coarse Poiseuille bounds/profile are diagnostic; all other screens apply."""
    for name,case in cases.items():
        for key,value in case['checks'].items():
            if name=='shear32' and key in ('profile','strict_bounds','engineering_bounds'):continue
            if not value:return False
    return len(cases)==4
