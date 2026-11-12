"""Unchanged native assessment with the separately refined continuum oracle."""
import json
from pathlib import Path
import numpy as np
from .reference_qualified import converged_reference

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

