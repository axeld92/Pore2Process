"""Read saved native data; no replacement pore-scale solver."""
import json
from pathlib import Path
import numpy as np
from .bridge import moments,call


def analyze(path):
    path=Path(path);raw=json.loads((path/'raw.json').read_text())
    flux=np.genfromtxt(path/'flux.csv',delimiter=',',names=True)
    inv=np.genfromtxt(path/'inventory.csv',delimiter=',',names=True)
    m=moments(flux['time_s'].tolist(),flux['outlet_flux_kg_s'].tolist())
    total=inv['fluid_tracer_kg']+inv['wall_population_tracer_kg']
    residual=total-inv['cumulative_inlet_kg']+inv['cumulative_outlet_kg']
    dt=raw['dt_s'];inc=np.cumsum(flux['inlet_flux_kg_s'])*dt;out=np.cumsum(flux['outlet_flux_kg_s'])*dt
    indices=inv['step'].astype(int)[1:]-1
    np.testing.assert_allclose(inc[indices],inv['cumulative_inlet_kg'][1:],rtol=1e-10,atol=1e-20)
    np.testing.assert_allclose(out[indices],inv['cumulative_outlet_kg'][1:],rtol=1e-10,atol=1e-20)
    np.testing.assert_allclose(residual,inv['mass_residual_kg'],rtol=1e-8,atol=1e-20)
    mass=float(np.max(abs(residual))/inv['cumulative_inlet_kg'][-1])
    checks=dict(mass=mass<=1e-8,positive=float(np.min(inv['min_concentration_kg_m3']))>=-1e-8,
        bounded=float(np.max(inv['max_concentration_kg_m3']))<=1+1e-6,
        final_signal=m['final_fraction']>=.999,tail_area=m['tail_area_fraction']<.01,
        tail_second_moment=m['tail_second_moment_fraction']<.01)
    transport_length=(raw['nx']-2)*raw['dx_m']
    data=dict(raw=raw,transport_length_m=transport_length,moments=m,checks=checks,normalized_mass_residual=mass,
        fluid_void_time_s=raw['fluid_volume_m3']/raw['flow_m3_s'],
        wall_population_fraction_final=float(inv['wall_population_tracer_kg'][-1]/total[-1]),
        final_inventory_kg=float(total[-1]),
        steady_exit_to_volume_flow_ratio=m['steady_exit_flux_kg_s']/raw['flow_m3_s'],
        scope='Synthetic inert-tracer EXIT-FLUX step with absorbing end; no adsorption or measured experimental dispersion.')
    if raw['uniform']:
        v=raw['flow_m3_s']/raw['area_m2']
        exact=call('absorbing_moments',length_m=transport_length,velocity_m_s=v,dispersion_m2_s=raw['diffusivity_m2_s'])
        data['analytical']=exact
        data['mean_relative_error']=abs(m['mean_s']/exact['mean_s']-1)
        data['variance_relative_error']=abs(m['variance_s2']/exact['variance_s2']-1)
        ref=call('solve_ade',length_m=transport_length,velocity_m_s=v,dispersion_m2_s=raw['diffusivity_m2_s'],
            t_end_s=float(flux['time_s'][-1]),dt_s=.002,cells=800)
        native=flux['outlet_flux_kg_s']/m['steady_exit_flux_kg_s']
        reference=np.interp(flux['time_s'],ref['time_s'],np.array(ref['outlet_flux_per_area_m_s'])/exact['steady_flux_per_area_m_s'])
        data['reference_RMSE']=float(np.sqrt(np.mean((native-reference)**2)))
        checks.update(analytical_mean=data['mean_relative_error']<=.02,analytical_variance=data['variance_relative_error']<=.05,
            reference_profile=data['reference_RMSE']<=.02)
    return data
