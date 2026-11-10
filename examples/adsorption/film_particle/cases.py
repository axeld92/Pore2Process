"""SI inputs with source values, geometry and fitted transport kept distinct."""
import math

from .bridge import call


def batch_inputs(diameter, diffusivity, film_ref, *, radial_cells=24, dt_s=30.):
    return dict(initial_concentration_kg_m3=.2,fluid_volume_m3=.00025,
        adsorbent_mass_kg=.001,particle_diameter_m=float(diameter),
        particle_density_kg_m3=740/(1-.427),capacity_kg_kg=.1115,affinity_m3_kg=139.,
        diffusivity_m2_s=float(diffusivity),film_m_s=float(film_ref*.0008/diameter),
        radial_cells=radial_cells,t_end_s=180000.,dt_s=dt_s,initial_loading_kg_kg=0.,
        linear_partition_m3_kg=None,sampling_times_s=[],sample_volume_m3=0.)


def column_inputs(diameter, diffusivity, *, multiplier=1.,interstitial=False,
                  dispersion_factor=0.,active_fraction=1.,axial_cells=80,radial_cells=16,dt_s=1.):
    A=math.pi*.011**2/4;u=43.8e-6/60/A
    film=call('film_coefficient',diameter_m=float(diameter),superficial_velocity_m_s=u,
        viscosity_Pa_s=1/1200,density_kg_m3=1000.,molecular_diffusivity_m2_s=1e-9,
        porosity=.427,interstitial=interstitial,multiplier=float(multiplier))
    p=dict(length_m=.1,area_m2=A,porosity=.427,superficial_velocity_m_s=u,
        inlet_concentration_kg_m3=.15,particle_diameter_m=float(diameter),
        particle_density_kg_m3=740/(1-.427),capacity_kg_kg=.1115,affinity_m3_kg=139.,
        diffusivity_m2_s=float(diffusivity),film_m_s=film['film_m_s'],
        axial_dispersion_m2_s=float(dispersion_factor*(u/.427)*diameter),
        active_inventory_fraction=float(active_fraction),axial_cells=axial_cells,
        radial_cells=radial_cells,t_end_s=18000.,dt_s=dt_s,
        initial_concentration_kg_m3=0.,initial_loading_kg_kg=0.,linear_partition_m3_kg=None)
    film['velocity_convention']='interstitial' if interstitial else 'superficial'
    film['multiplier']=multiplier
    film['dispersion_scenario_factor']=dispersion_factor
    return p,film
