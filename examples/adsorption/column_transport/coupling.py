"""Frozen particle/equilibrium inputs for a future qualified native closure.

Current native evidence is rejected. This module does not provide a replacement
dispersion, fit a column curve, or reinterpret printed density.
"""
import json
import math

from pore2process_demo.paths import ROOT
from pore2process_demo.io import sha
from examples.adsorption.film_particle.bridge import call,bed
from .bridge import qualified_dispersion

PARTICLE=ROOT/'examples/adsorption/film_particle/reference/film_particle_reference.json'
SOURCE=ROOT/'literature/sources/alhamed2009/film_particle_extraction.json'


def inputs(closure, *, particle_diameter_m,length_m,area_m2,porosity,flow_m3_s,
           dry_mass_kg,viscosity_Pa_s,density_kg_m3,molecular_diffusivity_m2_s,
           inlet_concentration_kg_m3=.15,axial_cells=400,radial_cells=32,
           t_end_s=18000.,dt_s=.5):
    D=qualified_dispersion(closure)
    for name,value in dict(particle_diameter_m=particle_diameter_m,length_m=length_m,
            area_m2=area_m2,porosity=porosity,flow_m3_s=flow_m3_s,dry_mass_kg=dry_mass_kg,
            viscosity_Pa_s=viscosity_Pa_s,density_kg_m3=density_kg_m3,
            molecular_diffusivity_m2_s=molecular_diffusivity_m2_s).items():
        if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value) or value<=0:
            raise ValueError(name+' must be positive finite SI input.')
    if porosity>=1:raise ValueError('Porosity must be below 1.')
    source=json.loads(SOURCE.read_text());reference=json.loads(PARTICLE.read_text())
    status=reference['physical_status']
    if not (status['source_equilibrium_audit_passed'] and status['batch_particle_size_comparison_passed']):
        raise ValueError('Frozen batch/equilibrium component is not qualified.')
    e=source['source_parameters'];fit=reference['batch_fit']
    u=flow_m3_s/area_m2
    point=closure.get('operating_point',{})
    for key,value in dict(superficial_velocity_m_s=u,porosity=porosity,
            particle_diameter_m=particle_diameter_m,molecular_diffusivity_m2_s=molecular_diffusivity_m2_s).items():
        saved=point.get(key)
        if isinstance(saved,bool) or not isinstance(saved,(float,int)) or not math.isfinite(saved) or not math.isclose(saved,value,rel_tol=1e-4,abs_tol=0.):
            raise ValueError('Native closure operating point does not match '+key)
    expected=(1-porosity)*e['particle_density_kg_m3']*area_m2*length_m
    if not math.isclose(dry_mass_kg,expected,rel_tol=1e-10,abs_tol=0.):
        raise ValueError('Declared dry mass must match the full geometry and frozen particle-envelope density.')
    film=call('film_coefficient',diameter_m=particle_diameter_m,superficial_velocity_m_s=u,
        viscosity_Pa_s=viscosity_Pa_s,density_kg_m3=density_kg_m3,
        molecular_diffusivity_m2_s=molecular_diffusivity_m2_s,porosity=porosity,
        interstitial=False,multiplier=1.)
    p=dict(particle_diameter_m=particle_diameter_m,length_m=length_m,area_m2=area_m2,
        porosity=porosity,superficial_velocity_m_s=u,inlet_concentration_kg_m3=inlet_concentration_kg_m3,
        particle_density_kg_m3=e['particle_density_kg_m3'],
        diffusivity_m2_s=fit['diffusivity_m2_s'],film_m_s=film['film_m_s'],
        capacity_kg_kg=e['capacity_kg_kg'],affinity_m3_kg=e['affinity_m3_kg'],
        axial_dispersion_m2_s=D,active_inventory_fraction=1.,axial_cells=axial_cells,
        radial_cells=radial_cells,t_end_s=t_end_s,dt_s=dt_s,
        initial_concentration_kg_m3=0.,initial_loading_kg_kg=0.,linear_partition_m3_kg=None)
    trace=dict(native_run_id=closure['native_run_id'],native_closure=closure,
        batch_run_id=reference['run_id'],particle_reference_sha256=sha(PARTICLE),
        source_equilibrium_sha256=sha(SOURCE),source_doi=source['doi'],
        source_printed_bulk_density_kg_m3=e['bulk_density_kg_m3'],
        source_printed_voidage=e['porosity'],
        scenario_dry_mass_kg=dry_mass_kg,scenario_bulk_density_kg_m3=(1-porosity)*e['particle_density_kg_m3'],
        scenario_pore_volume_m3=porosity*area_m2*length_m,
        particle_diffusivity_calibration='Figure 2 batch-only fit, frozen across diameters',
        column_film_role='Previously declared Wakao-form correlation scenario; superficial velocity, multiplier=1',
        batch_film_ref_m_s_not_applied_to_column=fit['batch_film_ref_m_s'],
        column_film=film,new_column_fitted_parameters=0,full_active_inventory=True,
        material_transfer='Synthetic packing with transferred batch-qualified particle properties; not experimental column validation')
    return p,trace


def predict(closure, **parameters):
    p,trace=inputs(closure,**parameters)
    result=bed(**p)
    return dict(inputs=p,provenance=trace,result=result)
