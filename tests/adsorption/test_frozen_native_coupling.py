"""Frozen-source contract, independently exercised with an explicit unit fixture."""
import json
import pytest
from pore2process_demo.paths import ROOT
from examples.adsorption.column_transport.bridge import REQUIRED_SCREENS
from examples.adsorption.column_transport.coupling import inputs,predict


def fixture(diameter):
    # Hypothetical passing closure tests plumbing only; no actual native pass.
    eps=.6;A=.001;L=.1;Q=1e-7
    p=dict(particle_diameter_m=diameter,length_m=L,area_m2=A,porosity=eps,flow_m3_s=Q,
        dry_mass_kg=(1-eps)*(740/(1-.427))*A*L,
        viscosity_Pa_s=.001,density_kg_m3=1000.,molecular_diffusivity_m2_s=1e-9,
        t_end_s=1.,dt_s=.1,axial_cells=8,radial_cells=4)
    c=dict(boundary='convective-column',native_run_id='unit-fixture-not-native-evidence',
        calibration_boundary='native-inert-tracer-only',axial_dispersion_m2_s=1e-6,
        checks={key:True for key in REQUIRED_SCREENS},operating_point=dict(
            superficial_velocity_m_s=Q/A,particle_diameter_m=diameter,porosity=eps,
            molecular_diffusivity_m2_s=1e-9))
    return c,p


def test_both_particle_sizes_share_frozen_equilibrium_diffusion_and_inventory_rule():
    records=[]
    for diameter in (.0008,.00045):
        c,p=fixture(diameter);r=predict(c,**p);records.append(r)
        assert r['result']['normalized_mass_residual']<1e-8
        assert r['result']['physical_adsorbent_mass_kg']==pytest.approx(p['dry_mass_kg'])
        assert r['result']['active_adsorbent_mass_kg']==pytest.approx(p['dry_mass_kg'])
        assert r['result']['fluid_volume_m3'][0]==pytest.approx(p['porosity']*p['area_m2']*p['length_m'])
        assert r['provenance']['new_column_fitted_parameters']==0
        assert r['provenance']['batch_run_id']=='20261110T083144Z-e33389b9'
    for key,value in dict(capacity_kg_kg=.1115,affinity_m3_kg=139.,
            diffusivity_m2_s=1.846035489738998e-12,active_inventory_fraction=1.).items():
        assert all(r['inputs'][key]==value for r in records)


def test_mismatched_native_point_and_reinterpreted_mass_are_rejected():
    c,p=fixture(.0008);c['operating_point']['particle_diameter_m']=.00045
    with pytest.raises(ValueError,match='operating point'):inputs(c,**p)
    c,p=fixture(.0008);p['dry_mass_kg']*=.573
    with pytest.raises(ValueError,match='dry mass'):inputs(c,**p)


def test_real_failed_packet_never_calls_frozen_particle_prediction():
    c=json.loads((ROOT/'examples/adsorption/column_transport/reference/transport_reference.json').read_text())['closure']
    _,p=fixture(.0008)
    with pytest.raises(ValueError,match='unqualified'):predict(c,**p)
