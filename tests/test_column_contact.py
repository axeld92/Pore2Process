"""Real Julia contact model and explicit separation of evidence and fitting."""
from concurrent.futures import ThreadPoolExecutor
import json
import numpy as np
import pytest
from adsorption.alhamed.source import ROOT
from adsorption.column_contact.audit import audit,frozen
from adsorption.column_contact.bridge import contact,module


def reference():
    return json.loads((ROOT/'docs/column_contact_reference.json').read_text())


def test_printed_inventory_is_a_derived_assumption_not_a_weighed_mass():
    a=audit()
    assert a['printed_bulk_density_kg_m3']==740.
    assert a['derived_bed_mass_kg']==pytest.approx(.007032455155060752)
    assert a['per_column_weighed_mass_kg'] is None
    assert a['density_method'] is None and a['per_fraction_voidage'] is None
    assert a['alternative_label_diagnostic']['forward_use'] is False
    assert a['geometry']['void_residence_time_s']==pytest.approx(5.5587899874323226)
    assert a['inlet_section']['liquid_residence_upper_bound_s']<7
    assert a['applied_previous_column_film_multiplier'] is False
    for r in a['paper_dispersion']:
        assert r['used_as_model_input'] is False
        assert 'not independently measured' in r['role']
    for r in a['observed_windows'].values():
        assert r['finite_window_retained_kg']<a['equilibrium']['full_equilibrium_inventory_kg']
        assert r['minimum_additional_monotone_tail_s']>36000


@pytest.mark.parametrize('index',range(8))
def test_real_julia_contact_and_packing_scenarios(index):
    record=reference()['cases'][index];r=contact(**record['inputs'])
    for key,value in record['result_scalars'].items():
        assert r[key]==pytest.approx(value,rel=1e-9,abs=1e-13)
    assert r['outlet_concentration_kg_m3'][-1]==pytest.approx(record['last_outlet_kg_m3'],rel=1e-9)
    assert r['mean_loading_kg_kg'][-1]==pytest.approx(record['last_loading_kg_kg'],rel=1e-9)
    q=np.array(r['final_shell_loading_kg_kg'])
    assert q.dtype==np.float64
    assert q.shape==(record['inputs']['sections'],record['inputs']['radial_cells'])
    assert r['physical_adsorbent_mass_kg']==r['active_adsorbent_mass_kg']
    assert r['normalized_mass_residual']<1e-8
    assert r['normalized_interface_flux_mismatch']<1e-8


def test_same_calibration_boundary_and_frozen_particles():
    old=frozen();ref=reference();fit=ref['selection']
    assert fit['fitted_curve']=='Figure 7, 0.8 mm only' and fit['holdout_used'] is False
    assert fit['particle_film_equilibrium_inventory_fitted'] is False
    selected=[r for r in ref['cases'] if r['kind']=='selected']
    assert len(selected)==2
    for r in ref['cases']:
        p=r['inputs']
        assert p['sections']==fit['selected_sections']
        assert p['diffusivity_m2_s']==old['batch_fit']['diffusivity_m2_s']
        assert p['particle_density_kg_m3']==740/(1-.427)
        assert 'active_inventory_fraction' not in p and 'axial_dispersion_m2_s' not in p
    assert ref['numerical_refinement_passed'] is True
    assert ref['physical_status']['common_column_comparison_passed'] is False


def test_nonadsorbing_single_contact_section_has_cstr_response():
    p=dict(sections=1,film_m_s=0.,t_end_s=20.,dt_s=.001,length_m=.01,
        area_m2=.001,porosity=.5,superficial_velocity_m_s=.001,inlet_concentration_kg_m3=.2)
    r=contact(**p);t=np.array(r['time_s']);fraction=np.array(r['outlet_concentration_kg_m3'])/.2
    np.testing.assert_allclose(fraction,1-np.exp(-t/5),rtol=0,atol=.0001)
    assert r['section_liquid_volume_m3']==pytest.approx(5e-6)
    assert r['liquid_residence_time_s']==pytest.approx(5.)


def test_runtime_guard_applies_after_loading():
    module()
    with ThreadPoolExecutor(max_workers=1) as pool:
        with pytest.raises(RuntimeError,match='main server thread'):
            pool.submit(contact,t_end_s=1.,dt_s=1.).result()


@pytest.mark.parametrize('kwargs',[{'sections':False},{'sections':0},{'porosity':1.},{'film_m_s':float('nan')}])
def test_bad_inputs_do_not_launch_a_valid_case(kwargs):
    with pytest.raises(Exception,match='sections|porosity|film'):
        contact(**kwargs)
