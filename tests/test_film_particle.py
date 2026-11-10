"""Real Julia forward checks and traceable calibration/holdout contracts."""
from concurrent.futures import ThreadPoolExecutor
import json

import numpy as np
import pytest

from examples.adsorption.film_particle.bridge import batch, bed, call
from examples.adsorption.film_particle.check import verify_hashes, inventories
from examples.adsorption.film_particle.data import ROOT, evidence, observations


def reference():
    return json.loads((ROOT/'docs/film_particle_reference.json').read_text())


def test_source_split_and_unreported_uncertainty_stay_explicit():
    e=evidence();obs=observations()
    assert sum(len(d['time_s']) for d in obs.values())==24
    assert {k for k,d in obs.items() if d['role']=='calibration'}=={'d800','d1500'}
    assert {k for k,d in obs.items() if d['role']=='holdout'}=={'d450','d225'}
    assert e['figure_2']['pointwise_measurement_uncertainty'] is None
    assert e['inventory']['mass_role'].startswith('density times')
    assert e['source_parameters']['capacity_kg_kg']==.1115
    assert e['source_parameters']['affinity_m3_kg']==139.
    for data in obs.values():
        assert np.all(np.diff(data['time_s'])>0)
        assert np.all(data['loading_kg_kg']>.03)


@pytest.mark.parametrize('index',range(8))
def test_actual_julia_replays_batch_and_column_reference(index):
    record=reference()['cases'][index]
    actual=(batch if record['id'].startswith('batch_') else bed)(**record['inputs'])
    for key,value in record['result_scalars'].items():
        assert actual[key]==pytest.approx(value,rel=1e-9,abs=1e-13)
    assert actual['mean_loading_kg_kg'][-1]==pytest.approx(record['last_loading_kg_kg'],rel=1e-9)
    assert actual['outlet_concentration_kg_m3'][-1]==pytest.approx(record['last_outlet_kg_m3'],rel=1e-9)
    inventories({**record,'result':actual})
    q=np.array(actual['final_shell_loading_kg_kg'])
    assert q.dtype==np.float64
    assert q.shape==(record['inputs'].get('axial_cells',1),record['inputs']['radial_cells'])


def test_finite_bath_endpoint_has_independent_mass_balance():
    r=call('bath_equilibrium',initial_concentration_kg_m3=.4,fluid_volume_m3=.00025,
           adsorbent_mass_kg=.001,capacity_kg_kg=.1115,affinity_m3_kg=139.)
    assert .00025*r['concentration_kg_m3']+.001*r['loading_kg_kg']==pytest.approx(.0001,rel=1e-12)
    assert r['loading_kg_kg']<.1


def test_withdrawal_carries_solute_without_a_mass_correction():
    p=dict(t_end_s=120.,dt_s=1.,sampling_times_s=[60.,120.],sample_volume_m3=1e-6)
    r=batch(**p)
    removed=r['cumulative_sampled_solute_kg'][-1]
    assert removed>0
    assert r['fluid_volume_m3'][-1]==pytest.approx(.000248)
    assert r['inventory_kg'][-1]+removed==pytest.approx(.00025*.2,rel=1e-9)
    assert r['normalized_mass_residual']<1e-8


def test_one_column_factor_does_not_refit_equilibrium_inventory_or_diffusivity():
    ref=reference();cases={c['id']:c for c in ref['cases']}
    a=cases['anchored_calibration']['inputs'];b=cases['anchored_holdout']['inputs']
    for key in ['diffusivity_m2_s','capacity_kg_kg','affinity_m3_kg',
                'particle_density_kg_m3','porosity','active_inventory_fraction']:
        assert a[key]==b[key]
    assert a['active_inventory_fraction']==1.
    assert set(ref['batch_fit']['fitted_diameter_ids'])=={'d800','d1500'}
    assert ref['endpoint_audit']['diagnostic_factor_applied_to_forward_model'] is False
    assert ref['physical_status']['batch_particle_size_comparison_passed'] is True
    assert ref['physical_status']['full_transport_comparison_passed'] is False
    assert ref['numerical_refinement_passed'] is True


def test_zero_exchange_is_a_conservative_transport_limit():
    p=dict(t_end_s=100.,dt_s=.5,length_m=.01,area_m2=.001,porosity=.5,
        superficial_velocity_m_s=.001,axial_cells=20,radial_cells=8)
    nofilm=bed(**p,film_m_s=0.)
    nodiffusion=bed(**p,diffusivity_m2_s=0.)
    np.testing.assert_allclose(nofilm['outlet_concentration_kg_m3'],nodiffusion['outlet_concentration_kg_m3'])
    assert all(x==0 for x in nofilm['adsorbed_inventory_kg'])


def test_background_thread_is_rejected_even_after_module_initialization():
    call('bath_equilibrium')
    with ThreadPoolExecutor(max_workers=1) as pool:
        with pytest.raises(RuntimeError,match='main server thread'):
            pool.submit(batch,t_end_s=1.,dt_s=1.).result()


@pytest.mark.parametrize('kwargs',[{'diffusivity_m2_s':-1.},{'film_m_s':float('nan')},{'radial_cells':False}])
def test_invalid_particle_inputs_fail_before_solving(kwargs):
    with pytest.raises(Exception,match='diffusivity|film|radial_cells'):
        batch(**kwargs)


def test_changed_saved_artifact_is_rejected(tmp_path):
    from examples.adsorption.alhamed.source import sha
    p=tmp_path/'field.csv';p.write_text('q_kg_kg\n0\n');digest=sha(p)
    p.write_text('q_kg_kg\n1\n')
    with pytest.raises(ValueError,match='Changed or invalid artifact'):
        verify_hashes(tmp_path,{'field.csv':digest})
