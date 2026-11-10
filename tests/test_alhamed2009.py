"""Source contracts and actual Julia replay for the separate date-stone carbon."""
import json

import numpy as np
import pytest

from adsorption.bridge import column
from adsorption.alhamed.source import ROOT, evidence, points
from adsorption.alhamed.study import crossing, compare


def reference():
    return json.loads((ROOT/'docs/alhamed2009_reference.json').read_text())


def test_source_si_units_bulk_density_and_paired_conditions():
    e=evidence();values={r['name']:r for r in e['measurements_and_fits']}
    assert values['Langmuir_capacity']['SI_value']==pytest.approx(.1115)
    assert values['Langmuir_affinity']['SI_value']==139.
    assert values['liquid_viscosity']['SI_value']==pytest.approx(8.333333333333333e-4)
    f=e['figure_7'];A=np.pi*.011**2/4
    assert f['flow_m3_s']==pytest.approx(7.3e-7)
    assert f['derived_dry_mass_kg']==pytest.approx(740*A*.1)
    assert f['pointwise_uncertainty'] is None
    assert e['unavailable']['pressure_drop_Pa'] is None
    rows=[r for r in e['column_table_5'] if r['length_cm']==10 and r['flow_mL_min']==43.8 and r['inlet_mg_L']==150]
    assert {r['diameter_mm'] for r in rows}=={.45,.8}
    assert all(abs(r['derived_Pe']-r['reported_Pe'])<max(.001,.02*r['reported_Pe']) for r in e['column_table_5'])
    assert min(r['fitted_R2'] for r in e['column_table_5'])==.919


def test_digitized_points_are_ordered_with_explicit_reading_sensitivity():
    p=points();f=evidence()['figure_7']
    assert p['calibration']['diameter_m']==.0008
    assert p['holdout']['diameter_m']==pytest.approx(.00045)
    for d in p.values():
        assert np.all(np.diff(d['time_s'])>0)
        assert np.all(np.diff(d['fraction'])>=0)
        assert np.all((d['fraction']>=0)&(d['fraction']<=1))
    assert f['reading_sensitivity_time_s']==pytest.approx(120/1.3075)
    assert 'not experimental uncertainty' in f['reading_sensitivity_role']


@pytest.mark.parametrize('index',[0,1,2])
def test_actual_julia_replays_separate_carbon_case(index):
    case=reference()['cases'][index];p=case['inputs'];actual=column(**p)
    for key,value in case['result'].items():
        if value is None:
            assert actual[key] is None
        else:
            assert actual[key]==pytest.approx(value,rel=1e-9,abs=1e-13)
    assert actual['adsorbent_mass_kg']==pytest.approx(740*np.pi*.011**2/4*.1)
    assert actual['normalized_mass_residual']<1e-8
    series='calibration' if index==0 else 'holdout'
    assert compare(actual,series)==case['metrics']


def test_holdout_changes_rate_without_refitting_inventory_or_equilibrium():
    r=reference();a,b,c=[s['inputs'] for s in r['cases']]
    assert b['ldf_s_inv']/a['ldf_s_inv']==pytest.approx((.8/.45)**2)
    assert b['axial_dispersion_m2_s']==a['axial_dispersion_m2_s']
    assert c['axial_dispersion_m2_s']/a['axial_dispersion_m2_s']==pytest.approx(.45/.8)
    for key in ['capacity_kg_kg','affinity_m3_kg','porosity','particle_density_kg_m3',
                'superficial_velocity_m_s','area_m2','length_m','inlet_concentration_kg_m3']:
        assert a[key]==b[key]==c[key]
    assert r['verification']['status']=='passed_evidence_and_replay_checks'
    assert r['verification']['physical_comparison_passed']==r['physical_comparison_passed']


def test_missing_crossing_does_not_become_zero():
    assert crossing([0.,1.],[0.,.04],.05) is None
    assert crossing([0.,1.],[0.,.1],.05)==pytest.approx(.5)


def test_numerical_replay_does_not_promote_failed_physical_comparison():
    r=reference()
    assert r['numerical_refinement_passed'] is True
    assert r['calibration_passed'] is False
    assert r['physical_comparison_passed'] is False
    holdout=next(iter(r['metrics']['holdout'].values()))
    assert holdout['RMSE']>.10
    assert holdout['t50']['relative_error']>.20
