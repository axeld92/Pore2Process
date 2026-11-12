"""Mechanical guards for native fixed-flow refinement; no surrogate solver."""
from copy import deepcopy

import numpy as np
import pytest

from examples.adsorption.transport_root_cause.hydraulics import run


def accepted_analysis():
    return dict(superficial_velocity_m_s=.008,
        plane_mass_relative_spread=1e-8,
        stationarity=dict(superficial_velocity_m_s=1e-9,mean_lattice_energy=1e-9),
        voxel_porosity=.6,
        raw=dict(max_Mach=.02,max_density_departure=1e-6,
            resolution=64,unit_cell_side_m=.0014))


def test_native_checks_do_not_accept_velocity_rescaling():
    a=accepted_analysis()
    assert all(run.native_checks(a,.008).values())
    a['superficial_velocity_m_s']*=1.001
    checks=run.native_checks(a,.008)
    assert not checks['flow']
    assert all(value for key,value in checks.items() if key!='flow')


def test_independent_stationarity_and_geometry_guards():
    a=accepted_analysis()
    a['stationarity']['mean_lattice_energy']=2e-7
    a['voxel_porosity']=.58
    checks=run.native_checks(a,.008)
    assert not checks['stationarity']
    assert not checks['porosity']
    assert checks['flow']


def test_source_point_carrier_scale_and_flow_contract(tmp_path,monkeypatch):
    a=accepted_analysis()
    monkeypatch.setattr(run,'analyze',lambda _:deepcopy(a))
    field=tmp_path/'field.csv'
    field.write_text('i,j,k,ux_m_s,uy_m_s,uz_m_s\n0,0,0,0.01,0,0\n1,0,0,0.03,0.04,0\n')
    (tmp_path/'raw.json').write_text('{}\n')
    carrier=dict(native_target_velocity_m_s=.008,carrier_viscosity_Pa_s=.00096,
        carrier_density_kg_m3=1000.,area_m2=1e-4)
    point=run.source_point(tmp_path,carrier)
    assert point['n']==64
    assert point['velocity_scale']==pytest.approx(.96)
    assert point['speed_m_s']==pytest.approx(.008*.96/.6)
    assert point['max_speed_m_s']==pytest.approx(.05*.96)
    assert point['actual_flow_m3_s']==pytest.approx(.008*.96*1e-4)
    assert point['field']==str(field.resolve())
    assert len(point['field_sha256'])==64


def test_source_point_rejects_unqualified_native_field(tmp_path,monkeypatch):
    a=accepted_analysis()
    a['raw']['max_Mach']=.1
    monkeypatch.setattr(run,'analyze',lambda _:a)
    with pytest.raises(ValueError,match='declared screens'):
        run.source_point(tmp_path,dict(native_target_velocity_m_s=.008))
