"""Independent convective-column limits and rejection of unqualified handoffs."""
import json

import numpy as np
import pytest

from pore2process_demo.paths import ROOT
from examples.adsorption.column_transport.bridge import (
    boundary_moments, qualified_dispersion, adsorption_column, REQUIRED_SCREENS)
from examples.adsorption.film_particle.bridge import bed


@pytest.mark.parametrize('Pe',[.05,2.,20.,1000.])
def test_inlet_condition_changes_mean_without_changing_liquid_inventory(Pe):
    p=dict(length_m=.1,velocity_m_s=.001,dispersion_m2_s=.0001/Pe)
    fixed=boundary_moments('dirichlet_neumann',**p)
    flux=boundary_moments('danckwerts',**p)
    tau=p['length_m']/p['velocity_m_s']
    assert flux['mean_s']==pytest.approx(tau)
    assert fixed['mean_s']==pytest.approx(tau*(1-(-np.expm1(-Pe))/Pe))
    assert fixed['mean_s']<flux['mean_s']
    assert flux['variance_s2']==pytest.approx(tau**2*(2/Pe-2*(-np.expm1(-Pe))/Pe**2))


@pytest.mark.parametrize('Pe',[2.,20.])
def test_frozen_particle_column_nonadsorbing_limit_has_physical_void_time(Pe):
    # An independent numerical-component test, not a native material closure.
    L=.1;eps=.4;A=.001;u=.0004;v=u/eps;tau=eps*L/u
    p=dict(length_m=L,area_m2=A,porosity=eps,superficial_velocity_m_s=u,
        inlet_concentration_kg_m3=1.,particle_diameter_m=.0008,
        particle_density_kg_m3=1291.44851658,diffusivity_m2_s=0.,film_m_s=0.,
        capacity_kg_kg=0.,affinity_m3_kg=0.,axial_dispersion_m2_s=v*L/Pe,
        active_inventory_fraction=1.,radial_cells=2,t_end_s=12*tau,dt_s=tau/2000,
        initial_concentration_kg_m3=0.,initial_loading_kg_kg=0.,linear_partition_m3_kg=None)
    exact=boundary_moments('danckwerts',length_m=L,velocity_m_s=v,dispersion_m2_s=v*L/Pe)
    errors=[]
    for n in (400,800):
        r=bed(**p,axial_cells=n)
        t=np.array(r['time_s']);survival=1-np.array(r['outlet_concentration_kg_m3'])
        mean=float(np.trapezoid(survival,t))
        variance=float(2*np.trapezoid(t*survival,t)-mean**2)
        assert mean==pytest.approx(tau,rel=.001)
        assert r['normalized_mass_residual']<1e-8
        assert r['physical_adsorbent_mass_kg']==pytest.approx((1-eps)*p['particle_density_kg_m3']*A*L)
        assert r['active_adsorbent_mass_kg']==r['physical_adsorbent_mass_kg']
        assert r['fluid_volume_m3'][0]==pytest.approx(eps*A*L)
        assert max(abs(np.array(r['adsorbed_inventory_kg'])))==0
        errors.append(abs(variance/exact['variance_s2']-1))
    assert errors[1]<errors[0]
    assert errors[1]<.025


def unit_closure():
    # A contract fixture, explicitly not evidence that a native run qualified.
    return dict(boundary='convective-column',native_run_id='unit-contract-fixture',
        calibration_boundary='native-inert-tracer-only',axial_dispersion_m2_s=1e-6,
        checks={k:True for k in REQUIRED_SCREENS})


@pytest.mark.parametrize('screen',REQUIRED_SCREENS)
def test_failed_native_screen_cannot_supply_adsorption_dispersion(screen):
    c=unit_closure();c['checks'][screen]=False
    with pytest.raises(ValueError,match='unqualified'):
        adsorption_column(c)


def test_absorbing_closure_and_per_curve_fits_are_rejected():
    c=unit_closure();c['boundary']='absorbing-end'
    with pytest.raises(ValueError,match='convective-column'):qualified_dispersion(c)
    c=unit_closure();c['calibration_boundary']='fit-each-column'
    with pytest.raises(ValueError,match='independently'):qualified_dispersion(c)
    c=unit_closure()
    with pytest.raises(ValueError,match='inventory'):adsorption_column(c,active_inventory_fraction=.5)


@pytest.mark.parametrize('bad',[None,float('nan'),float('inf'),-1.,True])
def test_invalid_or_unavailable_transport_is_not_a_valid_coefficient(bad):
    c=unit_closure();c['axial_dispersion_m2_s']=bad
    with pytest.raises(ValueError,match='positive finite'):qualified_dispersion(c)


def test_actual_failed_native_study_does_not_become_an_adsorption_fit():
    path=ROOT/'examples/adsorption/column_transport/reference/transport_reference.json'
    s=json.loads(path.read_text())
    assert s['native_runs']==19
    assert s['transport_qualified'] is False
    assert s['closure']['axial_dispersion_m2_s'] is None
    assert s['unqualified_handoff_rejected'] is True
    with pytest.raises(ValueError,match='unqualified'):
        adsorption_column(s['closure'])
