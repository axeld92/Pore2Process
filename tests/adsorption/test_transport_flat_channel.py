import numpy as np
import pytest

from examples.adsorption.transport_root_cause.flux_audit.flat_channel.reference import (
    mode, concentration, boundary_checks, converged_reference,
)


PHYSICAL=dict(height_m=.0013891725543141453,mean_speed_m_s=.012937295182166382,
              diffusivity_m2_s=1e-9,duration_s=.21478314576454161)


def test_neumann_oracle_zero_flow_matches_analytical_sine():
    p={**PHYSICAL,'mean_speed_m_s':0.}
    x=(np.arange(32)+.5)*p['height_m']/32;y=x[::-1]
    got=concentration(x,y,points=64,steps=32,**p)
    expected=.5+.25*np.exp(-p['diffusivity_m2_s']*(2*np.pi/p['height_m'])**2*p['duration_s'])*np.sin(2*np.pi*x/p['height_m'])
    np.testing.assert_allclose(got,expected,atol=1e-13,rtol=0.)


def test_neumann_oracle_initial_state():
    x=np.array([.1,.3,.7])*PHYSICAL['height_m'];y=np.array([0.,.4,1.])*PHYSICAL['height_m']
    got=concentration(x,y,points=64,steps=1,**{**PHYSICAL,'duration_s':0.})
    np.testing.assert_allclose(got,.5+.25*np.sin(2*np.pi*x/PHYSICAL['height_m']),atol=1e-13)


def test_even_extension_enforces_wall_neumann_condition():
    f=mode(points=128,steps=256,**PHYSICAL);checks=boundary_checks(f,PHYSICAL['height_m'])
    assert checks['even_reflection_error']<1e-12
    assert checks['max_wall_normal_mode_derivative_times_height']<1e-10


def test_independent_space_and_time_refinements_qualify():
    y=(np.arange(32)+.5)*PHYSICAL['height_m']/32
    _,record=converged_reference(y,y,**PHYSICAL)
    assert record['passed']


@pytest.mark.parametrize('change',[{'height_m':0.},{'mean_speed_m_s':-1.},{'diffusivity_m2_s':0.}])
def test_invalid_reference_inputs_rejected(change):
    with pytest.raises(ValueError):mode(**{**PHYSICAL,**change})


@pytest.mark.parametrize('n',[32,64,96])
def test_finer_reference_qualifies_near_wall_sample_locations(n):
    from examples.adsorption.transport_root_cause.flux_audit.flat_channel.reference_qualified import converged_reference as qualified
    y=(np.arange(n)+.5)*PHYSICAL['height_m']/n
    _,record=qualified(y,y,**PHYSICAL)
    assert record['passed']
