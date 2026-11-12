"""Independent shear oracle checks; native executions are separately archived."""
import numpy as np
import pytest
from examples.adsorption.transport_root_cause.channel.reference import mode,sample_mode,converged_reference
from examples.adsorption.transport_root_cause.channel.run import specifications,PHYSICAL


def test_uniform_mode_matches_exact_translation_and_diffusion():
    p=dict(side_m=.0014,speed_m_s=.013,shear_fraction=0.,diffusivity_m2_s=1e-9,duration_s=.215)
    f=mode(**p,points=64,steps=127);k=2*np.pi/p['side_m']
    exact=np.exp((-p['diffusivity_m2_s']*k*k-1j*k*p['speed_m_s'])*p['duration_s'])
    np.testing.assert_allclose(f,exact,atol=2e-14,rtol=2e-14)


def test_fourier_interpolation_preserves_cell_centred_phase():
    count=32;side=.004;y=(np.arange(count)+.5)*side/count
    f=np.exp(4j*np.pi*y/side);sample=np.linspace(.0001,.0039,47)
    np.testing.assert_allclose(sample_mode(f,sample,side),np.exp(4j*np.pi*sample/side),atol=1e-14)


def test_shear_reference_tends_to_exact_characteristic_transport():
    p=dict(side_m=.0014,speed_m_s=.013,shear_fraction=.75,diffusivity_m2_s=1e-16,duration_s=.215)
    count=128;y=(np.arange(count)+.5)*p['side_m']/count;k=2*np.pi/p['side_m']
    exact=np.exp(-1j*k*p['speed_m_s']*(1+.75*np.cos(k*y))*p['duration_s'])
    numerical=mode(**p,points=count,steps=512)
    np.testing.assert_allclose(numerical,exact,atol=5e-8,rtol=0)


def test_shear_reference_converges_and_respects_continuum_bounds():
    side=PHYSICAL['side_m'];x,y=np.meshgrid((np.arange(24)+.5)*side/24,(np.arange(24)+.5)*side/24,indexing='ij')
    exact,checks=converged_reference(x.ravel(),y.ravel(),side_m=side,speed_m_s=PHYSICAL['speed_m_s'],
        shear_fraction=.75,diffusivity_m2_s=1e-9,duration_s=PHYSICAL['end_s'])
    assert checks['passed']
    assert exact.min()>=.25 and exact.max()<=.75
    assert abs(exact.mean()-.5)<1e-14


def test_specifications_share_time_step_and_change_only_declared_controls():
    cases=specifications();shear=[s for s in cases if s['id'] in ('shear32','shear48','shear64')]
    assert [s['n'] for s in shear]==[32,48,64]
    assert all(s['shear_fraction']==.75 and s['amplitude']==.25 for s in shear)
    assert next(s for s in cases if s['id']=='constant_shear32')['amplitude']==0
    assert PHYSICAL['diffusivity_m2_s']==1e-9


def test_invalid_reference_parameters_fail():
    with pytest.raises(ValueError):
        mode(side_m=0,speed_m_s=.01,shear_fraction=.75,diffusivity_m2_s=1e-9,duration_s=.1)


def test_thin_controls_match_physical_time_and_bound_lattice_speed():
    from examples.adsorption.transport_root_cause.channel.thin.run import specifications as thin_specifications,DURATION
    for spec in thin_specifications():
        dt=PHYSICAL['dt_s']*spec['dt_factor'];steps=round(DURATION/dt)
        assert abs(steps*dt-DURATION)<1e-15
        assert abs(DURATION/dt-steps)<1e-6
        h=PHYSICAL['side_m']/spec['n']
        assert PHYSICAL['speed_m_s']*1.75*dt/h<.05
