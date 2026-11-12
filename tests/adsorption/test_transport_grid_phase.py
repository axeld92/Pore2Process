import numpy as np
import pytest
from examples.adsorption.transport_root_cause.grid_phase.geometry import sphere_distance


def test_periodic_translation_preserves_continuum_surface_gap_and_normal():
    side=.0013891725543141453;delta=side/64
    points=side*np.array([[.99,.12,.21],[.1,.2,.3],[.4,.62,.87]])
    old=sphere_distance(points,side_m=side)
    shifted=points.copy();shifted[:,0]=(shifted[:,0]+delta)%side
    new=sphere_distance(shifted,side_m=side,phase_x_m=delta)
    np.testing.assert_allclose(old[0],new[0],atol=1e-18)
    np.testing.assert_allclose(old[1],new[1],atol=2e-15)


def test_analytical_sphere_pole_and_periodic_image():
    side=.0014;radius=.0004;gap=.000006;delta=side/64
    points=np.array([[delta+radius+gap,0,0],[delta-radius-gap+side,0,0]])
    distance,normal,_=sphere_distance(points,side_m=side,phase_x_m=delta)
    np.testing.assert_allclose(distance,gap,atol=1e-18)
    np.testing.assert_allclose(normal,[[1,0,0],[-1,0,0]],atol=1e-15)


def test_nonfinite_surface_inputs_fail():
    with pytest.raises(ValueError):sphere_distance([[float('nan'),0,0]],side_m=.0014)
