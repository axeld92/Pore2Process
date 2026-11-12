"""Independent checks of the diagnostic initial profile, not a flow solver."""
import numpy as np
import pytest
import json

from examples.adsorption.transport_root_cause.initial_compatibility.reference import evaluate
from examples.adsorption.transport_root_cause.initial_compatibility.assessment import ordinary

SIDE=.0013891725543141453


def test_profile_range_and_periodicity():
    rng=np.random.default_rng(128)
    p=rng.uniform(-SIDE,2*SIDE,(1000,3))
    value,_=evaluate(p,SIDE)
    assert np.min(value)>=.25-1e-14
    assert np.max(value)<=.75+1e-14
    for dim in range(3):
        q=p.copy();q[:,dim]+=SIDE
        other,_=evaluate(q,SIDE)
        assert np.max(abs(other-value))<2e-14


def test_surface_trace_and_zero_smooth_normal_derivative():
    radius=SIDE/np.cbrt(2*np.pi/(3*.4))/2
    buffer=.4*(SIDE/np.sqrt(2)-2*radius)
    directions=np.array([[1.,0,0],[1,2,3],[-1,-3,2],[0,1,-1]])
    directions/=np.sqrt(np.sum(directions**2,axis=1))[:,None]
    wall=radius*directions
    c,_=evaluate(wall,SIDE)
    assert np.max(abs(c-(.5+.25*np.sin(2*np.pi*wall[:,0]/SIDE))))<1e-14
    h=buffer*1e-4
    plus,_=evaluate(wall+h*directions,SIDE)
    minus,_=evaluate(wall-h*directions,SIDE)
    assert np.max(abs(plus-minus)/(2*h))*buffer/.25<1e-7


def test_buffer_transition_matches_original_value_and_derivative():
    radius=SIDE/np.cbrt(2*np.pi/(3*.4))/2
    buffer=.4*(SIDE/np.sqrt(2)-2*radius)
    n=np.array([1.,2.,3.]);n/=np.linalg.norm(n)
    h=buffer*1e-4
    points=np.array([(radius+buffer+d)*n for d in [-h,0,h]])
    value,_=evaluate(points,SIDE)
    original=.5+.25*np.sin(2*np.pi*points[:,0]/SIDE)
    assert value[1]==pytest.approx(original[1],abs=1e-14)
    assert abs((value[2]-value[0]-original[2]+original[0])/(2*h))*buffer/.25<1e-7


def test_buffer_is_fixed_physical_nonoverlap_and_bulk_is_unchanged():
    point=np.array([[.25*SIDE,.25*SIDE,.25*SIDE]])
    value,data=evaluate(point,SIDE)
    assert 0<data['buffer']<(SIDE/np.sqrt(2)-2*data['radius'])/2
    assert data['signed_distance'][0]>data['buffer']
    assert value[0]==pytest.approx(.75,abs=1e-14)


def test_deep_solid_centre_is_finite_but_not_a_wall_projection():
    value,data=evaluate([[0.,0.,0.]],SIDE)
    assert value[0]==.5
    assert data['signed_distance'][0]==pytest.approx(-data['radius'])


def test_typed_assessment_preserves_values_and_strict_json():
    original=dict(checks=dict(passed=np.bool_(True)),values=[np.float64(.125),np.int64(64)])
    value=ordinary(original)
    assert type(value['checks']['passed']) is bool
    assert json.loads(json.dumps(value,allow_nan=False))==dict(checks=dict(passed=True),values=[.125,64])
    with pytest.raises(ValueError):
        json.dumps(ordinary(dict(value=np.float64(np.nan))),allow_nan=False)
