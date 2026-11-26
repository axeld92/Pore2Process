import math

import numpy as np
import pytest

from examples.adsorption.breakthrough_hypotheses.hydrodynamics import handoff, native
from .report import finite_pulse_bins, reference_bins, sphere_check
from .run import settings


@pytest.mark.parametrize("beads,diameter,period",[("fine",.001,5.25),("coarse",.004,7)])
def test_physical_flow_and_similarity(beads,diameter,period):
    case = settings(beads,period,31,8)
    speed = case["pore_speed_m_s"]
    assert speed*60_000 == pytest.approx(1.34350485,rel=1e-7)
    ratio = native.p.DIAMETER/diameter
    nu = native.p.VISCOSITY/native.p.DENSITY
    q = case["flow_m3_s"]*ratio*nu/case["nu_m2_s"]
    u = q/(math.pi*(ratio*case["tube_radius_m"])**2*case["pore_fraction"])
    assert u*native.p.DIAMETER/nu == pytest.approx(case["pore_Re"])
    assert speed*diameter/case["D_m2_s"] == pytest.approx(case["pore_Pe"])


def test_integrated_finite_pulse_and_reference_variance():
    # Unit-length pulse overlaps collection bins by different fractions.
    edges = np.arange(0.,5.)
    signal = finite_pulse_bins(np.array([.5,2.]),edges,1.)
    np.testing.assert_allclose(signal,[.25,.25,.5,0.])
    assert signal @ np.diff(edges) == pytest.approx(1.)
    edges = np.arange(0,60000,10.)
    mean, length, alpha, width = 8900.,200.,9.6,300.
    x = (edges[:-1]+edges[1:])/2
    mass = reference_bins(edges,width,mean,length,alpha)*np.diff(edges)/width
    assert mass.sum() == pytest.approx(1.,abs=1e-8)
    assert x @ mass == pytest.approx(mean+width/2,rel=1e-7)
    assert (x-x@mass)**2 @ mass == pytest.approx(2*alpha/length*mean**2+width**2/12,rel=1e-5)


def test_periodic_sphere_audit(tmp_path):
    path = tmp_path/"centres.tsv"
    np.savetxt(path,[[0,0,0],[2,1,0]])
    result = sphere_check(path,4,2)
    assert result["min_centre_distance_in_beads"] == pytest.approx(math.sqrt(5))
    overlap = tmp_path/"overlap.tsv"
    np.savetxt(overlap,[[.1,0,0],[3.9,0,0]])
    with pytest.raises(ValueError,match="overlap"):
        sphere_check(overlap,4,2)


def test_connected_component_and_historical_exclusion_limit():
    mask = np.zeros((4,6,6),dtype=np.uint8)
    mask[:,1,1] = mask[1,4,4] = 1
    ijk = np.argwhere(mask)
    velocity = np.zeros((len(ijk),3))
    velocity[0,0] = 1.
    assert handoff.connected_pore_component(ijk,mask,velocity).sum() == 4
    with pytest.raises(ValueError,match="0.1%"):
        handoff.connected_pores(ijk,mask,velocity)
