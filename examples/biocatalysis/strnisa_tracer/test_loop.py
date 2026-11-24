import tomllib

import numpy as np

from .loop import convolve, exponential_connection, measured, moments, rare_delay, tracer_kernel
from .native import bead_centers, config, dimensions, voxel_mask
from .handoff import matched_dt
from . import compare, loop
from .source_audit import audit
from . import storage_audit


def test_figure_readings_and_mass():
    t, b = measured("blank")
    mb = moments(t,b)
    a=moments(*measured("A"))
    assert abs(mb["observed_glucose_umol"]-.5)<.02
    assert a["mean_s"]>mb["mean_s"]
    assert a["variance_s2"]>mb["variance_s2"]
    np.testing.assert_allclose(t[-1],400.)
    assert 0.18 < (a["mean_s"]-mb["mean_s"])/60 < .22


def test_normalized_figure_a_matches_fig3b_full_loop():
    r=audit()
    assert abs(r["inferred_time_scale_s"]/r["full_measured_loop_mean_s"]-1)<.01
    assert r["normalized_time_max_error_using_full_loop_mean"]<.02
    assert r["normalized_concentration_max_error"]<.015
    assert r["normalized_time_max_error_using_printed_bed_mean"]>10


def test_convolution_preserves_dose_and_adds_moments():
    dt=.05
    t=np.arange(401)*dt
    b=np.maximum(0,1-abs(t-5))
    k=tracer_kernel(np.array([1.2,2.3,3.4]),dt)
    tp,pred=convolve(b,k,dt)
    a=moments(t,b)
    p=moments(tp,pred)
    np.testing.assert_allclose(p["observed_glucose_umol"],a["observed_glucose_umol"],rtol=1e-13)
    np.testing.assert_allclose(p["mean_s"]-a["mean_s"],2.3,rtol=1e-13)
    np.testing.assert_allclose(p["variance_s2"]-a["variance_s2"],np.var([1.2,2.3,3.4]),rtol=1e-12)
    # A censored marker must not silently be converted into recovered mass.
    k=tracer_kernel(np.array([1.2,2.3,np.nan]),dt)
    np.testing.assert_allclose(k.sum(),2/3)


def test_added_mixed_volume_moments():
    dt=.005
    k=exponential_connection(2.,dt)
    t=np.arange(len(k))*dt
    assert np.isclose(k.sum(),1)
    np.testing.assert_allclose(k@t,2.,rtol=1e-5)
    np.testing.assert_allclose(k@(t-2.)**2,4.,rtol=1e-5)


def test_rare_delay_adds_mean_and_variance_without_losing_mass():
    dt=.005
    for rate in (.045,.135):
        k=rare_delay(rate,40.,dt)
        t=np.arange(len(k))*dt
        m=rate*40.
        np.testing.assert_allclose(k.sum(),1.,rtol=1e-13)
        np.testing.assert_allclose(k@t,m,rtol=1e-6)
        np.testing.assert_allclose(k@(t-m)**2,2*rate*40**2,rtol=1e-6)


def test_long_prediction_is_forward_only(monkeypatch):
    plan={"protocol":config(),"reactors":{
        "A":{"extent_m":[.02,.003,.0004],"geometry":{"continuous_porosity":.57}},
        "C":{"extent_m":[.08,.003,.0004],"geometry":{"continuous_porosity":.59}}}}
    def fake_measurement(name,dt,bg):
        # Any attempt to access either reactor's measurements fails.
        assert name=="blank"
        t=np.arange(401)*dt
        return t,np.maximum(0.,1-abs(t-5.))
    monkeypatch.setattr(loop,"measured",fake_measurement)
    params={"extra_mixed_volume_mean_s":.5,
            "rare":{"available":True,"rate_A":.05,"mean_event_s":4.}}
    times=np.array([10.,12.,14.])
    pred=compare.predict_curves(times,"C",plan,params)
    def mean(k):
        return moments(np.arange(len(pred[k]))*compare.DT,pred[k])["mean_s"]
    assert abs(mean("measured_blank_bg0")-17.)<1e-10
    assert abs(mean("common_rare_delay_bg0")-17.2)<1e-5
    assert abs(mean("solid_volume_rare_delay_bg0")-
               (17.+.2*compare.solid_volume_ratio(plan,"C")))<1e-5


def test_accessible_storage_mean_and_variance_signs():
    p=config()
    plan={"protocol":p,"reactors":{"A":{"extent_m":[.02,.003,.0004],
              "geometry":{"continuous_porosity":.57}}}}
    s=dict(incremental_observed=dict(mean_s=8.,variance_s2=100.),
           incremental_predicted=dict(mean_s=9.,variance_s2=9.))
    r=storage_audit.check(plan,"A",s,storage_audit.parameters())
    assert not r["can_positive_storage_explain_mean"]
    assert r["implied_accessible_fraction_of_bead_envelope"]<0
    a,b=[r["cases"][name] for name in ("dry_sorption","aqueous_tracers")]
    assert 0<a["additional_mean_s"]<b["additional_mean_s"]
    assert 0<a["additional_variance_s2"]<b["additional_variance_s2"]
    np.testing.assert_allclose(a["additional_variance_s2"],2*a["additional_mean_s"]*r["diffusion_time_s"])


def test_packing_shape_and_time_step():
    p=config();g=p["geometry"]
    extent,e=dimensions("A")
    beads,r=bead_centers(extent,e,g["diameter_um"]*1e-6,4,2,.25,.002)
    assert abs(r["continuous_porosity"]-e)<.001
    assert r["minimum_gap_m"]>=g["diameter_um"]*1e-6*.002-1e-12
    h=extent[2]/12
    shape=tuple(np.rint(extent/h).astype(int))
    mask=voxel_mask(beads,shape,h)
    assert mask.shape==(shape[0],shape[1]+2,shape[2]+2)
    assert mask[0,1:-1,1:-1].all() and mask[-1,1:-1,1:-1].all()
    assert not mask[:,0,:].any() and not mask[:,:,-1].any()
    dt=matched_dt(p,20)
    assert 12*np.sqrt(2*p["fluid"]["glucose_diffusivity_m2_s"]*1.2*dt)/(.0004/20)<1
