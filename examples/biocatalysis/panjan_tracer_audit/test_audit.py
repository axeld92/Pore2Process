import numpy as np

from . import audit


def test_published_traces_match_at_raw_markers_but_R2_is_distinct():
    d=audit.duplicate_audit()
    assert d["blank"]["count_core"]>=10 and d["R1"]["count_core"]>=15
    assert d["blank"]["all_core_within_reading_allowance"]
    assert d["R1"]["all_core_within_reading_allowance"]
    assert not d["R2"]["all_core_within_reading_allowance"]


def test_last_rising_crossing_ignores_early_baseline_excursion():
    t=np.arange(7,dtype=float);c=np.array([0.,.2,0.,0.,.2,.4,0.])
    assert audit.last_rising_crossing(t,c,.1)==3.5


def test_single_exponential_necessary_bound_on_independent_reference():
    t=np.arange(0,31,1.); c=np.exp(-t/25.)
    r=audit.exponential_fall_bound(t,c,0.,.2,10.)
    np.testing.assert_allclose(r["tau_upper_s"],25.,rtol=1e-12)
    d=audit.detector_audit()["first_order_decay_necessary_bound"]["blank"]
    assert not d["compatible_with_reported_T90_as_single_exponential"]
    assert d["minimum_later_if_reported_first_order_mM"]>d["later_upper_mM"]


def test_visible_wide_window_is_not_a_complete_RTD():
    m=audit.moment_scenarios()
    for scenario in m.values():
        for name in audit.WIDE:
            r=scenario["wide_partial"][name]
            assert not r["full_moments_available"] and r["last_mM"]>.3
        assert scenario["wide_partial"]["R4"]["partial_glucose_umol"]>.5


def test_rescaling_weighted_Trt_kernel_preserves_mass_and_mean():
    k=np.array([0.,.1,.5,.4]);dt=.05
    v=audit.scale_kernel(k,1.7,dt)
    np.testing.assert_allclose(v.sum(),1.)
    np.testing.assert_allclose(audit.kernel_moments(v,dt)["mean_s"],1.7*audit.kernel_moments(k,dt)["mean_s"])


def test_predictions_do_not_access_thesis_reactor_concentrations(monkeypatch):
    def fake_pulse(name,bg,dt=.05):
        assert name=="blank"
        t=np.arange(81)*dt
        return t,np.maximum(0.,1-abs(t-2))
    monkeypatch.setattr(audit,"pulse",fake_pulse)
    old={"A":{"model_parameters":{"extra_mixed_volume_mean_s":0.,"rare":{"available":False}}},
         "independent_storage_bounds":{"A":{"bead_envelope_volume_ul":12.7}}}
    out,records=audit.predict_before_score(np.array([0.,.3,.7]),old)
    assert "R2__prior_trt_independent__bg0" in out
    for name in ("R1","R2"):
        np.testing.assert_allclose(records[name]["kernel_moments"]["prior_trt_independent"]["mean_s"],audit.geometry(name)["nominal_time_s"])


def test_archived_real_A_and_volume_feasibility():
    k,_,old=audit.prior_kernel()
    assert old["A"]["numerical_qualification"] and k.sum()>0.9999
    m=audit.moment_scenarios()
    models={name:{"kernel_moments":{"prior_trt_independent":{"variance_s2":10.}}} for name in ("R1","R2")}
    v=audit.volume_audit(m,models,old)["scenarios"]
    for s in v.values():
        assert not s["common_positive_delay"]["feasible"]
        assert not s["bead_volume_scaled_positive_delay"]["feasible"]
        assert s["narrow"]["R2"]["nominal_mean_s"]==21.
        assert s["narrow"]["R2"]["apparent_extra_hold_up_ul"]>15.
