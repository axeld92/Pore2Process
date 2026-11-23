import copy

import numpy as np
import pytest

from examples.biocatalysis.pentyl_propanoate_loading_evidence import run as r


def test_pentyl_nominal_concentrations_keep_acid_and_total_bases_separate():
    groups,basis,ref = r.pentyl_groups(r.read_json(r.PLAN))
    assert ref == pytest.approx(25.4428025)
    assert [len(groups[1,w]) for w in (2.5,5,10,20)] == [11,11,17,15]
    assert groups[1,20][0]["nominal_concentration_g_L"] == pytest.approx(80.527705)
    assert groups[1,20][0]["mass_g"] == pytest.approx(3.2211082)
    assert groups[1,20][0]["loading_pct_total_substrates"] < 10
    assert groups[2,20][0]["loading_pct_total_substrates"] < groups[1,20][0]["loading_pct_total_substrates"]
    assert basis[2]["nominal_liquid_density_g_L"] == pytest.approx(856.536125)


def test_reading_sensitivity_and_missing_latest_time():
    v, bracket, bounds = r.crossing_bounds([0.,5.,10.],[0.,.3,.6],[0.,.03,.03],.4)
    assert v == pytest.approx(20/3)
    assert bracket == [5,10]
    assert bounds[0] < v < bounds[1]
    assert r.crossing_bounds([0,5,10],[0,.3,.6],[0,.03,.03],.59)[2][1] is None


def test_external_figures_and_comparable_basis():
    plan=r.read_json(r.PLAN)
    groups,_,_=r.pentyl_groups(plan)
    e,i=r.all_progress(plan,groups)
    alternate,_=r.all_progress(plan,groups,5)
    at = lambda data,name,w,target: next(x for x in data[name] if x["loading_pct_reported"]==w
                                        and x["target_conversion"]==target)
    assert at(e,"richetti2010",20,.6)["relative_activity_per_gram"] < .6
    assert at(e,"paroul2010",10,.7)["time_min"] < at(alternate,"paroul2010",10,.7)["time_min"]
    assert all(x["loading_pct_total_substrates"] is None for x in e["kuperkar2014"])
    assert all(x["bracket_low_min"]>0 for rows in e.values() for x in rows)
    assert at(i,"lorenzo2026_1to1",20,.4)["relative_activity_per_gram"] > 1
    direct=r.direct_exposure_checks(plan)
    assert direct[0]["higher_exposure_pct_min"] == direct[0]["anchor_exposure_pct_min"]
    assert direct[0]["conversion_gap_pp"] == 46
    assert direct[0]["minimum_gap_with_reported_deviation_scenario_pp"] > 25
    for direction in (True,False):
        shifted = r.reading_progress_scenario(e,direction)
        assert len(shifted["richetti2010"])==len(e["richetti2010"])
        assert at(shifted,"richetti2010",10,.6)["time_min"]==at(e,"richetti2010",10,.6)["time_min"]


def test_crossing_predictions_do_not_read_held_out_response():
    # Analytic stand-in for the Julia total-rate model: this test isolates
    # the data split and intercept construction, not the rate implementation.
    class Model:
        @staticmethod
        def log_total_rate(law,shape,u):
            return shape*np.log(u)
    rows=[{"series":"example","target_conversion":.3,"relative_loading_for_model":u,
           "time_min":8*u**(-.8)} for u in (.5,1.,2.)]
    fit={"law":"power","shape":.8}
    anchors,test=rows[:1],rows[1:]
    p=r.crossing_predict(Model(),fit,anchors,test)
    changed=copy.deepcopy(test)
    changed[0]["time_min"]*=4
    q=r.crossing_predict(Model(),fit,anchors,changed)
    assert [v["predicted_time_min"] for v in p] == [v["predicted_time_min"] for v in q]
    assert r.log_time_scores(p)["log_time_rmse"] < 1e-14
