import copy

import numpy as np
import pytest

from examples.biocatalysis.pentyl_propanoate_loading.run import (
    PLAN, batch_rows, crossing_time, curve_alignment_check, liquid_sample_screen,
    matched_conversion, qualification, time_based_power, vectors,
)
from examples.biocatalysis.pentyl_propanoate.run import BATCH, read_rows
from pore2process_demo.io import read_json


def test_nominal_source_exposures_and_loadings():
    plan = read_json(PLAN)
    groups = batch_rows(read_rows(BATCH), plan)
    assert [len(rows) for rows in groups.values()] == [11, 11, 16, 19]
    assert groups[10][0]["mass_g"] == pytest.approx(1.0177121)
    assert groups[20][0]["nominal_concentration_g_L"] == pytest.approx(50.885605)
    assert groups[10][1]["exposure_kg_s_mol_acid"] == pytest.approx(0.1*0.0740785*5.3*60)
    assert groups[10][0]["printed_ratio"] == pytest.approx(.656/.344)
    _, train = vectors(groups, [2.5, 5., 20.], "printed", "nominal_ratio", True)
    assert 10. not in train[0]
    totals = [sum(train[-1][train[0] == loading]) for loading in (2.5, 5., 20.)]
    assert totals == pytest.approx([len(train[0])/3]*3)


def test_crossing_uses_increasing_bracket_without_extrapolation():
    value, interval = crossing_time([0, 2, 4, 8], [0, .2, .2, .6], .4)
    assert value == pytest.approx(6.)
    assert interval == [4., 8.]
    assert crossing_time([0, 2, 4], [0, .2, .1], .2)[0] == 2
    with pytest.raises(ValueError):
        crossing_time([0, 2], [.0, .1], .2)


def test_matched_progress_and_rate_proxy():
    plan = read_json(PLAN)
    groups = batch_rows(read_rows(BATCH), plan)
    records, summary = matched_conversion(groups, plan)
    assert len(records) == 48
    early = summary["printed"]["0.1"]
    later = summary["printed"]["0.4"]
    assert early["time_min_by_loading"]["10"] == pytest.approx(10.6166667)
    assert later["descriptive_total_rate_power_exponent"] > early["descriptive_total_rate_power_exponent"]
    assert later["time_based_speedup_10_to_20pct"] == pytest.approx(43.048/24.2293706)


def test_shape_free_holdout_uses_no_withheld_response():
    plan = read_json(PLAN)
    alpha = .82
    common_t = np.array([0., 5., 10., 20., 35., 60.])
    observed = [.0, .1, .2, .3, .4, .5]
    groups = {w: [{"time_min": t/(w/10)**alpha, "printed": y} for t,y in zip(common_t, observed)]
              for w in (2.5, 5., 10., 20.)}
    assert time_based_power(groups, plan["matched_conversion_levels"]) == pytest.approx(alpha)
    before, report = curve_alignment_check(groups,plan)
    assert report["equal_loading_rmse_pp_on_common_support"]["power"] < 1e-10
    modified = copy.deepcopy(groups)
    modified[10.][2]["printed"] += .02
    after, _ = curve_alignment_check(modified,plan)
    # Only scoring of this fold may change; its predictions cannot.
    fields = ("training_power_exponent", "proportional", "power")
    target = lambda rows: [[row[key] for key in fields] for row in rows if row["held_out_loading_pct"] == 10.]
    assert target(before) == target(after)


def test_qualification_requires_mean_improvement_and_no_bad_fold():
    plan = read_json(PLAN)
    scores = {"a": {"proportional": {"rmse_percentage_points": 2.}, "power": {"rmse_percentage_points": 1.}},
              "b": {"proportional": {"rmse_percentage_points": 3.}, "power": {"rmse_percentage_points": 2.}}}
    assert qualification(scores,"proportional","power",plan)["passes_screen"] is True
    scores["a"]["power"]["rmse_percentage_points"] = 2.3
    assert qualification(scores,"proportional","power",plan)["passes_screen"] is False


def test_liquid_sampling_inventory_sensitivity():
    plan = read_json(PLAN)
    result = liquid_sample_screen({10: [{"time_min": t} for t in (0., 10., 20., 30.)]},plan)["10"]
    assert result["max_volume_removed_before_final_sample_fraction"] == pytest.approx(.005)
    assert result["cumulative_proportional_exposure_increase_factor"] == pytest.approx(
        (1 + 1/.9975 + 1/.995)/3)
