"""Check all retained evidence, calibration inputs and numerical acceptance limits."""
from collections import Counter

from examples.biocatalysis.pentyl_propanoate_history import run
from pore2process_demo.io import read_json,sha
from pore2process_demo.paths import ROOT

REFERENCE=run.HERE/"reference"


def check(reference=REFERENCE):
    run.check_previous()
    s=read_json(reference/"summary.json"); plan=read_json(run.PLAN)
    for name,digest in s["source_sha256"].items():
        if sha(ROOT/name)!=digest:
            raise ValueError(f"Changed catalyst-history source: {name}")
    for name,digest in s["retained_sha256"].items():
        if sha(reference/name)!=digest:
            raise ValueError(f"Changed retained catalyst-history result: {name}")
    assert s["actual_mass_kg"]==.006 and s["history_count"]==560 and s["primary_history_count"]==280
    assert s["condition_count"]==40 and s["snapshot_count"]==560*len(plan["times_min"])
    limits=plan["numerics"]; n=s["numerical"]
    for scope in ("grid","time"):
        for k,limit in run.LIMIT_KEYS.items():
            assert n[f"{scope}_{k}"]<=limits[limit]
    for key,limit in (("max_stage_residual","max_stage_residual_fraction"),
        ("max_normalized_balance_error","max_normalized_cumulative_balance_error"),
        ("independent_normalized_balance_error","max_normalized_cumulative_balance_error"),
        ("max_film_error_mol_L_min","max_film_residual_mol_L_min"),
        ("max_positive_film_slope_ratio","max_sampled_film_source_slope_over_capacity"),
        ("sampled_global_film_slope_ratio","max_sampled_film_source_slope_over_capacity")):
        assert n[key]<=limits[limit]
    assert n["all_40_clock_free_180min_vs_steady_max_difference"]<=limits["max_sample_conversion_change_fraction"]
    assert s["uniform_activity_bound_check"]["maximum_activity_violation"]<=1e-10
    data=run.read_rows(reference/"bed_comparison.csv"); history=run.read_rows(reference/"history.csv")
    assert len(data)==560 and all(float(r["actual_mass_kg"])==.006 for r in data)
    assert len(history)==560*len(plan["times_min"])
    assert set(Counter(r["task_id"] for r in history).values())=={len(plan["times_min"])}
    assert set(Counter((r["scenario"],r["law"],r["feed_ratio"]) for r in data).values())=={20}
    for d in s["score_scenarios"].values():
        for score in d["scores"].values():
            assert score["joint_40"]["count"]==40
    replay=read_json(reference/"replay.json")
    assert replay["repeated_histories"]==36 and replay["all_refinement_decisions_identical"]
    assert replay["combined_summary_sha256"]==s["full_summary_sha256"]
    assert replay["max_absolute_numerical_scalar_difference"]<=1e-11
    assert max(replay["max_absolute_difference_by_array"].values())<=1e-11
    assert replay["replay_source_sha256"]==sha(run.HERE/"replay.py")


if __name__=="__main__":
    check()
    print("Spatial catalyst history: all frozen inputs, holdout coverage, numerical limits and replay hashes verified.")
