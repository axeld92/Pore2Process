"""Verify the dual-feed result counts, saved science boundary and immutable inputs."""
from collections import Counter

from examples.biocatalysis.pentyl_propanoate.run import read_rows
from examples.biocatalysis.pentyl_propanoate_equimolar.run import HERE,PLAN,check_prior
from pore2process_demo.io import read_json,sha
from pore2process_demo.paths import ROOT

REFERENCE=HERE/"reference"


def check():
    check_prior()
    record=read_json(REFERENCE/"summary.json")
    plan=read_json(PLAN)
    for name,digest in record["source_sha256"].items():
        if sha(ROOT/name)!=digest:
            raise ValueError(f"Changed dual-feed input: {name}")
    for name,digest in record["retained_sha256"].items():
        if sha(REFERENCE/name)!=digest:
            raise ValueError(f"Changed dual-feed retained evidence: {name}")
    limits=plan["numerical_acceptance"]
    maxima=record["numerical_maxima"]
    assert maxima["grid_change_absolute_conversion"]<=limits["max_absolute_grid_difference"]
    assert maxima["normalized_balance_error"]<=limits["max_normalized_reaction_flux_balance_error"]
    assert maxima["film_balance_error_mol_L_min"]<=limits["max_film_root_error_mol_L_min"]
    assert maxima["max_positive_rate_slope_over_film_capacity"]<limits["max_sampled_positive_slope_over_film_capacity"]
    assert record["catalyst_mass_kg"]==.006 and record["condition_count"]==40
    assert len(record["central_scores"])==22
    assert all(r["joint_40"]["count"]==40 for r in record["central_scores"].values())
    assert record["prior_2to1_replay"]["max_central_absolute_difference"]<=limits["prior_2to1_replay_tolerance_absolute_conversion"]
    assert record["prior_2to1_replay"]["max_clock_endpoint_absolute_difference"]<=limits["prior_2to1_replay_tolerance_absolute_conversion"]
    for name,count,keys in (("bed_comparison.csv",880,("case","feed_ratio")),
        ("clock_comparison.csv",560,("law","feed_ratio","catalyst_age_min")),
        ("transport_comparison.csv",800,("case","feed_ratio","transport_case"))):
        rows=read_rows(REFERENCE/name)
        assert len(rows)==count and all(float(r["actual_mass_kg"])==.006 for r in rows)
        assert set(Counter(tuple(r[k] for k in keys) for r in rows).values())=={20}
    replay=read_json(REFERENCE/"replay.json")
    assert replay["all_ledgers_byte_identical"] and replay["output_sha256"]==record["output_sha256"]
    sign=read_json(REFERENCE/"sign_summary.json")
    for name,digest in sign["source_sha256"].items():
        if sha(ROOT/name)!=digest:
            raise ValueError(f"Changed Table S3 sign audit input: {name}")
    assert sign["catalyst_mass_kg"]==.006 and sign["batch_count"]==450 and sign["bed_count"]==80
    assert sign["max_batch_step_refinement"]<1e-5
    assert sign["max_bed_grid_change"]<=limits["max_absolute_grid_difference"]
    assert sign["max_balance_error"]<=limits["max_normalized_reaction_flux_balance_error"]
    assert sign["max_film_root_error_mol_L_min"]<=limits["max_film_root_error_mol_L_min"]
    assert sign["previous_bed_comparison_sha256"]==record["output_sha256"]["bed_comparison.csv"]
    assert replay["sign_ledgers_byte_identical"] and replay["sign_output_sha256"]==sign["output_sha256"]


if __name__=="__main__":
    check()
    print("Dual-feed comparison: all frozen inputs, numerical limits and independent replay hashes verified.")
