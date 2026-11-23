"""Retain the compact, pointwise dual-feed audit and exact run/source hashes."""
from argparse import ArgumentParser
from pathlib import Path
from shutil import copyfile

from examples.biocatalysis.pentyl_propanoate.run import read_rows, write_rows
from pore2process_demo.io import read_json,sha,write_json

KEYS = ("feed_ratio","temperature_C","flow_ml_min","source_run_id","actual_mass_kg")
COLUMNS = {
    "bed_comparison.csv": ("case",*KEYS,"acid0_mol_L","eta_internal","activity","apparent_Kx",
        "equilibrium_conversion","observed_conversion","predicted_conversion","conversion_error_pp",
        "nominal_outlet_ester_g_h","observed_nominal_outlet_ester_g_h","inlet_backdiffusion_fraction",
        "grid_change_absolute_conversion","normalized_balance_error","film_balance_error_mol_L_min"),
    "clock_comparison.csv": ("law",*KEYS,"catalyst_age_min","inlet_activity","equilibrium_activity",
        "inlet_hazard_min_inv","equilibrium_hazard_min_inv","conversion_uniform_inlet_activity",
        "conversion_uniform_equilibrium_activity","observed_conversion","grid_change_absolute_conversion"),
    "transport_comparison.csv": ("case",*KEYS,"transport_case","diffusivity_used_m2_s","eta_internal",
        "observed_conversion","predicted_conversion","change_from_central_pp",
        "grid_change_absolute_conversion","normalized_balance_error"),
}


def verified_summary(path, name):
    record=read_json(path/name)
    for ledger,digest in record["output_sha256"].items():
        if sha(path/ledger) != digest:
            raise ValueError(f"Run ledger changed: {ledger}")
    return record


def retain(run,replay_run,sign_run,sign_replay,output):
    record=verified_summary(run,"summary.json")
    replay=verified_summary(replay_run,"summary.json")
    sign=verified_summary(sign_run,"sign_summary.json")
    sign_repeat=verified_summary(sign_replay,"sign_summary.json")
    if record!=replay or sign!=sign_repeat:
        raise ValueError("A cold independent repeat changed the numerical summaries or ledgers.")
    output.mkdir(parents=True,exist_ok=False)
    for filename,columns in COLUMNS.items():
        write_rows(output/filename,[{k:r[k] for k in columns} for r in read_rows(run/filename)])
    for filename in ("sign_bed_comparison.csv","sign_batch_comparison.csv","sign_summary.json"):
        copyfile(sign_run/filename,output/filename)
    write_json(output/"replay.json",{"all_ledgers_byte_identical":True,
        "sign_ledgers_byte_identical":True,"central_predictions_repeated":record["central_count"],
        "clock_intervals_repeated":record["clock_count"],"transport_predictions_repeated":record["transport_count"],
        "sign_batch_points_repeated":sign["batch_count"],"sign_bed_points_repeated":sign["bed_count"],
        "output_sha256":record["output_sha256"],"sign_output_sha256":sign["output_sha256"]})
    compact=record | {"full_summary_sha256":sha(run/"summary.json"),
                      "retained_sha256":{p.name:sha(p) for p in sorted(output.iterdir())}}
    write_json(output/"summary.json",compact)


if __name__=="__main__":
    parser=ArgumentParser(description=__doc__)
    parser.add_argument("--run",type=Path,required=True)
    parser.add_argument("--replay-run",type=Path,required=True)
    parser.add_argument("--sign-run",type=Path,required=True)
    parser.add_argument("--sign-replay",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    retain(args.run.resolve(),args.replay_run.resolve(),args.sign_run.resolve(),args.sign_replay.resolve(),args.output.resolve())
