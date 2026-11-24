"""Retain compact pointwise history and numerical proof without bulk field files."""
from argparse import ArgumentParser
from pathlib import Path
from shutil import copyfile

from examples.biocatalysis.pentyl_propanoate_history import run
from pore2process_demo.io import read_json,sha,write_json


HISTORY_FIELDS=("task_id","scenario","law","feed_ratio","temperature_C","flow_ml_min","time_min",
                "predicted_conversion","mean_liquid_progress","mean_activity","minimum_activity","maximum_activity",
                "first_quarter_activity","last_quarter_activity","mean_cumulative_hazard",
                "generated_mol","effluent_product_mol","backdiffused_product_mol","liquid_product_mol")


def retain(combined,replay,output):
    record=read_json(combined/"summary.json"); repeat=read_json(replay)
    for name,digest in record["output_sha256"].items():
        if sha(combined/name)!=digest:
            raise ValueError(f"Changed full output ledger: {name}")
    if repeat["combined_summary_sha256"]!=sha(combined/"summary.json"):
        raise ValueError("Replay does not refer to the same joined run.")
    output.mkdir(parents=True,exist_ok=False)
    for name in ("bed_comparison.csv","numerics.csv","uniform_envelope_check.csv","selected_axial_profiles.csv","reaction_zone.csv"):
        copyfile(combined/name,output/name)
    run.write_rows(output/"history.csv",[{k:r[k] for k in HISTORY_FIELDS} for r in run.read_rows(combined/"predictions_before_comparison.csv")])
    copyfile(replay,output/"replay.json")
    compact=record|{"full_summary_sha256":sha(combined/"summary.json"),
                    "retained_sha256":{p.name:sha(p) for p in sorted(output.iterdir())}}
    write_json(output/"summary.json",compact)


if __name__=="__main__":
    p=ArgumentParser(description=__doc__)
    p.add_argument("--combined",type=Path,required=True)
    p.add_argument("--replay",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args(); retain(args.combined.resolve(),args.replay.resolve(),args.output.resolve())
