"""Verify the retained new audit and its provenance without reopening private PDFs."""
from pathlib import Path

from pore2process_demo.io import read_json,sha
from pore2process_demo.paths import ROOT

REFERENCE=Path(__file__).resolve().parent/"reference"
FILES=("matched_progress.csv","external_fold_scores.csv","external_progress_predictions.csv","batch_fits.csv",
       "predictions_before_comparison.csv","bed_comparison.csv")


def check():
    summary=read_json(REFERENCE/"summary.json")
    for name,digest in summary["source_sha256"].items():
        if sha(ROOT/name)!=digest:
            raise ValueError(f"Changed independent loading evidence: {name}")
    for name in FILES:
        if sha(REFERENCE/name)!=summary["output_sha256"][name]:
            raise ValueError(f"Changed retained numerical result: {name}")
    if summary["transfer"]["catalyst_mass_kg"] != .006:
        raise ValueError("Keep the measured six gram inventory.")
    if summary["bed_comparison"]["film_mass_proportional"]["count"] != 20:
        raise ValueError("The reactor comparison must retain every observation.")


if __name__=="__main__":
    check()
    print("The independent loading audit matches its source and result hashes.")
