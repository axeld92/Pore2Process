"""Verify the compact loading audit against its frozen sources and run outputs."""

from pathlib import Path

from pore2process_demo.io import read_json, sha
from pore2process_demo.paths import ROOT


REFERENCE = Path(__file__).resolve().parent/"reference"
FILES = (
    "matched_conversion.csv", "shape_free_holdout.csv", "batch_fits.csv",
    "batch_fold_scores.csv", "predictions_before_comparison.csv", "bed_comparison.csv",
)


def check():
    summary = read_json(REFERENCE/"summary.json")
    for name,digest in summary["source_sha256"].items():
        if sha(ROOT/name) != digest:
            raise ValueError(f"Loading-audit input or model changed: {name}")
    for name in FILES:
        if sha(REFERENCE/name) != summary["output_sha256"][name]:
            raise ValueError(f"Retained loading-audit output changed: {name}")
    scope = summary["bed_transfer"]
    if scope["catalyst_mass_kg"] != .006:
        raise ValueError("The solid inventory must remain 6 g.")


if __name__ == "__main__":
    check()
    print("The retained loading audit matches its recorded sources and output hashes.")
