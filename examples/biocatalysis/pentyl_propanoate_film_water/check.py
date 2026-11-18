"""Check the retained film and batch comparisons against their input hashes."""

from pathlib import Path

from pore2process_demo.io import read_json, sha
from pore2process_demo.paths import ROOT


HERE = Path(__file__).resolve().parent
RECORDS = (
    ("film", "bed_comparison.csv", "bed_comparison_sha256"),
    ("water_activity", "comparison.csv", "comparison_sha256"),
)


def check():
    for stage, table, digest_key in RECORDS:
        reference = HERE / stage / "reference"
        record = read_json(reference / "summary.json")
        for path, expected in record["source_sha256"].items():
            if sha(ROOT / path) != expected:
                raise ValueError(f"The {stage} comparison source changed: {path}")
        if sha(reference / table) != record[digest_key]:
            raise ValueError(f"The {stage} comparison table changed.")


if __name__ == "__main__":
    check()
    print("Both comparisons match the retained input and output hashes.")
