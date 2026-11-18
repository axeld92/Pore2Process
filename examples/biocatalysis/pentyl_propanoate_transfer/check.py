"""Check the curated transfer evidence without running the solvers."""

from pathlib import Path

from pore2process_demo.io import read_json, sha
from pore2process_demo.paths import ROOT


HERE = Path(__file__).resolve().parent
RECORDS = (
    ("mixing/reference/summary.json", "source_hashes",
     "mixing/reference/bed_comparison.csv", "bed_comparison_sha256"),
    ("published/reference_summary.json", "source_sha256",
     "published/reference_comparison.csv", "comparison_sha256"),
    ("reference/summary.json", "inputs_sha256",
     "reference/bed_comparison.csv", "bed_comparison_sha256"),
)


def check():
    for filename, sources_key, table, table_key in RECORDS:
        record = read_json(HERE / filename)
        for name, digest in record[sources_key].items():
            if sha(ROOT / name) != digest:
                raise ValueError(f"{filename} differs from its source {name}.")
        if sha(HERE / table) != record[table_key]:
            raise ValueError(f"{table} differs from its recorded bytes.")


if __name__ == "__main__":
    check()
    print("The three retained transfer comparisons match their source records.")
