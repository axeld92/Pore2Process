"""Check the compact conditioning records against their frozen inputs."""

from pathlib import Path

from pore2process_demo.io import read_json, sha
from pore2process_demo.paths import ROOT


FILES = (
    "batch_activity.csv", "batch_storage.csv", "hydration_summary.csv",
    "predictions_before_comparison.csv", "bed_comparison.csv",
)


def check():
    reference = Path(__file__).resolve().parent/"reference"
    summary = read_json(reference/"summary.json")
    for name, digest in summary["source_sha256"].items():
        if sha(ROOT/name) != digest:
            raise ValueError(f"A conditioning source or model changed: {name}")
    for name in FILES:
        if sha(reference/name) != summary["output_sha256"][name]:
            raise ValueError(f"A retained conditioning result changed: {name}")


if __name__ == "__main__":
    check()
    print("The retained conditioning outputs match their recorded input and output hashes.")
