"""Verify the retained activity comparison against its frozen inputs."""

from pathlib import Path

from pore2process_demo.io import read_json, sha
from pore2process_demo.paths import ROOT


def check():
    reference = Path(__file__).resolve().parent / "reference"
    summary = read_json(reference / "summary.json")
    for name, digest in summary["source_sha256"].items():
        if sha(ROOT / name) != digest:
            raise ValueError(f"An activity comparison input changed: {name}")
    for name in ("headspace.csv", "representative.csv"):
        if sha(reference / name) != summary["output_sha256"][name]:
            raise ValueError(f"The retained activity comparison changed: {name}")


if __name__ == "__main__":
    check()
    print("The retained comparison matches its input and output hashes.")
