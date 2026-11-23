"""Check compact verification and exact retained input hashes."""
from pathlib import Path
from pore2process_demo.io import read_json,sha
from pore2process_demo.paths import ROOT

REFERENCE=Path(__file__).resolve().parent/"reference"


def check():
    record=read_json(REFERENCE/"summary.json")
    for name,digest in record["source_sha256"].items():
        if sha(ROOT/name)!=digest:
            raise ValueError(f"Changed composition-audit source: {name}")
    for name,digest in record["retained_sha256"].items():
        if sha(REFERENCE/name)!=digest:
            raise ValueError(f"Changed retained composition-audit result: {name}")
    if record["bed_numerics"]["actual_mass_kg"]!=.006 or any(
        v["count"]!=20 for v in record["bed_comparison"].values()):
        raise ValueError("The measured catalyst inventory and all twenty reactor observations must remain unchanged.")
    if (REFERENCE/"replay.json").exists():
        replay=read_json(REFERENCE/"replay.json")
        if replay["primary_cases_repeated"]!=10 or replay["models_per_case"]!=13 or any(
            v>1e-9 for k,v in replay.items() if k.startswith("max_")):
            raise ValueError("Invalid primary reproducibility record.")


if __name__=="__main__":
    check()
    print("Composition audit: immutable input and retained-result hashes verified.")
