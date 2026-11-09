"""Recheck immutable OpenLB exports and rerun Julia translation without solving."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import numpy as np

from run import ROOT, analyze, write_json


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("summary",type=Path,nargs="?")
    selected=parser.parse_args().summary
    if selected is None:
        selected=Path(json.loads((ROOT / "results/hydraulic/latest.json").read_text())["summary"])
    summary_path=selected.resolve()
    summary=json.loads(summary_path.read_text())
    if summary["status"] != "passed" or not all(summary["checks"].values()):
        raise SystemExit("This snapshot did not pass the hydraulic acceptance suite.")
    from pore2process_demo.bridge import runtime
    jl=runtime()
    checked=[]
    for case in summary["cases"]:
        directory=summary_path.parent / case["name"]
        for name,digest in case["artifacts_sha256"].items():
            if hashlib.sha256((directory / name).read_bytes()).hexdigest() != digest:
                raise RuntimeError(f"Solver artifact changed after simulation: {directory / name}")
        raw=json.loads((directory / "raw.json").read_text())
        if raw != case["inputs_and_measurements"]:
            raise RuntimeError("Input/measurement snapshot differs from raw solver output.")
        metrics=analyze(directory,raw,jl)
        for name,value in metrics.items():
            if not np.isclose(value,case["metrics"][name],rtol=1e-12,atol=1e-30):
                raise RuntimeError(f"Recomputed metric differs: {case['name']} {name}")
        checked.append(case["name"])
    names=("hydraulic/run.py","hydraulic/check_exports.py","julia/src/Hydraulics.jl",
           "julia/Project.toml","julia/Manifest.toml","uv.lock")
    record={"timestamp_utc":datetime.now(timezone.utc).isoformat(),"run_id":summary["run_id"],
            "status":"passed","checked_cases":checked,"solver_artifact_hashes_verified":True,
            "si_coordinates_and_int64_indices_verified":True,"julia_translation_recomputed":True,
            "julia_version":str(jl.seval("VERSION")),
            "implementation_sha256":{name:hashlib.sha256((ROOT / name).read_bytes()).hexdigest() for name in names}}
    write_json(summary_path.parent / "analysis_verification.json",record)
    print(f"Verified hashes, SI coordinates/indexing, and Julia-derived metrics for {len(checked)} actual OpenLB cases.")


if __name__ == "__main__":
    main()
