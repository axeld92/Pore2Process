"""Verify stored artifacts and recompute scale translation with the real Julia core."""
from pathlib import Path
import argparse
import json
import math

from coupling.run import sha, validate_process
from hydraulic.run import ROOT, analyze, write_json


def verify(path):
    path=Path(path).resolve()
    summary=json.loads(path.read_text())
    if summary["status"]!="passed" or not all(summary["checks"].values()):
        raise RuntimeError("The selected run did not pass acceptance.")
    for relative,digest in summary["artifacts_sha256"].items():
        file=path.parent / relative
        if sha(file)!=digest:
            raise RuntimeError(f"Recorded artifact changed: {relative}")
    from pore2process_demo.bridge import runtime
    jl=runtime()
    for design in summary["designs"]:
        design_path=path.parent / design["name"]
        if design["dwsim"]!=json.loads((design_path / "dwsim/summary.json").read_text()):
            raise RuntimeError("DWSIM summary differs from the process snapshot.")
        if design["dwsim"]["inputs"]!=json.loads((design_path / "dwsim/inputs.json").read_text()):
            raise RuntimeError("DWSIM input snapshot differs from its actual input file.")
        for case in design["pore_cases"]:
            directory=path.parent / design["name"] / f"n{case['resolution']}"
            raw=json.loads((directory / "raw.json").read_text())
            if raw!=case["inputs_and_measurements"]:
                raise RuntimeError("Raw measurement and input snapshots disagree.")
            metrics=analyze(directory,raw,jl)
            for name,value in metrics.items():
                if not math.isclose(value,case["metrics"][name],rel_tol=1e-12,abs_tol=1e-30):
                    raise RuntimeError(f"Recomputed Julia metric differs: {design['name']}/{name}")
        t=design["translation"]
        if t!=json.loads((design_path / "translation.json").read_text()):
            raise RuntimeError("Translation snapshot differs from its recorded file.")
        k=design["pore_cases"][-1]["metrics"]["permeability_cell_m2"]
        settings=summary["process_inputs"]
        dp=float(jl.predict_column_pressure_drop(k,t["viscosity_Pa_s"],settings["column_length_m"],
                                                settings["column_area_m2"],settings["volume_flow_m3_s"]))
        if not math.isclose(dp,t["pressure_drop_Pa"],rel_tol=1e-12):
            raise RuntimeError("Measured permeability did not reproduce the recorded pressure load.")
        if not math.isclose(design["dwsim"]["inputs"]["Outlet_pressure_Pa"],settings["inlet_absolute_pressure_Pa"]+dp,rel_tol=1e-12):
            raise RuntimeError("DWSIM absolute pressure does not contain the Julia load.")
        validate_process(design["dwsim"],design["dwsim"]["inputs"],t["density_kg_m3"],t["viscosity_Pa_s"],dp)
    result={"status":"passed","run_id":summary["run_id"],
            "artifact_count":len(summary["artifacts_sha256"]),
            "julia_translation_recomputed":True,"stored_DWSIM_balances_rechecked":True,
            "scope":"Artifact/translation verification; DWSIM is not re-solved by this checker. Each original adapter run saved, reloaded and solved the actual flowsheet."}
    write_json(path.parent / "verification.json",result)
    print(f"Verified {result['artifact_count']} saved artifacts and both Julia → DWSIM input contracts.")
    return result


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("summary",type=Path,nargs="?")
    selected=parser.parse_args().summary
    if selected is None:
        selected=Path(json.loads((ROOT / "results/coupling/latest.json").read_text())["summary"])
    verify(selected)
