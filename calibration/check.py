"""Replay saved inputs through Julia and verify numerical/export integrity."""
import argparse
import csv
import math
from pathlib import Path

from calibration.run import analyze, read_json, sha256


def close_tree(actual, expected, path="result"):
    if isinstance(expected,dict):
        if not isinstance(actual,dict) or actual.keys()!=expected.keys():
            raise ValueError(f"Different fields at {path}")
        for key in expected:
            close_tree(actual[key],expected[key],f"{path}.{key}")
    elif isinstance(expected,list):
        if not isinstance(actual,list) or len(actual)!=len(expected):
            raise ValueError(f"Different array length at {path}")
        for i,(a,e) in enumerate(zip(actual,expected,strict=True)):
            close_tree(a,e,f"{path}[{i}]")
    elif isinstance(expected,(float,int)) and not isinstance(expected,bool):
        if isinstance(actual,bool) or not isinstance(actual,(float,int)) or not math.isclose(actual,expected,rel_tol=1e-11,abs_tol=0.):
            raise ValueError(f"Different numerical value at {path}: {actual!r} vs {expected!r}")
    elif actual != expected:
        raise ValueError(f"Different value at {path}")


def check(directory,require_screen_pass=False):
    directory=Path(directory)
    summary=read_json(directory/"summary.json")
    required={"inputs.json","predictions.csv","permeability_validation.png","conditional_energy.png"}
    if set(summary["artifact_sha256"]) != required:
        raise ValueError("Missing or unexpected artifact declarations")
    for name,digest in summary["artifact_sha256"].items():
        if sha256(directory/name)!=digest:
            raise ValueError(f"Artifact hash mismatch: {name}")
    snapshot=read_json(directory/"inputs.json")
    replay=analyze(snapshot)
    for field in ("fit","cases","metrics","leave_one_calibration_out","primary_screen","conditional_process"):
        close_tree(summary[field],replay[field],field)
    if summary["protocol_id"] != snapshot["protocol"]["protocol_id"]:
        raise ValueError("Summary protocol differs from saved input")
    with (directory/"predictions.csv").open(newline="") as stream:
        records=list(csv.DictReader(stream))
    if len(records)!=len(replay["cases"]):
        raise ValueError("CSV row count differs")
    for record,expected in zip(records,replay["cases"],strict=True):
        for key,value in expected.items():
            actual=None if value is None and record[key]=="" else float(record[key]) if isinstance(value,(int,float)) else record[key]
            close_tree(actual,value,f"CSV.{key}")
    if require_screen_pass and not replay["primary_screen"]["passed"]:
        raise ValueError(f"Execution/artifacts verified, but scientific screen failed: {replay['primary_screen']}")
    return {"integrity_and_numerical_replay":"passed","cases":len(replay["cases"]),"primary_screen":replay["primary_screen"]}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory",type=Path)
    parser.add_argument("--require-screen-pass",action="store_true")
    args=parser.parse_args()
    print(check(args.directory,args.require_screen_pass))


if __name__=="__main__":
    main()
