"""Resolve the shared-contact 0.80 mm midpoint on the first 1,200 s."""
import argparse
from pathlib import Path

import numpy as np

from examples.adsorption.contact_column.run import coordinates
from examples.adsorption.transport_accuracy.run import validate_record
from pore2process_demo.io import read_json, sha, write_json

from .run import HERE, OUTPUT, conservation, crossings, solve

END_S = 1200.
LEVELS = ((28,.5,.006),(42,1/3,.004))


def summarize(results,feed):
    times = [crossings(r,feed)["t50"] for r in results]
    if any(t is None for t in times):
        raise ValueError("Short window does not cover the predicted midpoint")
    t,c = coordinates(results[0],feed)
    t,c = t[t<=END_S],c[t<=END_S]
    changes = [float(max(abs(c-np.interp(t,*coordinates(r,feed))))) for r in results[1:]]
    error = abs(times[-2]/times[-1]-1)
    return dict(t50_s=times,finest_pair_t50_change=error,passed_at_128_shells=error<=.02,
                max_fraction_change_to_full=changes)


def run(path):
    validate_record(path)
    base = path.parent
    frozen = read_json(base/"frozen_predictions.json")
    full_file = base/frozen["cases"]["d800"]["wall"]["full"][-1]
    full = read_json(full_file)
    model = read_json(base/"model.json")
    point = model["sizes"]["d800"]
    closure = model["closure"]
    output = OUTPUT/"plateau"/base.name
    output.mkdir(parents=True)
    results = [full["result"]]
    files = {}
    for i,(sub,dt,early) in enumerate(LEVELS):
        p = {**full["inputs"],"subcells":sub,"dt_s":dt,"startup_dt_s":early,"t_end_s":END_S}
        r = solve(closure,point["mix"],point["weights"]["wall"],p)
        conservation(r,p,point["trace"],point["weights"]["wall"])
        results.append(r)
        name = f"refinement{i}.json"
        write_json(output/name,dict(inputs=p,result=r))
        files[name] = sha(output/name)
        print(i,crossings(r,p["inlet_concentration_kg_m3"]),r["total_s"],flush=True)
    check = summarize(results,full["inputs"]["inlet_concentration_kg_m3"])
    write_json(output/"summary.json",dict(parent_sha256=sha(path),source_sha256=sha(HERE/"plateau.py"),
                                          files=files,analysis=check))
    print(check,flush=True)


def check(path):
    validate_record(path)
    base = path.parent
    saved_path = OUTPUT/"plateau"/base.name
    saved = read_json(saved_path/"summary.json")
    if sha(path)!=saved["parent_sha256"] or sha(HERE/"plateau.py")!=saved["source_sha256"]:
        raise ValueError("Plateau source or parent changed")
    frozen = read_json(base/"frozen_predictions.json")
    parent = read_json(base/frozen["cases"]["d800"]["wall"]["full"][-1])
    results = [parent["result"]]
    for file,digest in saved["files"].items():
        if sha(saved_path/file)!=digest:
            raise ValueError("Plateau history changed")
        results.append(read_json(saved_path/file)["result"])
    comparison = summarize(results,parent["inputs"]["inlet_concentration_kg_m3"])
    if comparison!=saved["analysis"]:
        raise ValueError("Plateau comparison differs on replay")
    print(comparison)


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("summary",type=Path)
    parser.add_argument("--check",action="store_true")
    args = parser.parse_args()
    check(args.summary.resolve()) if args.check else run(args.summary.resolve())
