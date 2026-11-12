"""Resolve the numerically sensitive 0.80 mm midpoint before source scoring."""
import argparse
from pathlib import Path

import numpy as np

from examples.adsorption.exchange_attribution import kernel
from examples.adsorption.exchange_attribution.run import allocation_check
from pore2process_demo.io import read_json, sha, write_json

from .run import HERE, OUTPUT, coordinates, crossings

LIMIT_S = 600.
LEVELS = ((12,128,.75,.01),(18,128,.5,.0075))


def comparison(results, feed):
    reference = results[0]
    t,c = coordinates(reference,feed)
    window = t<=LIMIT_S
    times = [crossings(case,feed)["t50"] for case in results]
    if any(x is None for x in times):
        raise ValueError("Plateau midpoint was not covered in the short run")
    profiles = [float(max(abs(c[window]-np.interp(t[window],*coordinates(case,feed))))) for case in results[1:]]
    return dict(t50_s=times,maximum_fraction_changes_to_full=profiles,
                finest_pair_t50_relative=abs(times[-2]/times[-1]-1),
                qualified_at_128_shells=abs(times[-2]/times[-1]-1)<=.02)


def run(base):
    model_file = base/"model.json"
    archived_file = base/"d800_advection_rank_level2.json"
    model = read_json(model_file)
    archived = read_json(archived_file)
    points = model["sizes"]["d800"]
    result = OUTPUT/"plateau"/base.name
    result.mkdir(parents=True)
    generated = result/"AllocationColumn.jl"
    generated.write_text(kernel.source())
    params = archived["inputs"]
    closure = model["closure"]
    allocation = dict(closure=closure,allocations=points["allocations"])
    shares = allocation["allocations"]["advection_rank"]
    histories = [archived["result"]]
    print("Refining 0.80 mm midpoint through 600 s",result,flush=True)
    files = {}
    for n,(subcells,radial,dt,startup) in enumerate(LEVELS):
        p = {**params,"subcells":subcells,"radial_cells":radial,"dt_s":dt,"startup_dt_s":startup,"t_end_s":LIMIT_S}
        r = kernel.solve(generated,closure["flow_fractions"],closure["normalized_transit_times"],points["mix"],shares,**p)
        allocation_check(r,p,archived["trace"],allocation,"advection_rank")
        history = result/f"plateau{n}.json"
        write_json(history,dict(inputs=p,result=r))
        files[history.name] = sha(history)
        histories.append(r)
        print(n,crossings(r,p["inlet_concentration_kg_m3"]),r["total_s"],flush=True)
    report = comparison(histories,params["inlet_concentration_kg_m3"])
    write_json(result/"summary.json",dict(parent_model_sha256=sha(model_file),parent_full_sha256=sha(archived_file),
        study_source_sha256=sha(HERE/"plateau.py"),generated_kernel_sha256=sha(generated),files=files,comparison=report))
    print(report,flush=True)


def check(base):
    result = OUTPUT/"plateau"/base.name
    record = read_json(result/"summary.json")
    parent = base/"d800_advection_rank_level2.json"
    if (sha(base/"model.json")!=record["parent_model_sha256"] or sha(parent)!=record["parent_full_sha256"]
            or sha(HERE/"plateau.py")!=record["study_source_sha256"]
            or sha(result/"AllocationColumn.jl")!=record["generated_kernel_sha256"]
            or (result/"AllocationColumn.jl").read_text()!=kernel.source()):
        raise ValueError("Plateau inputs changed")
    archived = read_json(parent)
    histories = [archived["result"]]
    for filename,digest in record["files"].items():
        if sha(result/filename)!=digest:
            raise ValueError("Plateau history changed")
        histories.append(read_json(result/filename)["result"])
    if comparison(histories,archived["inputs"]["inlet_concentration_kg_m3"])!=record["comparison"]:
        raise ValueError("Plateau comparison differs on replay")
    print(record["comparison"])


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cohort",type=Path)
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    check(args.cohort.resolve()) if args.check else run(args.cohort.resolve())
