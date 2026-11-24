"""Cold process repeat of the full refinement decision on prespecified histories."""
from argparse import ArgumentParser
from pathlib import Path
from time import perf_counter

import numpy as np

from examples.biocatalysis.pentyl_propanoate_history import aggregate,run
from pore2process_demo.io import read_json,sha,write_json
from pore2process_demo.paths import ROOT


def selected(c):
    # Repeat every law/scenario at a hot and fast condition on both feeds.
    # Include severe slow-bed exposure and a moderate representative condition.
    if c["temperature_C"]==80 and c["flow_ml_min"]==3:
        return True
    return (c["scenario"]=="fresh_bulk" and c["law"] in ("acid_both","protected_residual_decay") and
            (c["temperature_C"],c["flow_ml_min"]) in ((80,.7),(60,3)))


def replay(combined,output):
    record=read_json(combined/"summary.json")
    for name,digest in record["output_sha256"].items():
        if sha(combined/name)!=digest:
            raise ValueError(f"Combined prediction ledger changed: {name}")
    old={r["task_id"]:aggregate.typed(r) for r in run.read_rows(combined/"numerics.csv")}
    profiles=read_json(combined/"profile_index.json")
    plan=read_json(run.PLAN); jl=run.load_runtime(); tasks,laws=run.prepare(jl,plan)
    chosen=[c for c in tasks if selected(c)]
    if len(chosen)!=36:
        raise ValueError("Unexpected cold-replay boundary.")
    diffs={k:0. for k in run.ARRAYS}; scalar_diff=0.; decisions=[]; start=perf_counter()
    for c in chosen:
        r,d=run.solve_refined(jl,laws[c["law"]],c,plan)
        prior=old[c["task_id"]]
        for k in ("fine_cells","fine_max_step_min","coarse_cells","coarse_max_step_min","numerical_solve_count","accepted_steps","rejected_steps"):
            if d[k]!=prior[k]:
                raise RuntimeError(f"Changed refinement or solver decision: {c['task_id']} {k}")
        p=profiles[c["task_id"]]; path=ROOT/p["path"]
        if sha(path)!=p["sha256"]:
            raise RuntimeError("Saved profile changed before replay.")
        with np.load(path,allow_pickle=False) as oldp:
            for k in run.ARRAYS:
                v=float(np.max(np.abs(r[k]-oldp[k])))
                diffs[k]=max(diffs[k],v)
        scalar_diff=max(scalar_diff,max(abs(float(v)-float(prior[k])) for k,v in d.items()))
        decisions.append({"task_id":c["task_id"],"fine_cells":d["fine_cells"],"max_step_min":d["fine_max_step_min"],
                          "numerical_solve_count":d["numerical_solve_count"]})
        print(f"Replayed {len(decisions)}/{len(chosen)} histories ({perf_counter()-start:.1f} s)",flush=True)
    if max(diffs.values())>1e-11 or scalar_diff>1e-11:
        raise RuntimeError(f"Transient replay changed: {diffs}, scalar={scalar_diff}")
    result={"repeated_histories":len(chosen),"all_refinement_decisions_identical":True,
            "max_absolute_difference_by_array":diffs,"max_absolute_numerical_scalar_difference":scalar_diff,
            "exact_array_replay":max(diffs.values())==0.,"replayed_decisions":decisions,
            "combined_summary_sha256":sha(combined/"summary.json"),"replay_source_sha256":sha(Path(__file__))}
    if output.exists():
        raise FileExistsError(output)
    output.parent.mkdir(parents=True,exist_ok=True)
    write_json(output,result)
    return result


if __name__=="__main__":
    p=ArgumentParser(description=__doc__)
    p.add_argument("--combined",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args(); replay(args.combined.resolve(),args.output.resolve())
