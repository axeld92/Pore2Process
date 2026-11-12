"""Run flux-weighted wall tracers on a frozen, matched-Péclet field."""
import argparse
from pathlib import Path

import numpy as np

from examples.adsorption.dispersion_holdout.run import execute
from examples.adsorption.repeated_geometry.diagnostics import path_persistence
from examples.adsorption.repeated_geometry.predict import observe
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.provenance import source_manifest

from .mesh import COMPACT, choose_cut, compact
from .native import HERE, ORDERED
from .packing import PERIODS, SIDE_M

COUNT = 8192
SEED = 3722747
STOP = 12.
THREADS = 5


def tracer_plan(base,label,n,cut):
    directory = base/f"{label}_N{n}"
    info = read_json(directory/"faces/summary.json")["info"]
    config = read_json(base/"plan.json")
    ordered = read_json(ORDERED/"frozen_prediction.json")
    julia, _, *flags = ordered["parity_command"][:4]
    runtime = [julia,f"--threads={THREADS}",*flags]
    carrier_u = info["interstitial_velocity_m_s"][0]
    ratio = carrier_u/config["reference_interstitial_m_s"]
    cases = {}
    d = directory/"compact"
    settings = read_json(d/"summary.json")
    if settings["cut_index"] != cut:
        raise ValueError("Compact field cut changed")
    ends = ",".join(map(str,settings["period_end_counts"]))
    for name,old in ordered["predictions"].items():
        if n==48 and name.endswith("_half_step"):
            continue
        dt = old["dt_s"]
        diffusion = old["diffusivity_m2_s"]*ratio
        pe = carrier_u*SIDE_M/diffusion
        target_pe = config["reference_interstitial_m_s"]*SIDE_M/old["diffusivity_m2_s"]
        np.testing.assert_allclose(pe,target_pe,rtol=1e-14,atol=0.)
        output = directory/(name+".f64")
        cmd = [*runtime,str(HERE/"track.jl"),str(d),str(output),
               *map(str,(info["raw"]["dx_m"],diffusion,dt,STOP,COUNT,SEED,n,ends))]
        cases[name] = dict(command=cmd,diffusivity_m2_s=diffusion,dt_s=dt,
                           Pe=pe,target_Pe=target_pe)
    arrays = [name+(".i32" if dt=="<i4" else ".f64") for name,dt in COMPACT.items()]
    return dict(label=label,grid=n,compact=settings,compact_sha256={name:sha(d/name) for name in arrays},
                cases=cases,marker_count=COUNT,seed=SEED,stop_s=STOP,
                source_sha256=source_manifest([HERE/"track.jl",HERE/"tracer.py",
                    HERE.parent/"repeated_geometry/track.jl",HERE/"protocol.md"]))


def observe_cases(directory,plan):
    result = {}
    for name,request in plan["cases"].items():
        p = directory/(name+".f64")
        data = np.fromfile(p,dtype="<f8").reshape(plan["marker_count"],2*PERIODS)
        tracked = observe(data)
        recovered = tracked["recovery"][-1] >= .999 and tracked["maximum_unresolved_fraction"] <= .001
        result[name] = dict(observed=tracked,inert_recovery_passed=recovered,persistence=path_persistence(data),
                            Pe=request["Pe"],diffusivity_m2_s=request["diffusivity_m2_s"],
                            archive_sha256=sha(p))
    return result


def run(base,label,n,requested_cut):
    directory = base/f"{label}_N{n}"
    if (directory/"tracer_summary.json").exists():
        raise ValueError("Completed tracer archive already exists")
    if requested_cut is None:
        cut = 0 if label=="ordered" else choose_cut(base)*(n//16)
    else:
        cut = requested_cut*(n//16)
    compact(directory,cut)
    plan = tracer_plan(base,label,n,cut)
    path = directory/"tracer_plan.json"
    write_json(path,plan)
    print("Frozen matched tracer",directory,sha(path),flush=True)
    for name,case in plan["cases"].items():
        elapsed = execute(case["command"],directory/(name+".log"),timeout=3600)
        write_json(directory/(name+"_execution.json"),elapsed)
        print(label,n,name,elapsed,flush=True)
    observed = observe_cases(directory,plan)
    write_json(directory/"tracer_summary.json",dict(tracer_plan_sha256=sha(path),cases=observed))
    print("Matched tracer complete",directory,flush=True)


if __name__=="__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("base",type=Path)
    p.add_argument("label",choices=("ordered","disordered"))
    p.add_argument("grid",type=int)
    p.add_argument("--cut",type=int,help="Axial inlet in units of one sixteenth of the reference length")
    a = p.parse_args()
    run(a.base.resolve(),a.label,a.grid,a.cut)
