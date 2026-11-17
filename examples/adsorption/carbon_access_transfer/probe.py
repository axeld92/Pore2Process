"""Probe the saved TRT velocity shape at the prospective flow Péclet numbers."""
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from uuid import uuid4

import numpy as np

from examples.adsorption.disordered_exchange.run import PARENT_SHA, parent_inputs
from examples.adsorption.disordered_section_exchange.run import errors, verify_saved_field
from examples.adsorption.disordered_travel_memory.run import HERE as MEMORY
from examples.adsorption.dispersion_holdout.run import execute
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT
from pore2process_demo.provenance import source_manifest

HERE = Path(__file__).resolve().parent
OUTPUT = ROOT/"results/adsorption/carbon_access_transfer"
REFERENCE_FLOW = 43.8
NEW_FLOWS = {"q233":23.3,"q625":62.5}


def settings(base,plan):
    ref = plan["cases"]["d800_primary"]
    settings = {}
    n = plan["marker_count"]
    for name,flow in NEW_FLOWS.items():
        ratio = flow/REFERENCE_FLOW
        command = ref["command"].copy()
        command[6] = str(base/name/"markers.f64")
        command[8] = str(ref["diffusivity_m2_s"]/ratio)
        command[11:13] = [str(3*n),str(plan["seed"]+8*n)]
        settings[name] = dict(flow_mL_min=flow,velocity_ratio=ratio,target_Pe=ref["Pe"]*ratio,
            archive_coordinate_diffusivity_m2_s=float(command[8]),archive_coordinate_dt_s=ref["dt_s"],
            training_count=n,operator_probe_count=2*n,command=command)
    return settings


def trial(base,name,spec):
    case = base/name
    case.mkdir()
    n = spec["training_count"]
    timing = execute(spec["command"],case/"tracker.log",timeout=2400)
    markers = np.fromfile(case/"markers.f64",dtype="<f8").reshape(3*n,16)
    # The archived group supplies calibration sections 1–2 only.
    markers[:n,:2].copy(order="C").tofile(case/"training.f64")
    markers[n:,:].copy(order="C").tofile(case/"operator.f64")
    markers[n:,:8].copy(order="C").tofile(case/"passages.f64")
    solve = [*spec["command"][:4],str(MEMORY/"solve.jl"),
             str(case/"training.f64"),str(case/"passages.f64"),str(case)]
    compute = execute(solve,case/"solve.log",timeout=120)
    frozen = {item:sha(case/item) for item in ("training.f64","operator.f64","rules.f64","four.csv")}
    write_json(case/"frozen.json",dict(hashes=frozen,tracker=timing,coarsen=compute))
    arrivals = markers[:n,:8]
    if not np.isfinite(arrivals).all() or np.any(np.diff(arrivals,axis=1)<=0):
        raise ValueError(f"Incomplete independent arrivals: {name}")
    observed = np.column_stack((arrivals[:,2:].mean(axis=0),arrivals[:,2:].var(axis=0)))
    scored = errors(observed,np.loadtxt(case/"four.csv",delimiter=","))
    write_json(case/"holdout.json",dict(frozen_sha256=sha(case/"frozen.json"),score=scored))
    print(name,"tracer variance max error",scored["maximum_variance_error"],flush=True)
    return scored


def run():
    parent,plan = parent_inputs()
    verify_saved_field(parent,plan)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")+"-"+uuid4().hex[:8]
    base = OUTPUT/"probes"/stamp
    base.mkdir(parents=True)
    specs = settings(base,plan)
    sources = source_manifest([HERE/"probe.py",MEMORY/"solve.jl",MEMORY/"model.jl"])
    write_json(base/"plan.json",dict(parent_sha256=PARENT_SHA,cases=specs))
    print("Probes started on saved TRT field",base,flush=True)
    with ThreadPoolExecutor(max_workers=2) as pool:
        jobs = {name:pool.submit(trial,base,name,spec) for name,spec in specs.items()}
        results = {name:job.result() for name,job in jobs.items()}
    write_json(base/"summary.json",dict(parent_sha256=PARENT_SHA,source_sha256=sources,
        artifacts_sha256={p.relative_to(base).as_posix():sha(p) for p in sorted(base.rglob("*")) if p.is_file()},
        cases=results))
    write_json(OUTPUT/"latest_probes.json",dict(summary=(base/"summary.json").relative_to(ROOT).as_posix()))
    print(base/"summary.json",flush=True)


if __name__=="__main__":
    run()
