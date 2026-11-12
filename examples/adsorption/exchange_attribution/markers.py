"""Match diffusion-to-advection ratios in the archived synthetic velocity field."""
from datetime import datetime, timezone
from uuid import uuid4

import numpy as np

from examples.adsorption.dispersion_holdout.run import execute
from examples.adsorption.exchange_column.closure import ARCHIVE, DT, MARKERS, closure, group
from examples.adsorption.exchange_column.run import source_inputs
from examples.adsorption.path_exchange.predict import summary
from examples.adsorption.path_exchange.run import COUNT, HERE as TRACKER, SEED, STOP_S
from examples.adsorption.transport_accuracy.run import validate_record
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT

HERE = ROOT / "examples/adsorption/exchange_attribution"
OUTPUT = ROOT / "results/adsorption/exchange_attribution"
COLUMN = ROOT / "results/adsorption/exchange_column/runs/20261112T080827Z-1b52e000"
COLUMN_SHA = "5d3fb79f18f2ed0e2499c8ec7ab1e9fea310ed1406ea7ad85540a706d9073b45"


def inputs():
    path = COLUMN / "summary.json"
    if sha(path)!=COLUMN_SHA:
        raise ValueError("Original column experiment changed")
    validate_record(path)
    validate_record(MARKERS/"summary.json")
    cases = read_json(COLUMN/"frozen_predictions.json")["cases"]
    baselines = {k: {law: COLUMN/v[role] for law,role in (("synthetic","primary"),("none","control"))}
                 for k,v in cases.items()}
    points = {size: read_json(paths["synthetic"])["trace"] for size,paths in baselines.items()}
    _, _, sp, _ = source_inputs()
    return baselines,points,sp


def marker_specs(base, trace, sp):
    archive = read_json(MARKERS / "summary.json")
    dx = archive["mesh"]["dx_m"] * 3
    runtime = archive["requests"]["primary"][:4]
    results = {}
    for size, point in trace.items():
        d = sp["molecular_diffusivity_m2_s"] * point["reference_period_Peclet"] / point["period_Peclet"]
        for variant, dt in (("primary",DT),("half_step",DT/2)):
            name = size+"_"+variant
            results[name] = dict(D_m2_s=d,dt_s=dt,command=[*runtime,str(TRACKER/"track.jl"),
                *map(str,(MARKERS/"compact",base/(name+".f64"),dx,d,dt,STOP_S,COUNT,SEED))])
    return results


def execute_pair_laws(base, specs):
    # Each trajectory job uses ten threads. Run them sequentially.
    for name, spec in specs.items():
        timing = execute(spec["command"],base/(name+".log"),timeout=1800)
        write_json(base/(name+"_execution.json"),timing)
        print(name,timing,flush=True)


def paired_laws(base, specs, closure):
    delays, weights = np.loadtxt(ARCHIVE/"paths.csv",delimiter=",",skiprows=1,unpack=True)
    shares = np.array(closure["flow_fractions"])
    # Same 40 classes as the frozen closure, reconstructed from their flow shares.
    edges = np.r_[0.,np.cumsum(shares)]
    tau = closure["reference_liquid_transit_s"]
    time = DT*np.arange(21711)
    cases = {}
    curves = {}
    for name,spec in specs.items():
        pairs = np.fromfile(base/(name+".f64"),dtype="<f8").reshape(-1,2)
        valid = pairs[np.isfinite(pairs).all(axis=1)]
        recovery = len(valid)/COUNT
        if len(pairs)!=COUNT or recovery<.999:
            raise ValueError("Incomplete matched-Peclet marker sample")
        grouped_shares, grouped_delays, counts = group(valid,delays,weights,edges)
        np.testing.assert_allclose(grouped_shares,shares,rtol=0,atol=1e-15)
        np.testing.assert_allclose(grouped_delays/tau,closure["normalized_transit_times"],rtol=1e-12,atol=0)
        times = (grouped_delays[:,None]+grouped_delays).ravel()
        order = np.argsort(times)
        curves[name], metrics = summary(time,times[order],(counts/len(valid)).ravel()[order])
        cases[name] = dict(D_m2_s=spec["D_m2_s"],dt_s=spec["dt_s"],recovery=recovery,
                           t99=metrics["t99"],variance=metrics["variance"],counts=counts.tolist())
    results = {}
    for size in dict.fromkeys(k.split("_")[0] for k in specs):
        a,b = (cases[size+"_"+k] for k in ("primary","half_step"))
        change = dict(max_cdf=float(max(abs(curves[size+"_primary"]-curves[size+"_half_step"]))),
                      t99=abs(a["t99"]/b["t99"]-1),variance=abs(a["variance"]/b["variance"]-1))
        passed = change["max_cdf"]<=.005 and max(change["t99"],change["variance"])<=.05
        if not passed:
            raise ValueError(f"Matched-Peclet marker refinement failed: {size}, {change}")
        results[size] = dict(primary=a,half_step=b,refinement=change)
    return results


def prepare():
    _, trace, sp = inputs()
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")+"-"+uuid4().hex[:8]
    base = OUTPUT/"runs"/stamp
    base.mkdir(parents=True)
    specs = marker_specs(base,trace,sp)
    write_json(base/"marker_plan.json",dict(cases=specs,protocol_sha256=sha(HERE/"protocol.md"),
                                          preparation_source_sha256=sha(HERE/"markers.py")))
    print("Peclet marker preparation",base,flush=True)
    execute_pair_laws(base,specs)
    write_json(base/"marker_laws.json",paired_laws(base,specs,closure()))
    print("Marker laws frozen",base,flush=True)


if __name__=="__main__":
    prepare()
