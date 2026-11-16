"""Compare frozen mixed and path-preserving repeat joins at fixed carbon mass."""
import argparse
from datetime import datetime, timezone
from functools import lru_cache
from pathlib import Path
from uuid import uuid4

import numpy as np

from examples.adsorption.alhamed.study import compare
from examples.adsorption.contact_column.run import crossings
from examples.adsorption.disordered_memory_column import run as parent
from examples.adsorption.exchange_column.run import CASES, PIXELS, source_inputs
from examples.adsorption.film_particle.metrics import column_passed
from examples.adsorption.shared_contact.run import full_refinement
from examples.adsorption.transport_accuracy.run import validate_record
from pore2process_demo.bridge import runtime, runtime_versions
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT
from pore2process_demo.provenance import source_manifest

HERE = Path(__file__).resolve().parent
OUTPUT = ROOT/"results/adsorption/disordered_repeat_join"
PARENT = parent.OUTPUT/"runs/20261116T155928Z-7a767fbb"
PARENT_SHA = "f9759fa3232d8b44a81e55e3f3f80afa8c480205d3d17c308d2d38d805ba31ff"
JOINS = ("mixed","preserved")


@lru_cache(maxsize=1)
def julia_model():
    parent.julia_models()
    jl = runtime()
    jl.include(str(HERE/"model.jl"))
    return jl.RepeatJoin


def solve(size,join,p):
    name = size+"_primary"
    if join=="mixed":
        return parent.solve(name,4,True,p)
    m = parent.transport(name,4)
    r = julia_model().solve_preserved(m.sections,m.flow,m.area,**p)
    jl = runtime()
    return {str(f):np.array(v,dtype=float).tolist() if bool(jl.isa(v,jl.AbstractArray)) else float(v)
            for f in jl.fieldnames(jl.typeof(r)) for v in (getattr(r,str(f)),)}


def archive():
    path = PARENT/"summary.json"
    if sha(path)!=PARENT_SHA:
        raise ValueError("Frozen mixed column archive changed")
    validate_record(path)
    return read_json(PARENT/"frozen_predictions.json")


def inert(size,join,p):
    p = dict(p,film_m_s=0.,inlet_concentration_kg_m3=1.,radial_cells=1,
             startup_s=0.,dt_s=.002,t_end_s=20.)
    r = solve(size,join,p)
    parent.verify(r,p)
    cdf = np.array(r["outlet_concentration_kg_m3"])[1:]
    weights = np.diff(cdf,prepend=0.)
    if min(weights)<-1e-9 or abs(cdf[-1]-1)>1e-8:
        raise ValueError("The inert full-bed transit distribution did not resolve")
    edges = np.array(r["time_s"])
    centres = (edges[1:]+edges[:-1])/2
    times = np.r_[centres[0]/2,(centres[1:]+centres[:-1])/2]
    mean = float(weights@times)
    var = float(weights@times**2-mean**2)
    model = parent.transport(size+"_primary",4)
    tau = p["liquid_volume_m3"]/p["flow_m3_s"]
    expected = julia_model().repeat_moments(model.sections,model.flow,p["blocks"],tau,join=="preserved")
    a,b = float(expected.mean),float(expected.variance)
    np.testing.assert_allclose(a,tau,rtol=1e-12,atol=0.)
    np.testing.assert_allclose(mean,a,rtol=.01,atol=0.)
    np.testing.assert_allclose(var,b,rtol=.10,atol=0.)
    return dict(mean_s=mean,variance_s2=var,analytical_mean_s=a,analytical_variance_s2=b,
                recovery=float(cdf[-1]),mass_residual=r["normalized_mass_residual"],dt_s=p["dt_s"])


def archived_parity(size,file,result):
    original = read_json(PARENT/file)["result"]
    for field in ("time_s","outlet_concentration_kg_m3","adsorbed_inventory_kg","liquid_inventory_kg"):
        np.testing.assert_allclose(result[field],original[field],rtol=1e-12,atol=1e-14,
                                   err_msg=f"{size} mixed-join archive {field}")


def refinement(size,base,p,reference):
    files = {key:[] for key in JOINS}
    reviews = {key:[] for key in JOINS}
    previous = {}
    feed = p["inlet_concentration_kg_m3"]
    for i,(radial,dt,early) in enumerate(parent.LEVELS):
        parameters = dict(p,radial_cells=radial,dt_s=dt,startup_dt_s=early)
        for join in JOINS:
            r = solve(size,join,parameters)
            parent.verify(r,parameters)
            if join=="mixed" and i<len(reference["refinement_files"]):
                item = reference["refinement_files"][i]
                if parameters!=read_json(PARENT/item)["inputs"]:
                    raise ValueError("The new control no longer matches the archived mixed inputs")
                archived_parity(size,item,r)
            path = base/f"{size}_{join}_{i}.json"
            write_json(path,dict(inputs=parameters,result=r))
            files[join].append(path.name)
            if i>0:
                reviews[join].append(full_refinement(previous[join],r,feed))
            previous[join] = r
            print(path.name,crossings(r,feed),"runtime s",round(r["wall_s"],2),flush=True)
        if i>0 and all(v[-1]["passed"] and v[-1]["startup"]["passed"] for v in reviews.values()):
            break
    return {join:dict(selected=files[join][-1],files=files[join],refinement=reviews[join],
        numerically_qualified=reviews[join][-1]["passed"] and reviews[join][-1]["startup"]["passed"])
        for join in JOINS}


def score(base,frozen):
    _,source,_,_ = source_inputs()
    if sha(PIXELS)!=source["pixel_csv_sha256"]:
        raise ValueError("The Figure 7 marker extraction changed")
    scores = {}
    for size,(series,_) in CASES.items():
        cases = {}
        for join,item in frozen["cases"][size].items():
            record = read_json(base/item["selected"])
            r,p = record["result"],record["inputs"]
            m = compare(r,series)
            cases[join] = dict(metrics=m,experimental_screen_passed=column_passed(m),
                numerical_screen_passed=item["numerically_qualified"],
                concentration_at_60_s=float(np.interp(60.,r["time_s"],r["outlet_concentration_kg_m3"])/
                                             p["inlet_concentration_kg_m3"]),
                normalized_mass_residual=r["normalized_mass_residual"],
                normalized_interface_flux_mismatch=r["normalized_interface_flux_mismatch"],
                dry_mass_kg=r["physical_adsorbent_mass_kg"],delay_scale=r["delay_scale"])
        base_mid = cases["mixed"]["metrics"]["t50"]["predicted_s"]
        path_mid = cases["preserved"]["metrics"]["t50"]["predicted_s"]
        exp_mid = cases["mixed"]["metrics"]["t50"]["observed_interpolated_s"]
        scores[size] = dict(**cases,midpoint_shift_s=path_mid-base_mid,
                            fraction_of_mixed_midpoint_error_removed=(base_mid-path_mid)/(base_mid-exp_mid))
    return scores


def run():
    reference = archive()
    bed,_ = parent.inputs()
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")+"-"+uuid4().hex[:8]
    base = OUTPUT/"runs"/stamp
    base.mkdir(parents=True)
    sources = source_manifest([HERE/name for name in ("run.py","model.jl","test_model.jl","protocol.md")]+[
        parent.HERE/name for name in ("run.py","model.jl","handoff.jl")])
    print(base,flush=True)
    cases,transport_checks = {},{}
    for size,setup in bed["cases"].items():
        p = setup["parameters"]
        transport_checks[size] = {join:inert(size,join,p) for join in JOINS}
        print(size,"inert",transport_checks[size],flush=True)
        cases[size] = refinement(size,base,p,reference["cases"][size])
    frozen = dict(cases=cases,parent_sha256=PARENT_SHA,inert=transport_checks,
                  curve_sha256={path.name:sha(path) for path in sorted(base.glob("d*.json"))})
    saved = base/"frozen_predictions.json"
    write_json(saved,frozen)
    digest = sha(saved)
    print("Curves frozen before reading Figure 7 markers",digest,flush=True)
    scores = score(base,read_json(saved))
    summary = base/"summary.json"
    write_json(summary,dict(schema="pore2process.disordered_repeat_join.v1",run_id=stamp,source_sha256=sources,
        artifacts_sha256={path.name:sha(path) for path in sorted(base.iterdir())},
        frozen_prediction_sha256=digest,figure7=scores,figure7_fitted_parameters=0,
        parent_sha256=PARENT_SHA,versions=runtime_versions()))
    write_json(OUTPUT/"latest.json",dict(summary=summary.relative_to(ROOT).as_posix()))
    print(summary,flush=True)


def check(path):
    archive()
    record = validate_record(path)
    base = path.parent
    saved = base/"frozen_predictions.json"
    frozen = read_json(saved)
    if (sha(saved)!=record["frozen_prediction_sha256"] or frozen["parent_sha256"]!=PARENT_SHA
            or any(sha(base/name)!=digest for name,digest in frozen["curve_sha256"].items())
            or record["figure7"]!=score(base,frozen)):
        raise ValueError("Frozen comparison changed")
    for size,case in frozen["cases"].items():
        for join,item in case.items():
            files = [read_json(base/name) for name in item["files"]]
            for r in files:
                parent.verify(r["result"],r["inputs"])
            feed = files[0]["inputs"]["inlet_concentration_kg_m3"]
            if [full_refinement(a["result"],b["result"],feed)
                    for a,b in zip(files,files[1:])]!=item["refinement"]:
                raise ValueError(f"Refinement changed for {size}/{join}")
    print(dict(sources=len(record["source_sha256"]),artifacts=len(record["artifacts_sha256"])))


if __name__=="__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--check",type=Path)
    args = p.parse_args()
    check(args.check.resolve()) if args.check else run()
