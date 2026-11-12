"""Freeze contact-based exchange predictions and replay the Figure 7 comparison."""
import argparse
from datetime import datetime, timezone
from pathlib import Path
import shutil
from uuid import uuid4

import numpy as np

from examples.adsorption.alhamed.study import compare, crossing
from examples.adsorption.exchange_attribution import kernel
from examples.adsorption.exchange_attribution.markers import inputs as baseline_inputs
from examples.adsorption.exchange_attribution.run import allocation_check
from examples.adsorption.exchange_column.run import CASES, PIXELS, source_inputs
from examples.adsorption.film_particle.metrics import column_passed
from examples.adsorption.transport_accuracy.run import validate_record
from pore2process_demo.bridge import runtime_versions
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT
from pore2process_demo.provenance import source_manifest

from .inputs import model

HERE = Path(__file__).parent
OUTPUT = ROOT/"results/adsorption/contact_column"
LEVELS = ((2,48,4.,.04),(4,96,2.,.02),(6,128,1.5,.015))
STARTUP_LEVELS = ((8,128,.01),(12,128,.0075),(16,128,.005))
SHORT_S = 12.
THRESHOLDS = {"t05":.05,"t50":.5}


def coordinates(result,feed):
    return np.array(result["time_s"]),np.array(result["outlet_concentration_kg_m3"])/feed


def crossings(result,feed):
    t,y = coordinates(result,feed)
    return {k:crossing(t,y,f) for k,f in THRESHOLDS.items()}


def startup_refinement(coarse,fine,feed):
    a,b = (crossings(r,feed) for r in (coarse,fine))
    changes = {k:(abs(a[k]/b[k]-1) if a[k] is not None and b[k] is not None else None)
               for k in THRESHOLDS if min(a[k] or np.inf,b[k] or np.inf)<=SHORT_S}
    return dict(relative=changes,passed=all(v is not None and v<=.02 for v in changes.values()))


def refine(coarse,fine,feed):
    (t,c),(tf,cf) = (coordinates(r,feed) for r in (coarse,fine))
    mask = t>=24
    change = abs(c-np.interp(t,tf,cf))
    errors = dict(max_after_startup=float(max(change[mask])),max_startup=float(max(change[~mask])))
    a,b = (crossings(r,feed) for r in (coarse,fine))
    for k in THRESHOLDS:
        if min(a[k] or np.inf,b[k] or np.inf)>SHORT_S:
            errors[k+"_relative"] = abs(a[k]/b[k]-1) if a[k] is not None and b[k] is not None else None
    passed = errors["max_after_startup"]<=.005 and all(v is not None and v<=.02 for k,v in errors.items()
                                                       if k.endswith("_relative"))
    return dict(errors=errors,passed=passed,startup=startup_refinement(coarse,fine,feed))


def join_difference(full,short,feed):
    t,c = coordinates(full,feed)
    _,cs = coordinates(short,feed)
    return float(abs(np.interp(SHORT_S,t,c)-cs[-1]))


def scoring_curve(result,early):
    if early is None:
        return result
    # Only the concentration profile is used for observation scoring.
    # Both mass histories remain those of their individual solver runs.
    fields = ("time_s","outlet_concentration_kg_m3")
    tail = np.array(result["time_s"])>SHORT_S
    return {**result,**{k:early[k]+np.array(result[k])[tail].tolist() for k in fields}}


def score(base, frozen, prepared, pixels_hash):
    if sha(PIXELS)!=pixels_hash:
        raise ValueError("Figure 7 digitization changed")
    summary = {}
    for size,(series,_) in CASES.items():
        out = {}
        for rank,item in frozen["cases"][size].items():
            case = read_json(base/item["selected"])
            short = read_json(base/item["startup_files"][-1])["result"] if item["startup_files"] else None
            metrics = compare(scoring_curve(case["result"],short),series)
            out[rank] = dict(metrics=metrics,physical_screen_passed=column_passed(metrics),
                             numerical_refinement_passed=item["numerically_qualified"],
                             total_mass_kg=case["result"]["physical_adsorbent_mass_kg"],
                             film_flux_residual=case["result"]["normalized_interface_flux_mismatch"])
        out["void_proportional_control"] = prepared["sizes"][size]["void_proportional_control"]
        summary[size] = out
    return summary


def run():
    prepared = model()
    baselines,points,_ = baseline_inputs()
    _, source, _, _ = source_inputs()
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")+"-"+uuid4().hex[:8]
    base = OUTPUT/"runs"/stamp
    base.mkdir(parents=True)
    files = [HERE/name for name in ("run.py","inputs.py","protocol.md")]
    files += [kernel.ORIGINAL,kernel.JULIA/"FilmParticle.jl",HERE.parent/"exchange_attribution/kernel.py",
              HERE.parent/"exchange_attribution/run.py",HERE.parent/"carbon_contact/analyze.py",
              HERE.parent/"carbon_contact/geometry.py",HERE.parent/"alhamed/study.py",
              HERE.parent/"exchange_column/closure.py",HERE.parent/"exchange_column/run.py"]
    manifest = source_manifest(files)
    for name in manifest:
        target = base/"source_snapshot"/name
        target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(ROOT/name,target)
    generated = base/"AllocationColumn.jl"
    generated.write_text(kernel.source())
    verified = kernel.verify(generated)
    write_json(base/"model.json",prepared)
    frozen = dict(model_sha256=sha(base/"model.json"),cases={})
    cl = prepared["closure"]
    print("Contact-allocation runs",base,verified,flush=True)
    for size, setup in prepared["sizes"].items():
        prior = read_json(baselines[size]["synthetic"])["inputs"]
        frozen["cases"][size] = {}
        for rank,shares in setup["allocations"].items():
            levels,refinements = [],[]
            previous = None
            allocation = dict(closure=cl,allocations={rank:shares})
            for index,(subcells,radial,dt,startup) in enumerate(LEVELS):
                p = {**prior,"subcells":subcells,"radial_cells":radial,"dt_s":dt,"startup_dt_s":startup}
                r = kernel.solve(generated,cl["flow_fractions"],cl["normalized_transit_times"],
                                 setup["mix"],shares,**p)
                allocation_check(r,p,points[size],allocation,rank)
                filename = f"{size}_{rank}_level{index}.json"
                write_json(base/filename,dict(inputs=p,trace=points[size],allocation=rank,result=r))
                levels.append(filename)
                print(size,rank,index,"seconds",r["total_s"],"mass residual",r["normalized_mass_residual"],flush=True)
                if previous is not None:
                    refinements.append(refine(previous,r,p["inlet_concentration_kg_m3"]))
                    print(size,rank,refinements[-1],flush=True)
                    if refinements[-1]["passed"]:
                        break
                previous = r
            full = r
            feed = p["inlet_concentration_kg_m3"]
            startups = []
            startup_changes = [refinements[-1]["startup"]]
            join = None
            if not startup_changes[-1]["passed"]:
                previous = full
                for index,(subcells,radial,dt) in enumerate(STARTUP_LEVELS):
                    short_p = {**p,"subcells":subcells,"radial_cells":radial,"startup_dt_s":dt,
                               "startup_s":SHORT_S,"t_end_s":SHORT_S}
                    r = kernel.solve(generated,cl["flow_fractions"],cl["normalized_transit_times"],setup["mix"],shares,**short_p)
                    allocation_check(r,short_p,points[size],allocation,rank)
                    filename = f"{size}_{rank}_startup{index}.json"
                    write_json(base/filename,dict(inputs=short_p,trace=points[size],allocation=rank,result=r))
                    startups.append(filename)
                    startup_changes.append(startup_refinement(previous,r,feed))
                    join = join_difference(full,r,feed)
                    print(size,rank,"startup",index,startup_changes[-1],"join",join,flush=True)
                    if startup_changes[-1]["passed"]:
                        break
                    previous = r
            numerical = refinements[-1]["passed"] and startup_changes[-1]["passed"] and (join is None or join<=.005)
            frozen["cases"][size][rank] = dict(selected=levels[-1],refinement=refinements,
                startup_files=startups,startup_refinement=startup_changes,join_fraction_change=join,
                numerically_qualified=numerical,files={name:sha(base/name) for name in levels+startups})
    write_json(base/"frozen_predictions.json",frozen)
    frozen_hash = sha(base/"frozen_predictions.json")
    print("Predictions frozen before Figure 7 scoring",frozen_hash,flush=True)
    scored = score(base,frozen,prepared,source["pixel_csv_sha256"])
    record = dict(schema="pore2process.contact-column.v1",run_id=stamp,source_sha256=manifest,
        artifacts_sha256={str(p.relative_to(base)):sha(p) for p in base.rglob("*") if p.is_file()},
        frozen_prediction_sha256=frozen_hash,source_pixel_sha256=source["pixel_csv_sha256"],
        solver_checks=verified,figure7=scored,versions=runtime_versions(),new_column_fitted_parameters=0,
        new_OpenLB_runs=0,new_DWSIM_runs=0)
    path = base/"summary.json"
    write_json(path,record)
    write_json(OUTPUT/"latest.json",dict(summary=str(path.relative_to(ROOT))))
    print("Complete",path,flush=True)


def check(path):
    record = validate_record(path)
    base = path.parent
    prepared = model()
    frozen = read_json(base/"frozen_predictions.json")
    if (sha(base/"frozen_predictions.json")!=record["frozen_prediction_sha256"]
            or read_json(base/"model.json")!=prepared or sha(base/"model.json")!=frozen["model_sha256"]
            or (base/"AllocationColumn.jl").read_text()!=kernel.source()):
        raise ValueError("The frozen model changed")
    for size,cases in frozen["cases"].items():
        for rank,item in cases.items():
            results = {}
            allocation = dict(closure=prepared["closure"],allocations=prepared["sizes"][size]["allocations"])
            for file,digest in item["files"].items():
                if sha(base/file)!=digest:
                    raise ValueError("Frozen prediction changed")
                case = read_json(base/file)
                r,p = case["result"],case["inputs"]
                allocation_check(r,p,case["trace"],allocation,rank)
                results[file] = r
            feed = p["inlet_concentration_kg_m3"]
            levels = [r for name,r in results.items() if name not in item["startup_files"]]
            refinements = [refine(a,b,feed) for a,b in zip(levels,levels[1:])]
            startup = [results[name] for name in item["startup_files"]]
            seq = [levels[-1],*startup]
            changes = [refinements[-1]["startup"]]+[startup_refinement(a,b,feed) for a,b in zip(seq,seq[1:])]
            join = join_difference(levels[-1],startup[-1],feed) if startup else None
            if refinements!=item["refinement"] or changes!=item["startup_refinement"] or join!=item["join_fraction_change"]:
                raise ValueError("Refinement comparison differs on replay")
    if score(base,frozen,prepared,record["source_pixel_sha256"])!=record["figure7"]:
        raise ValueError("Figure 7 scores differ on replay")
    print(dict(sources=len(record["source_sha256"]),artifacts=len(record["artifacts_sha256"])))


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check",type=Path)
    args = parser.parse_args()
    check(args.check.resolve()) if args.check else run()
