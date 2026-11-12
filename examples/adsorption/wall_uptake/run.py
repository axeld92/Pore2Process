"""Freeze wall-uptake films and score both shared-carbon column predictions."""
import argparse
from pathlib import Path
import shutil

import numpy as np

from examples.adsorption.alhamed.study import compare
from examples.adsorption.contact_column.run import SHORT_S, crossings, scoring_curve, startup_refinement
from examples.adsorption.exchange_column.run import CASES, PIXELS, source_inputs, validate_solution
from examples.adsorption.film_particle.metrics import column_passed
from examples.adsorption.shared_contact.run import LEVELS, STARTUP, full_refinement, review
from examples.adsorption.transport_accuracy.run import validate_record
from pore2process_demo.bridge import runtime_versions
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT
from pore2process_demo.provenance import source_manifest

from . import kernel
from .estimate import model
from .prepare import HERE, OUTPUT


def conservation(r,p,setup):
    validate_solution(r,p,setup["trace"])
    film = setup["coefficients"]
    area = film["external_area_m2"]
    np.testing.assert_allclose(r["external_particle_area_m2"],area,rtol=1e-12,atol=0)
    np.testing.assert_allclose(r["contact_area_by_path_m2"],area*np.array(setup["weights"]),rtol=1e-12,atol=0)
    np.testing.assert_allclose(r["film_coefficient_by_path_m_s"],film["path_film_m_s"],rtol=1e-12,atol=0)
    np.testing.assert_allclose(r["film_conductance_by_path_m3_s"],film["path_bed_conductance_m3_s"],rtol=1e-12,atol=0)
    np.testing.assert_allclose(r["total_film_conductance_m3_s"],film["total_bed_conductance_m3_s"],rtol=1e-12,atol=0)
    np.testing.assert_allclose(r["cumulative_film_uptake_kg"]-r["cumulative_film_release_kg"],
                               r["adsorbed_inventory_kg"][-1],rtol=1e-8,atol=1e-16)


def score(base,frozen,prepared,pixels_hash):
    if sha(PIXELS)!=pixels_hash:
        raise ValueError("Figure 7 digitization changed")
    result = {}
    for size,(series,_) in CASES.items():
        item = frozen["cases"][size]
        full = read_json(base/item["full"][-1])["result"]
        short = read_json(base/item["short"][-1])["result"] if item["short"] else None
        metrics = compare(scoring_curve(full,short),series)
        result[size] = dict(metrics=metrics,physical_screen_passed=column_passed(metrics),
            numerical_screen_passed=item["review"]["numerically_qualified"],
            physical_carbon_mass_kg=full["physical_adsorbent_mass_kg"],
            film_flux_residual=full["normalized_interface_flux_mismatch"],
            correlation_controls=prepared["sizes"][size]["controls"])
    return result


def run(base):
    if (base/"summary.json").exists():
        raise ValueError("A completed cohort exists. Use --check.")
    prepared = model(base)
    _,source,_,_ = source_inputs()
    files = [HERE/name for name in ("run.py","estimate.py","prepare.py","kernel.py","track_uptake.jl",
                                   "test_track_uptake.jl","test_estimate.py","protocol.md")]
    files += [kernel.JULIA/name for name in ("SharedContactColumn.jl","FilmParticle.jl","ExchangeColumn.jl")]
    files += [HERE.parent/name for name in ("carbon_contact/track_contacts.jl","path_exchange/track.jl",
                                           "shared_contact/run.py","contact_column/inputs.py",
                                           "exchange_column/run.py","alhamed/study.py")]
    manifest = source_manifest(files)
    for name in manifest:
        target = base/"source_snapshot"/name
        target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(ROOT/name,target)
    generated = base/"FilmByPathColumn.jl"
    generated.write_text(kernel.source())
    kernel_checks = kernel.verify(generated)
    write_json(base/"model.json",prepared)
    frozen = dict(model_sha256=sha(base/"model.json"),cases={})
    cl = prepared["closure"]
    print("Wall-uptake column predictions",base,flush=True)
    for size,setup in prepared["sizes"].items():
        item = dict(full=[],short=[])
        p0 = setup["parameters"]
        feed = p0["inlet_concentration_kg_m3"]
        print(size,setup["coefficients"]["ratio_to_correlation"],setup["step_changes"],
              setup["uptake"]["primary"]["errors"],flush=True)
        def save(p,group,index):
            r = kernel.solve(generated,cl,setup["mix"],setup["weights"],setup["coefficients"]["path_film_m_s"],p)
            conservation(r,p,setup)
            name = f"{size}_{group}{index}.json"
            write_json(base/name,dict(inputs=p,result=r))
            item[group].append(name)
            print(size,group,index,crossings(r,feed),"seconds",r["total_s"],flush=True)
            return r
        previous = None
        for i,(sub,nr,dt,early) in enumerate(LEVELS):
            p = {**p0,"subcells":sub,"radial_cells":nr,"dt_s":dt,"startup_dt_s":early}
            full = save(p,"full",i)
            if previous is not None:
                changes = full_refinement(previous,full,feed)
                print(size,changes,flush=True)
                if changes["passed"]:
                    break
            previous = full
        if not changes["startup"]["passed"]:
            previous = full
            for i,(sub,nr,dt) in enumerate(STARTUP):
                early = {**p,"subcells":sub,"radial_cells":nr,"startup_dt_s":dt,
                         "startup_s":SHORT_S,"t_end_s":SHORT_S}
                short = save(early,"short",i)
                if startup_refinement(previous,short,feed)["passed"]:
                    break
                previous = short
        item["review"] = review(base,item)
        item["hashes"] = {name:sha(base/name) for name in item["full"]+item["short"]}
        frozen["cases"][size] = item
    write_json(base/"frozen_predictions.json",frozen)
    frozen_hash = sha(base/"frozen_predictions.json")
    print("Predictions frozen before Figure 7 scoring",frozen_hash,flush=True)
    scored = score(base,frozen,prepared,source["pixel_csv_sha256"])
    record = dict(schema="pore2process.wall-uptake-film.v1",run_id=base.name,source_sha256=manifest,
        artifacts_sha256={str(p.relative_to(base)):sha(p) for p in base.rglob("*") if p.is_file()},
        frozen_prediction_sha256=frozen_hash,source_pixel_sha256=source["pixel_csv_sha256"],
        solver_checks=kernel_checks,figure7=scored,versions=runtime_versions(),
        new_column_fitted_parameters=0,new_OpenLB_runs=0,new_DWSIM_runs=0)
    path = base/"summary.json"
    write_json(path,record)
    write_json(OUTPUT/"latest.json",dict(summary=str(path.relative_to(ROOT))))
    print("Complete",path,flush=True)


def check(path):
    record = validate_record(path)
    base = path.parent
    prepared = model(base)
    frozen = read_json(base/"frozen_predictions.json")
    if (read_json(base/"model.json")!=prepared or sha(base/"model.json")!=frozen["model_sha256"]
            or sha(base/"frozen_predictions.json")!=record["frozen_prediction_sha256"]
            or (base/"FilmByPathColumn.jl").read_text()!=kernel.source()):
        raise ValueError("Frozen model or prediction manifest changed")
    for size,item in frozen["cases"].items():
        for name,digest in item["hashes"].items():
            if sha(base/name)!=digest:
                raise ValueError("Frozen prediction changed")
            run = read_json(base/name)
            conservation(run["result"],run["inputs"],prepared["sizes"][size])
        if review(base,item)!=item["review"]:
            raise ValueError("Numerical refinement differs on replay")
    if score(base,frozen,prepared,record["source_pixel_sha256"])!=record["figure7"]:
        raise ValueError("Figure 7 scores differ on replay")
    print(dict(sources=len(record["source_sha256"]),artifacts=len(record["artifacts_sha256"])))


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--prepared",type=Path)
    modes.add_argument("--check",type=Path)
    args = parser.parse_args()
    check(args.check.resolve()) if args.check else run(args.prepared.resolve())
