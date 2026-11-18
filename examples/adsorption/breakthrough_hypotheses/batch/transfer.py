"""Test the qualified particle diffusivities in the unchanged N60 columns."""

import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from itertools import product
from pathlib import Path
import sys
from uuid import uuid4

import numpy as np

from examples.adsorption.breakthrough_hypotheses.batch import run as bath
from examples.adsorption.figure7_retention.run import checked_model, integrate
from examples.adsorption.porosity_matched_column import run as column
from examples.adsorption.transport_accuracy.run import validate_record
from pore2process_demo.bridge import runtime_versions
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT
from pore2process_demo.provenance import source_manifest

HERE = Path(__file__).resolve().parent
REFERENCE = column.HERE / "reference_summary.json"


def archived_inputs():
    record = read_json(REFERENCE)
    files = {key: ROOT / record[key]["path"] for key in ("column_inputs", "frozen_predictions")}
    for key, path in files.items():
        if sha(path) != record[key]["sha256"]:
            raise ValueError(f"Archived N60 column record changed: {key}")
    setup, freeze = (read_json(files[key]) for key in ("column_inputs", "frozen_predictions"))
    if freeze["plan_sha256"] != sha(files["column_inputs"]):
        raise ValueError("Archived N60 column plan differs from its freeze")
    if column.kernel.source() != (files["column_inputs"].parent / "solver.jl").read_text():
        raise ValueError("Current column kernel differs from the archived N60 kernel")
    return setup, freeze, files["column_inputs"].parent


def plans(batch_base):
    batch = bath.check(batch_base)
    if batch["figure7_fitted_parameters"] != 0:
        raise ValueError("Expected a batch-only qualification")
    original, freeze, _ = archived_inputs()
    result = {}
    for name in batch["selected"]:
        diffusivity = batch["qualification"][name]["diffusivity_m2_s"]
        cases = {}
        for flow, case in original["cases"].items():
            archived_levels = freeze["cases"][flow][column.PRIMARY]
            if not archived_levels["numerically_qualified"] or len(archived_levels["files"]) != 2:
                raise ValueError("Expected the qualified 48-shell N60 reference")
            cases[flow] = dict(case, parameters=dict(case["parameters"], diffusivity_m2_s=diffusivity),
                               variants={column.PRIMARY: case["variants"][column.PRIMARY]})
        result[name] = dict(cases=cases, primary=column.PRIMARY,
                            qualified_batch_sha256=sha(batch_base / "summary.json"),
                            archived_column_sha256=sha(REFERENCE))
    return result


def predict(base, batch_base):
    setup = plans(batch_base)
    for name, plan in setup.items():
        folder = base / name
        folder.mkdir()
        write_json(folder / "plan.json", plan)
    (base / "solver.jl").write_text(column.kernel.source())
    print("Frozen fixed-inventory column inputs", base, flush=True)

    def solve(job):
        name, flow = job
        folder = base / name
        command = [sys.executable, "-m", column.__name__, str(folder), "--variant", flow, column.PRIMARY]
        timing = column.execute(command, folder / f"{flow}.log", timeout=3600)
        print(name, flow, "seconds", round(timing["wall_s"], 1), flush=True)

    with ThreadPoolExecutor(max_workers=2) as workers:
        list(workers.map(solve, product(setup, column.FLOWS)))
    seal(base, batch_base, setup)


def seal(base, batch_base, setup):
    sealed = dict(batch_run=batch_base.relative_to(ROOT).as_posix(), solver_sha256=sha(base / "solver.jl"),
                  alternatives={})
    for name in setup:
        folder = base / name
        files = [folder / "plan.json"]
        for flow in column.FLOWS:
            levels = folder / f"{flow}_{column.PRIMARY}_levels.json"
            files += [levels, *(folder / f for f in read_json(levels)["files"])]
        sealed["alternatives"][name] = {p.name: sha(p) for p in files}
    write_json(base / "frozen.json", sealed)
    print("Predictions frozen", sha(base / "frozen.json"), flush=True)


def verified_curves(base):
    sealed = read_json(base / "frozen.json")
    setup = plans(ROOT / sealed["batch_run"])
    original, freeze, archived = archived_inputs()
    if (sealed["solver_sha256"] != sha(base / "solver.jl")
            or (base / "solver.jl").read_text() != column.kernel.source()
            or set(sealed["alternatives"]) != set(setup)):
        raise ValueError("Frozen solver or column alternatives changed")
    curves, numerical = {}, {}
    for flow in column.FLOWS:
        parent = freeze["cases"][flow][column.PRIMARY]["selected"]
        if sha(archived / parent) != freeze["curve_sha256"][parent]:
            raise ValueError("Archived primary prediction changed")
        curves[flow] = {"primary": read_json(archived / parent)}
        numerical[flow] = {"primary": freeze["cases"][flow][column.PRIMARY]["refinement"][-1]}
    for name, plan in setup.items():
        folder = base / name
        for filename, digest in sealed["alternatives"][name].items():
            if sha(folder / filename) != digest:
                raise ValueError(f"Frozen transfer artifact changed: {filename}")
        if read_json(folder / "plan.json") != plan:
            raise ValueError("A fixed-inventory transfer input changed")
        for flow, case in plan["cases"].items():
            p, spec = case["parameters"], case["variants"][column.PRIMARY]
            levels = read_json(folder / f"{flow}_{column.PRIMARY}_levels.json")
            if (not levels["numerically_qualified"] or not 2 <= len(levels["files"]) <= len(column.LEVELS)
                    or levels["selected"] != levels["files"][-1]):
                raise ValueError("Numerical refinement did not qualify")
            checks, previous = [], None
            for filename, (shells, step, early) in zip(levels["files"], column.LEVELS):
                saved = read_json(folder / filename)
                expected = dict(p, film_m_s=spec["film_m_s"], radial_cells=shells, dt_s=step,
                                startup_dt_s=early)
                if saved["inputs"] != expected or saved["scenario"] != spec["access"]:
                    raise ValueError("Saved column inputs differ from the frozen case")
                result = saved["result"]
                column.transfer.verify_run(result, expected, spec["access"], case)
                np.testing.assert_allclose(result["effective_conductance_m3_s"],
                    spec["film_matched_clean_bed_conductance_m3_s"], rtol=1e-12)
                if previous is not None:
                    checks.append(column.full_refinement(previous, result, p["inlet_concentration_kg_m3"]))
                previous = result
            if checks != levels["refinement"] or not (checks[-1]["passed"] and checks[-1]["startup"]["passed"]):
                raise ValueError("Frozen numerical refinement differs on replay")
            curves[flow][name] = saved
            numerical[flow][name] = checks[-1]
    return original, curves, numerical


def compare(base):
    original, curves, numerical = verified_curves(base)
    source = column.source
    data = source.markers()
    keys = {"f4_q233_l25": "q233", "f4_q625_l25": "q625", "f5_q233_l25": "q233"}
    parameters = {flow: case["parameters"] for flow, case in original["cases"].items()}
    feed = {flow: p["flow_m3_s"] * p["inlet_concentration_kg_m3"] for flow, p in parameters.items()}
    common_dose = min(data[key]["time"][-1] * feed[flow] / parameters[flow]["dry_mass_kg"]
                      for key, flow in keys.items() if key.startswith("f4"))
    comparisons = {}
    for key, flow in keys.items():
        p, markers = parameters[flow], data[key]
        mass, rate = p["dry_mass_kg"], feed[flow]
        common_time = min(float(common_dose * mass / rate), markers["time"][-1])
        observed = source.observation(markers, rate, mass, common_time)
        mid, stop = observed["t50_s"], observed["end_s"]
        time = np.r_[0., markers["time"]]
        fraction = np.r_[0., markers["fraction"]]
        observed["late_retained_kg_kg"] = (integrate(time, fraction, rate, stop)
                                           - integrate(time, fraction, rate, mid)) / mass
        scores = {}
        for name, saved in curves[flow].items():
            basis = dict(p, feed_rate_kg_s=rate)
            metrics = column.curve_metrics(saved, basis, markers, observed, common_time)
            model = checked_model(saved, basis)
            retained = np.interp([mid, stop], model["time_s"], model["retained_kg"])
            scores[name] = dict(**metrics, late_retained_kg_kg=float(retained[1] - retained[0]) / mass,
                diffusivity_m2_s=saved["inputs"]["diffusivity_m2_s"],
                radial_cells=saved["inputs"]["radial_cells"], time_step_s=saved["inputs"]["dt_s"])
        comparisons[key] = dict(observed=observed, predictions=scores)
    return dict(comparisons=comparisons, numerical=numerical, common_dose_kg_kg=float(common_dose))


def score(base):
    scores = compare(base)
    source = column.source
    files = [Path(__file__), Path(bath.__file__), HERE / "protocol.md", REFERENCE, Path(column.__file__),
             Path(column.kernel.__file__), Path(source.__file__), source.PIXELS, source.EXTRACTION,
             ROOT / "julia/adsorption/FilmParticle.jl", column.kernel.access.ORIGINAL]
    write_json(base / "summary.json", dict(schema="pore2process.breakthrough-particle-transfer.v1",
        frozen_sha256=sha(base / "frozen.json"), sources=dict(pixels=sha(source.PIXELS), paper=sha(bath.PDF)),
        source_sha256=source_manifest(files), versions=runtime_versions(), figure_parameters_fitted=0,
        artifacts_sha256={p.relative_to(base).as_posix(): sha(p) for p in base.rglob("*") if p.is_file()},
        score=scores))
    print(base / "summary.json", flush=True)


def run(batch_base):
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ") + "-" + uuid4().hex[:8]
    base = bath.OUTPUT / "transfer" / stamp
    base.mkdir(parents=True)
    predict(base, batch_base)
    score(base)


def check(base):
    summary = validate_record(base / "summary.json")
    if (summary["frozen_sha256"] != sha(base / "frozen.json")
            or summary["score"] != compare(base)):
        raise ValueError("Archived Diffusivity transfer differs on replay")
    print(dict(sources=len(summary["source_sha256"]), artifacts=len(summary["artifacts_sha256"])))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("base", type=Path, help="Frozen batch run, or column transfer to verify with --check")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.check:
        check(args.base.resolve())
    else:
        run(args.base.resolve())
