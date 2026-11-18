"""Calibrate one shared equilibrium weight on two depths and test three holdouts."""

import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from functools import lru_cache
from math import exp, log
from pathlib import Path
import sys
from uuid import uuid4

import numpy as np

from examples.adsorption.breakthrough_hypotheses.batch import freundlich, freundlich_transfer, transfer
from examples.adsorption.breakthrough_hypotheses.batch import run as bath
from examples.adsorption.figure7_retention.run import checked_model, integrate
from examples.adsorption.film_particle.bridge import batch, module
from examples.adsorption.porosity_matched_column import run as column
from examples.adsorption.transport_accuracy.run import validate_record
from pore2process_demo.bridge import runtime, runtime_versions
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT
from pore2process_demo.provenance import source_manifest

HERE = Path(__file__).resolve().parent
OUTPUT = bath.OUTPUT / "conditional"
REFERENCE = HERE.parent / "reference_summary.json"
WEIGHTS = {"w000": 0., "w025": .25, "w050": .50, "w075": .75, "w100": 1.}
CALIBRATION = {"q233_l10": "f5_q233_l10", "q233_l15": "f5_q233_l15"}
HOLDOUT = {"q233_l20": "f5_q233_l20", "q233_l25": "f4_q233_l25", "q625_l25": "f4_q625_l25"}
SURFACE_DISPATCH = """
function SharedContactColumn.surface_value(rhs,beta,qmax,b,law::EffectiveEquilibrium.Blend)
    upper = SharedContactColumn.surface_value(rhs,beta*(1-law.weight),qmax,b,nothing)
    EffectiveEquilibrium.surface_concentration(rhs,beta,qmax,b,upper,law)
end
"""


@lru_cache(maxsize=2)
def load(kind):
    if kind == "bath":
        module()
        freundlich.load()
    else:
        freundlich_transfer.load_column_law()
    jl = runtime()
    jl.include(str(HERE / "equilibrium.jl"))
    if kind == "column":
        jl.seval(SURFACE_DISPATCH)
    return jl.EffectiveEquilibrium


def law_at(weight, model):
    if weight in (0., 1.):
        return None if weight == 0. else freundlich.load().TABLE1
    return model.Blend(weight)


def equilibrium_loading(c, p, weight):
    bc = p["affinity_m3_kg"] * c
    return (1-weight) * p["capacity_kg_kg"] * bc/(1+bc) + weight * freundlich_transfer.loading_at_feed(c)


def verify_inversion(model):
    jl = runtime()
    p = dict(affinity_m3_kg=139., capacity_kg_kg=.1115)
    for weight in WEIGHTS.values():
        law = law_at(weight, model)
        for c in (0., 1e-250, 1e-9, .001, .0075, .075, .15):
            for beta in (0., .01, 1., 30., 1000., 1e6):
                rhs = c + beta * equilibrium_loading(c, p, weight)
                guess = jl.SharedContactColumn.surface_value(rhs, beta, p["capacity_kg_kg"],
                                                             p["affinity_m3_kg"], law)
                np.testing.assert_allclose(guess, c, rtol=1e-11)


def setup_profiles():
    refs = read_json(REFERENCE)["runs"]
    paths = {name: ROOT / refs[name]["path"] for name in ("langmuir_batch", "freundlich_batch")}
    for name, path in paths.items():
        if sha(path) != refs[name]["sha256"]:
            raise ValueError(f"Independently qualified endpoint input changed: {name}")
    langmuir = bath.check(paths["langmuir_batch"].parent)["qualification"]["primary"]
    friend = freundlich.check(paths["freundlich_batch"].parent)["qualification"]
    if not langmuir["qualifies"] or not friend["qualifies"]:
        raise ValueError("The two required source equilibrium fits must first pass the batch checks")
    variables = ("diffusivity_m2_s", "batch_film_ref_m_s")
    result = {}
    for name, w in WEIGHTS.items():
        values = {key: (langmuir[key] if w == 0. else friend[key] if w == 1. else
                        exp((1-w)*log(langmuir[key]) + w*log(friend[key]))) for key in variables}
        result[name] = dict(weight=w, **values)
    return result


def qualify_bath(base, profiles):
    directory = base / "batch"
    directory.mkdir()
    points = bath.data.observations()
    sizes = {key: dict(diameter_m=row["diameter_m"], role=row["role"]) for key, row in points.items()}
    audit = bath.inventory_and_equilibrium(bath.data.evidence())
    plan = dict(profiles=profiles, sizes=sizes, equilibrium=audit)
    write_json(directory / "plan.json", plan)
    model = load("bath")
    artifacts = {}
    for name, profile in profiles.items():
        law = law_at(profile["weight"], model)
        for size, data in sizes.items():
            for grid, (shells, step) in bath.RESOLUTIONS.items():
                p = bath.cases.batch_inputs(data["diameter_m"], profile["diffusivity_m2_s"],
                    profile["batch_film_ref_m_s"], radial_cells=shells, dt_s=step)
                result = batch(**dict(p, linear_partition_m3_kg=law))
                bath.verify_bath(result)
                path = directory / f"{name}_{size}_{grid}.json"
                write_json(path, dict(inputs=p, weight=profile["weight"], result=result))
                artifacts[path.name] = sha(path)
        print("Batch history", name, flush=True)
    write_json(directory / "frozen.json", dict(plan_sha256=sha(directory / "plan.json"), curve_sha256=artifacts))
    scores = bath.compare(directory, plan, points)
    write_json(base / "batch_scores.json", scores)
    print("Batch-qualified weights", [name for name, row in scores.items() if row["qualifies"]], flush=True)
    return scores


def depth_conditions():
    original, _, _ = transfer.archived_inputs()
    prior = column.bypass.COLUMN
    if sha(prior) != column.bypass.COLUMN_SHA:
        raise ValueError("The archived depth study changed")
    old_path = prior.parent / "plan.json"
    if sha(old_path) != read_json(prior)["artifacts_sha256"]["plan.json"]:
        raise ValueError("The archived depth plans changed")
    archived = read_json(old_path)["cases"]
    cases = {}
    for name in CALIBRATION | HOLDOUT:
        flow = name.split("_", 1)[0]
        parent = original["cases"][flow]
        depth = archived[name]
        source_spec = parent["variants"][column.PRIMARY]
        parameters = dict(parent["parameters"])
        for key in ("blocks", "terminal_sections", "dry_mass_kg", "liquid_volume_m3"):
            parameters[key] = depth["parameters"][key]
        fraction = parameters["dry_mass_kg"] / parent["parameters"]["dry_mass_kg"]
        np.testing.assert_allclose([depth["Pe"], depth["flow_mL_min"],
                                    depth["scenarios"]["f50_slow01"]["added_resistance_s_m"]],
                                   [parent["Pe"], parent["flow_mL_min"], parent["added_resistance_s_m"]], rtol=1e-13)
        np.testing.assert_allclose(fraction, depth["sections"]/parent["sections"], rtol=1e-13)
        spec = dict(source_spec, film_matched_clean_bed_conductance_m3_s=
                    fraction * source_spec["film_matched_clean_bed_conductance_m3_s"])
        cases[name] = dict(parameters=parameters, sections=depth["sections"],
                           flow_mL_min=parent["flow_mL_min"], spec=spec)
    return cases


def verify_column(result, p, spec, condition, weight):
    feed = p["inlet_concentration_kg_m3"]
    loading = equilibrium_loading(feed, p, weight)
    bc = p["affinity_m3_kg"] * feed
    # Reuse the historical checks with the stated feed loading as ceiling.
    ceiling_basis = dict(p, capacity_kg_kg=loading * (1 + bc) / bc)
    column.transfer.verify_run(result, ceiling_basis, spec["access"], condition)
    if result["loading_range_kg_kg"][1] > loading + 1e-10:
        raise ValueError("A particle exceeded the mixture's equilibrium loading at feed")
    np.testing.assert_allclose(result["effective_conductance_m3_s"],
                               spec["film_matched_clean_bed_conductance_m3_s"], rtol=1e-12)


def solve(base, name, case):
    plan = read_json(base / "plan.json")
    profile = plan["qualified_profiles"][name]
    weight = profile["weight"]
    condition = plan["conditions"][case]
    spec = condition["spec"]
    p = dict(condition["parameters"], diffusivity_m2_s=profile["diffusivity_m2_s"])
    model = load("column")
    verify_inversion(model)
    liquid_equilibrium = law_at(weight, model)
    previous, files = None, []
    for index, level in enumerate(column.LEVELS):
        inputs = freundlich_transfer.level_inputs(p, spec, level)
        result = column.kernel.solve(condition["flow_mL_min"], ROOT / spec["probe_directory"],
            dict(inputs, linear_partition_m3_kg=liquid_equilibrium), spec["access"])
        verify_column(result, inputs, spec, condition, weight)
        filename = f"{name}_{case}_{index}.json"
        write_json(base / filename, dict(inputs=inputs, result=result, scenario=spec["access"], weight=weight))
        files.append(filename)
        print(case, name, index, "wall seconds", round(result["wall_s"], 1), flush=True)
        if previous is not None:
            refinement = column.full_refinement(previous, result, inputs["inlet_concentration_kg_m3"])
            if refinement["passed"] and refinement["startup"]["passed"]:
                write_json(base / f"{name}_{case}_levels.json",
                           dict(files=files, selected=filename, refinement=refinement))
                return
        previous = result
    raise ValueError(f"{name} {case} did not meet the numerical refinement screens")


def predict(base, profiles, cases, stage):
    def worker(job):
        profile, case = job
        command = [sys.executable, "-m", __spec__.name, str(base), "--case", profile, case]
        timing = column.execute(command, base / f"{profile}_{case}.log", timeout=3600)
        print(profile, case, "seconds", round(timing["wall_s"], 1), flush=True)

    jobs = [(name, case) for name in profiles for case in cases]
    with ThreadPoolExecutor(max_workers=3) as workers:
        list(workers.map(worker, jobs))
    files = [base / "plan.json"]
    for name, case in jobs:
        levels = base / f"{name}_{case}_levels.json"
        files.extend([levels, *(base / file for file in read_json(levels)["files"])])
    if stage == "holdout":
        files.append(base / "selection.json")
    write_json(base / f"frozen_{stage}.json", dict(artifacts={p.name: sha(p) for p in files}))
    print("Sealed", stage, sha(base / f"frozen_{stage}.json"), flush=True)


def score_one(base, candidate, case, markers, plan):
    condition = plan["conditions"][case]
    levels = read_json(base / f"{candidate}_{case}_levels.json")
    previous = None
    spec = condition["spec"]
    profile = plan["qualified_profiles"][candidate]
    w = profile["weight"]
    initial = dict(condition["parameters"], diffusivity_m2_s=profile["diffusivity_m2_s"])
    if not 2 <= len(levels["files"]) <= len(column.LEVELS) or levels["selected"] != levels["files"][-1]:
        raise ValueError("Incomplete spatial and temporal refinement")
    for filename, level in zip(levels["files"], column.LEVELS):
        saved = read_json(base / filename)
        p = freundlich_transfer.level_inputs(initial, spec, level)
        if saved["weight"] != w or saved["inputs"] != p or saved["scenario"] != spec["access"]:
            raise ValueError("An equilibrium or fixed column input changed")
        result = saved["result"]
        verify_column(result, p, spec, condition, w)
        refinement = column.full_refinement(previous, result, p["inlet_concentration_kg_m3"]) if previous else None
        previous = result
    if refinement != levels["refinement"] or not (refinement["passed"] and refinement["startup"]["passed"]):
        raise ValueError("Frozen refinement differs on replay")
    rate, mass = p["flow_m3_s"] * p["inlet_concentration_kg_m3"], p["dry_mass_kg"]
    end = markers["time"][-1]
    observed = column.source.observation(markers, rate, mass, end)
    basis = dict(p, feed_rate_kg_s=rate)
    predicted = column.curve_metrics(saved, basis, markers, observed, end)
    model = checked_model(saved, basis)
    mid = observed["t50_s"]
    if mid is None:
        observed["late_retained_kg_kg"] = predicted["late_retained_kg_kg"] = None
    else:
        time, concentration = np.r_[0., markers["time"]], np.r_[0., markers["fraction"]]
        observed["late_retained_kg_kg"] = (integrate(time, concentration, rate, end)
                                           - integrate(time, concentration, rate, mid)) / mass
        r = np.interp([mid, end], model["time_s"], model["retained_kg"])
        predicted["late_retained_kg_kg"] = float(r[1] - r[0]) / mass
    return dict(observed=observed, predicted=predicted, refinement=levels["refinement"],
                radial_cells=p["radial_cells"], time_step_s=p["dt_s"])


def compare(base, names, stage):
    sealed = read_json(base / f"frozen_{stage}.json")["artifacts"]
    if any(sha(base / file) != digest for file, digest in sealed.items()):
        raise ValueError(f"Frozen {stage} predictions changed")
    plan = read_json(base / "plan.json")
    if plan["conditions"] != depth_conditions():
        raise ValueError("The fixed depth assumptions changed")
    markers = column.source.markers()
    cases = CALIBRATION if stage == "calibration" else HOLDOUT
    scores = {name: {case: score_one(base, name, case, markers[key], plan)
                     for case, key in cases.items()} for name in names}
    if stage == "holdout":
        for name in names:
            scores[name]["q233_l25_Figure5_duplicate"] = score_one(
                base, name, "q233_l25", markers["f5_q233_l25"], plan)
    return scores


def select(calibration):
    errors = {key: sum((row["predicted"]["retained_kg_kg"] - row["observed"]["retained_kg_kg"])**2
                       for row in beds.values()) for key, beds in calibration.items()}
    winner = min(errors, key=errors.get)
    return dict(candidate=winner, weight=WEIGHTS[winner], squared_retention_error_by_weight=errors,
                calibration=["f5_q233_l10", "f5_q233_l15"])


def run():
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ") + "-" + uuid4().hex[:8]
    base = OUTPUT / stamp
    base.mkdir(parents=True)
    print(base, flush=True)
    profiles = setup_profiles()
    scores = qualify_bath(base, profiles)
    eligible = {name: profile for name, profile in profiles.items() if scores[name]["qualifies"]}
    if not eligible:
        raise ValueError("No proposed mixture passed the independent batch screens")
    write_json(base / "plan.json", dict(conditions=depth_conditions(), qualified_profiles=eligible,
                                       batch_scores_sha256=sha(base / "batch_scores.json"),
                                       batch_frozen_sha256=sha(base / "batch/frozen.json")))
    predict(base, eligible, CALIBRATION, "calibration")
    calibration = compare(base, eligible, "calibration")
    selected = select(calibration)
    write_json(base / "selection.json", selected)
    print("Frozen calibration", selected, flush=True)
    predict(base, [selected["candidate"]], HOLDOUT, "holdout")
    holdout = compare(base, [selected["candidate"]], "holdout")
    files = [Path(__file__), HERE / "equilibrium.jl", HERE / "protocol.md", REFERENCE, bath.PDF,
             Path(transfer.__file__), Path(freundlich_transfer.__file__), freundlich.JULIA, bath.data.PIXELS,
             column.source.PIXELS, column.source.EXTRACTION, Path(column.kernel.__file__),
             column.kernel.access.ORIGINAL]
    write_json(base / "summary.json", dict(schema="pore2process.breakthrough-effective-equilibrium.v1",
        calibration=calibration, selection=selected, holdout=holdout, figure7_fitted_parameters=0,
        fitted_shared_figure5_parameters=1, versions=runtime_versions(), source_sha256=source_manifest(files),
        artifacts_sha256={p.relative_to(base).as_posix(): sha(p) for p in base.rglob("*") if p.is_file()}))
    print(base / "summary.json", flush=True)


def check(base):
    result = validate_record(base / "summary.json")
    plan = read_json(base / "plan.json")
    batch_plan = read_json(base / "batch/plan.json")
    scores = bath.compare(base / "batch", batch_plan, bath.data.observations())
    if (scores != read_json(base / "batch_scores.json") or
            plan["batch_scores_sha256"] != sha(base / "batch_scores.json") or
            plan["batch_frozen_sha256"] != sha(base / "batch/frozen.json")):
        raise ValueError("Conditional equilibrium batch scores differ on replay")
    calibration = compare(base, plan["qualified_profiles"], "calibration")
    selection = select(calibration)
    holdout = compare(base, [selection["candidate"]], "holdout")
    if (calibration != result["calibration"] or selection != read_json(base / "selection.json") or
            selection != result["selection"] or holdout != result["holdout"]):
        raise ValueError("Calibrated predictions or scores differ on replay")
    print(dict(sources=len(result["source_sha256"]), artifacts=len(result["artifacts_sha256"])))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("base", type=Path, nargs="?")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--case", nargs=2)
    mode.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if (args.case or args.check) and args.base is None:
        parser.error("Pass the conditional equilibrium run directory")
    if args.case:
        solve(args.base.resolve(), *args.case)
    elif args.check:
        check(args.base.resolve())
    else:
        run()
