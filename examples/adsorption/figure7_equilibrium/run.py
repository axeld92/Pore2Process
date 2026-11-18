"""Test the fixed equilibrium cases against both Figure 7 particle sizes."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from itertools import product
from pathlib import Path
import sys
from uuid import uuid4

import numpy as np

from examples.adsorption.alhamed.source import PDF, PDF_SHA
from examples.adsorption.alhamed.study import compare
from examples.adsorption.breakthrough_hypotheses.batch import run as bath
from examples.adsorption.breakthrough_hypotheses.batch.conditional import run as equilibrium
from examples.adsorption.carbon_access import kernel, run as access
from examples.adsorption.disordered_memory_column import run as column
from examples.adsorption.dispersion_holdout.run import execute
from examples.adsorption.exchange_column.run import CASES, PIXELS
from examples.adsorption.figure7_retention import run as retention
from examples.adsorption.film_particle.metrics import batch_metrics, column_passed
from examples.adsorption.shared_contact.run import full_refinement
from examples.adsorption.transport_accuracy.run import validate_record
from pore2process_demo.bridge import runtime_versions
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT
from pore2process_demo.provenance import source_manifest

HERE = Path(__file__).resolve().parent
OUTPUT = ROOT / "results/adsorption/figure7_equilibrium"
BATCH_REFERENCE = equilibrium.REFERENCE
BLEND_REFERENCE = equilibrium.HERE / "reference_summary.json"
ACCESS_REFERENCE = access.HERE / "reference_summary.json"
PUMP_REFERENCE = HERE.parent / "column_transport/reference/process_reference.json"
GAC4 = dict(capacity_kg_kg=.103, affinity_m3_kg=173.7820329925531)


def recorded(path):
    return dict(path=path.relative_to(ROOT).as_posix(), sha256=sha(path))


def gac4_batch(profile):
    """Replay the fixed-kinetics check for the GAC4-informed equilibrium."""
    scores = {}
    for size, observed in bath.data.observations().items():
        levels = {}
        for grid, (shells, step) in bath.RESOLUTIONS.items():
            p = bath.cases.batch_inputs(observed["diameter_m"], profile["diffusivity_m2_s"],
                profile["batch_film_ref_m_s"], radial_cells=shells, dt_s=step)
            result = bath.batch(**dict(p, **GAC4))
            bath.verify_bath(result)
            levels[grid] = result
        nominal = levels.pop("base")
        target = bath.loading_at(nominal, observed["time_s"])
        scale = nominal["mean_loading_kg_kg"][-1]
        change = {grid: float(max(abs(bath.loading_at(result, observed["time_s"]) - target)) / scale)
                  for grid, result in levels.items()}
        metrics = batch_metrics(nominal, observed)
        if not metrics["RMSE_screen_passed"] or max(change.values()) > .02:
            raise ValueError(f"GAC4 particle check failed for {size}")
        scores[size] = dict(metrics, relative_refinements=change)
    return scores


def materials():
    profiles = equilibrium.setup_profiles()
    descriptor = read_json(BLEND_REFERENCE)
    archive = ROOT / descriptor["run"]["path"]
    if sha(archive) != descriptor["run"]["sha256"]:
        raise ValueError("The frozen Figure 5 equilibrium selection changed")
    frozen = archive.parent / "selection.json"
    if sha(frozen) != descriptor["run"]["frozen_selection_sha256"]:
        raise ValueError("The archived selection changed")
    selected = read_json(frozen)["candidate"]
    unblended = dict(capacity_kg_kg=.1115, affinity_m3_kg=139.)
    result = {
        "gac4": dict(profiles["w000"], **GAC4),
        "langmuir": dict(profiles["w000"], **unblended),
        "blend": dict(profiles[selected], **unblended),
        "freundlich": dict(profiles["w100"], **unblended),
    }
    for p in result.values():
        p["loading_at_feed_kg_kg"] = equilibrium.equilibrium_loading(.15, p, p["weight"])
    np.testing.assert_allclose(result["blend"]["loading_at_feed_kg_kg"],
                               descriptor["selection"]["loading_at_feed_kg_kg"], rtol=1e-13)
    return result, recorded(archive)


def pump_inputs():
    """Use the actual saved DWSIM shaft powers at the Figure 7 flow."""
    export = read_json(PUMP_REFERENCE)
    path = ROOT / export["source_summary"]
    if sha(path) != export["source_summary_sha256"] or not export["native_to_Julia_to_DWSIM_executed"]:
        raise ValueError("Expected the archived actual DWSIM process run")
    summary = read_json(path)
    native = ROOT / export["native_summary"]
    if sha(native) != export["native_summary_sha256"] or summary["native_summary"] != export["native_summary"]:
        raise ValueError("The archived hydraulic geometry changed")
    geometry = read_json(native)["inputs"]
    pumps = {}
    for size, row in export["cases"].items():
        saved = path.parent / size / "dwsim/summary.json"
        pump = read_json(saved)
        if (pump != summary["cases"][size]["dwsim"] or not all(row["checks"].values())
                or not pump["save_reload_solve_checked"] or row["shaft_power_W"] != pump["shaft_power_W"]):
            raise ValueError(f"Unverified DWSIM pump output for {size}")
        pumps[size] = dict(source=recorded(saved), flow_m3_s=row["julia_output"]["flow_m3_s"],
            bed_length_m=row["julia_inputs"]["length_m"], shaft_power_W=pump["shaft_power_W"],
            pressure_drop_Pa=row["julia_output"]["pressure_drop_Pa"], dry_mass_kg=geometry["dry_mass_kg"],
            nominal_porosity=geometry["nominal_porosity"])
    return pumps, recorded(path), recorded(native)


def plan():
    if sha(PDF) != PDF_SHA:
        raise ValueError("The published Figure 7 source changed")
    bed, _ = column.inputs()
    reference = read_json(ACCESS_REFERENCE)
    archive = ROOT / reference["path"]
    if sha(archive) != reference["sha256"] or reference["figure7_fitted_parameters"] != 0:
        raise ValueError("Expected the archived unfitted access scenarios")
    old = archive.parent / "scenario_plan.json"
    if sha(old) != reference["scenario_plan_sha256"] or read_json(old)["bed"] != bed:
        raise ValueError("The archived bed or access assumptions changed")
    initial_access = read_json(old)["scenarios"]["f50_slow01"]
    reference_film = bed["cases"]["d800"]["parameters"]["film_m_s"]
    resistance = (1 / initial_access["slow_ratio"] - 1) / reference_film
    profiles, material_record = materials()
    pumps, pump_record, hydraulic_record = pump_inputs()
    cases = bed["cases"]
    basis = read_json(retention.REFERENCE)["basis"]
    for size, condition in cases.items():
        p = condition["parameters"]
        for quantity in ("dry_mass_kg", "liquid_volume_m3", "flow_m3_s", "inlet_concentration_kg_m3"):
            if p[quantity] != basis[quantity]:
                raise ValueError(f"Figure 7 basis changed for {size}/{quantity}")
        np.testing.assert_allclose([p["flow_m3_s"], p["dry_mass_kg"], condition["column_section_Pe"]],
            [pumps[size]["flow_m3_s"], pumps[size]["dry_mass_kg"], condition["probe_Pe"]], rtol=1e-11)
        condition["access"] = dict(direct_fraction=initial_access["direct_fraction"],
                                   slow_ratio=1/(1 + resistance*p["film_m_s"]))
    return dict(cases=cases, basis=basis, profiles=profiles, pumps=pumps,
        added_resistance_s_m=resistance, archived_access=recorded(archive),
        archived_tracer=recorded(column.MEMORY / "summary.json"), archived_pumps=pump_record,
        archived_hydraulic_geometry=hydraulic_record,
        archived_blend=material_record, figure7_retention_reference=recorded(retention.REFERENCE),
        figure7_fitted_parameters=0, numerical_levels=list(map(list, column.LEVELS)))


def inputs(condition, profile, level):
    radial, step, early = level
    changes = {k: profile[k] for k in ("capacity_kg_kg", "affinity_m3_kg", "diffusivity_m2_s")}
    return dict(condition["parameters"], **changes, radial_cells=radial, dt_s=step, startup_dt_s=early)


def verify(result, p, condition, profile):
    feed = p["inlet_concentration_kg_m3"]
    loading = equilibrium.equilibrium_loading(feed, p, profile["weight"])
    np.testing.assert_allclose(loading, profile["loading_at_feed_kg_kg"], rtol=1e-13)
    bc = p["affinity_m3_kg"] * feed
    # The existing cohort check takes a Langmuir input. Match its feed ceiling
    # to the actual equilibrium law without changing the simulation input.
    access.verify(result, dict(p, capacity_kg_kg=loading*(1+bc)/bc), condition["access"])
    if result["loading_range_kg_kg"][1] > loading + 1e-10 or result["sections"] != 8*p["blocks"]:
        raise ValueError("Loading or section count is inconsistent with the frozen case")


def refinement(base, size, name, saved=None):
    setup = read_json(base / "plan.json")
    condition, profile = setup["cases"][size], setup["profiles"][name]
    transport = condition["access"]
    law = None
    if saved is None and profile["weight"]:
        model = equilibrium.load("column")
        equilibrium.verify_inversion(model)
        law = equilibrium.law_at(profile["weight"], model)
    files, checks = [], []
    previous = None
    for i, level in enumerate(setup["numerical_levels"]):
        filename = f"{size}_{name}_{i}.json"
        p = inputs(condition, profile, level)
        basis = dict(inputs=p, access=transport, profile=profile)
        if saved is None:
            r = kernel.solve(size, dict(p, linear_partition_m3_kg=law), transport)
            verify(r, p, condition, profile)
            write_json(base / filename, dict(basis, result=r))
            print(size, name, i, round(r["wall_s"], 2), flush=True)
        else:
            if i >= len(saved["files"]) or saved["files"][i] != filename:
                raise ValueError("Frozen refinement levels changed")
            data = read_json(base / filename)
            if any(data[key] != value for key, value in basis.items()):
                raise ValueError("Frozen input changed")
            r = data["result"]
            verify(r, p, condition, profile)
        files.append(filename)
        if previous is not None:
            check = full_refinement(previous, r, p["inlet_concentration_kg_m3"])
            checks.append(check)
            if check["passed"] and check["startup"]["passed"]:
                break
        previous = r
    return dict(files=files, selected=files[-1], refinement=checks,
                numerically_qualified=checks[-1]["passed"] and checks[-1]["startup"]["passed"])


def service(seconds, flow, shaft_power):
    """Evaluate the volume and shaft work through the first 5% crossing."""
    if seconds is None:
        return None
    if not np.isfinite([seconds, flow, shaft_power]).all() or seconds < 0 or min(flow, shaft_power) <= 0:
        raise ValueError("Expected a finite service time, positive flow and positive shaft power")
    return dict(time_s=seconds, volume_m3=flow*seconds, shaft_energy_J=shaft_power*seconds)


def evaluate(base, setup, frozen):
    audit = read_json(retention.REFERENCE)
    if sha(PIXELS) != audit["source"]["pixel_sha256"]:
        raise ValueError("The Figure 7 digitization changed")
    basis = setup["basis"]
    scores = {}
    for size, (series, _) in CASES.items():
        experiment = audit["cases"][size]
        mid, end = (experiment[k] for k in ("experimental_midpoint", "last_marker"))
        reading = experiment["reading_scenarios_retention_kg"]
        flow = basis["flow_m3_s"]
        power = setup["pumps"][size]["shaft_power_W"]
        observed = {}
        predictions = {}
        for name, item in frozen["cases"][size].items():
            saved = read_json(base / item["selected"])
            eq = {k: saved["inputs"][k] for k in ("capacity_kg_kg", "affinity_m3_kg")}
            model = retention.checked_model(saved, dict(basis, **eq))
            m = compare(saved["result"], series)
            middle, retained = (float(np.interp(point["time_s"], model["time_s"], model["retained_kg"]))
                                for point in (mid, end))
            within_reading = reading[0] <= retained <= reading[1]
            predictions[name] = dict(comparison=m, retained_kg=retained, post_midpoint_retained_kg=retained-middle,
                retention_screen_passed=within_reading, numerical_screen_passed=item["numerically_qualified"],
                passed=column_passed(m) and within_reading and item["numerically_qualified"],
                service=service(m["t05"]["predicted_s"], flow, power))
            if not observed:
                timing = m["t05"]
                observed = dict(t05_s=timing["observed_interpolated_s"], t50_s=mid["time_s"],
                    retained_kg=end["figure7_retained_kg"], retention_reading_kg=reading,
                    post_midpoint_retained_kg=end["figure7_retained_kg"]-mid["figure7_retained_kg"],
                    window_end_s=end["time_s"], service=service(timing["observed_interpolated_s"], flow, power),
                    service_readings=[service(t, flow, power) for t in timing["reading_sensitivity_s"]])
        scores[size] = dict(observed=observed, predictions=predictions,
                            shaft_energy_per_volume_Wh_m3=power/flow/3600)
    qualified = {name: all(row["predictions"][name]["passed"] for row in scores.values())
                 for name in setup["profiles"]}
    return dict(cases=scores, profile_passed_both_sizes=qualified, adsorption_qualified=any(qualified.values()))


def frozen_curves(base, index):
    return {name: sha(base / name) for variants in index.values() for item in variants.values()
            for name in item["files"]}


def predict(workers):
    setup = plan()
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ") + "-" + uuid4().hex[:8]
    base = OUTPUT / "runs" / stamp
    base.mkdir(parents=True)
    write_json(base / "plan.json", setup)
    write_json(base / "gac4_batch.json", gac4_batch(setup["profiles"]["gac4"]))
    (base / "solver.jl").write_text(kernel.source())
    print(base, flush=True)

    def process(job):
        size, name = job
        args = [sys.executable, "-m", __spec__.name, str(base), "--worker", size, name]
        print(size, name, execute(args, base / f"{size}_{name}.log", timeout=3600), flush=True)

    with ThreadPoolExecutor(max_workers=workers) as pool:
        list(pool.map(process, product(setup["cases"], setup["profiles"])))
    index = {size: {name: read_json(base / f"{size}_{name}_levels.json") for name in setup["profiles"]}
             for size in setup["cases"]}
    frozen = dict(plan_sha256=sha(base / "plan.json"), cases=index, curve_sha256=frozen_curves(base, index))
    write_json(base / "frozen.json", frozen)
    print("Frozen Figure 7 forecasts", sha(base / "frozen.json"), flush=True)
    files = [Path(__file__), HERE / "protocol.md", Path(kernel.__file__), kernel.ORIGINAL, access.HERE / "cohorts.jl",
             Path(equilibrium.__file__), equilibrium.HERE / "equilibrium.jl", BATCH_REFERENCE, BLEND_REFERENCE,
             ACCESS_REFERENCE, PUMP_REFERENCE, retention.REFERENCE, PIXELS]
    write_json(base / "summary.json", dict(schema="pore2process.figure7-equilibrium.v1",
        sources=source_manifest(files), frozen_sha256=sha(base / "frozen.json"),
        solver_sha256=sha(base / "solver.jl"), gac4_batch_sha256=sha(base / "gac4_batch.json"),
        figure7_fitted_parameters=0, versions=runtime_versions(), comparison=evaluate(base, setup, frozen)))
    print(base / "summary.json", flush=True)


def check(base):
    setup, frozen, report = (read_json(base / filename) for filename in ("plan.json", "frozen.json", "summary.json"))
    if (setup != plan() or report["frozen_sha256"] != sha(base / "frozen.json")
            or frozen["plan_sha256"] != sha(base / "plan.json")
            or report["solver_sha256"] != sha(base / "solver.jl")
            or (base / "solver.jl").read_text() != kernel.source()
            or report["gac4_batch_sha256"] != sha(base / "gac4_batch.json")
            or any(sha(ROOT / file) != digest for file, digest in report["sources"].items())
            or frozen_curves(base, frozen["cases"]) != frozen["curve_sha256"]):
        raise ValueError("A recorded input or forecast changed")
    if gac4_batch(setup["profiles"]["gac4"]) != read_json(base / "gac4_batch.json"):
        raise ValueError("GAC4 batch qualification differs on replay")
    for size, variants in frozen["cases"].items():
        for name, saved in variants.items():
            if refinement(base, size, name, saved) != saved:
                raise ValueError("Refinement differs on replay")
    if evaluate(base, setup, frozen) != report["comparison"]:
        raise ValueError("The Figure 7 comparison differs on replay")
    print("Checked", len(frozen["curve_sha256"]), "frozen column curves", flush=True)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("base", type=Path, nargs="?")
    parser.add_argument("--workers", type=int, default=4)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--check", action="store_true")
    mode.add_argument("--worker", nargs=2)
    args = parser.parse_args()
    if args.worker:
        size, profile = args.worker
        write_json(args.base / f"{size}_{profile}_levels.json", refinement(args.base, size, profile))
    elif args.check:
        check(args.base)
    else:
        predict(args.workers)
