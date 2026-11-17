"""Freeze fixed-mass Julia adsorption predictions from the dense TRT packing."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from itertools import product
from pathlib import Path
import sys
from uuid import uuid4

import numpy as np

from examples.adsorption.alhamed.study import crossing
from examples.adsorption.carbon_access_transfer import kernel, run as transfer, source
from examples.adsorption.disordered_early_bypass import run as bypass
from examples.adsorption.disordered_memory_column.run import LEVELS
from examples.adsorption.dispersion_holdout.run import execute
from examples.adsorption.figure7_retention.run import checked_model
from examples.adsorption.porosity_matched_film import run as dense
from examples.adsorption.shared_contact.run import full_refinement
from examples.adsorption.transport_accuracy.run import validate_record
from pore2process_demo.bridge import runtime_versions
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT
from pore2process_demo.provenance import source_manifest

HERE = Path(__file__).resolve().parent
OUTPUT = ROOT / "results/adsorption/porosity_matched_column"
FLOWS = bypass.FLOWS
GRIDS = (45, 60)
VARIANTS = {f"film{film}_path{path}": (film, path) for film, path in product(GRIDS, repeat=2)}
PRIMARY = "film60_path60"
BASELINE = "archived_correlation"


def parent_inputs():
    provenance = read_json(dense.HERE / "reference_summary.json")["run"]
    path = ROOT / provenance["path"]
    if sha(path) != provenance["sha256"] or sha(bypass.COLUMN) != bypass.COLUMN_SHA:
        raise ValueError("The fixed packing or baseline column record changed")
    native = validate_record(path)
    old = validate_record(bypass.COLUMN)
    if native["figure7_fits"] or old["figure_parameters_fitted"]:
        raise ValueError("Expected unfitted parent records")
    prior = bypass.COLUMN.parent
    manifest = prior / "frozen.json"
    if sha(manifest) != old["frozen_sha256"]:
        raise ValueError("Previous column freeze changed")
    return path, native["analysis"]["levels"], read_json(prior / "plan.json"), read_json(manifest)


def plan():
    path, levels, old_plan, old_frozen = parent_inputs()
    cases = {}
    for flow in FLOWS:
        case = old_plan["cases"][flow + "_l25"]
        original = case["parameters"]
        access = case["scenarios"]["f50_slow01"]
        direct, resistance = access["inputs"]["direct_fraction"], access["added_resistance_s_m"]
        base = old_frozen["cases"][flow + "_l25"]["f50_slow01"]["selected"]
        old_curve = bypass.COLUMN.parent / base
        baseline = dict(path=old_curve.relative_to(ROOT).as_posix(), sha256=sha(old_curve))
        if baseline["sha256"] != old_frozen["curves_sha256"][base]:
            raise ValueError("Archived baseline curve changed")
        variants = {}
        for key, (film, rules) in VARIANTS.items():
            estimate = levels[f"N{film}_20us"][flow]
            geometry = levels[f"N{rules}_20us"][flow]["column"]
            for row in (estimate["column"], geometry):
                values = (row["dry_carbon_mass_kg"], row["liquid_volume_m3"], row["flow_m3_s"],
                          row["direct_fraction"], row["added_resistance_s_m"], row["Pe"], row["sections"])
                expected = (original["dry_mass_kg"], original["liquid_volume_m3"], original["flow_m3_s"],
                            direct, resistance, case["Pe"], case["sections"])
                np.testing.assert_allclose(values, expected, rtol=1e-13)
            coefficient = estimate["film_coefficient_m_s"][-1]
            variants[key] = dict(film_grid=film, transport_grid=rules, film_m_s=coefficient,
                film_matched_clean_bed_conductance_m3_s=estimate["clean_film"]["carbon_area"]
                    ["effective_conductance_m3_s"],
                probe_directory=(path.parent / f"tracer_N{rules}_20us").relative_to(ROOT).as_posix(),
                access=dict(direct_fraction=direct, slow_ratio=1 / (1 + resistance * coefficient)))
        cases[flow] = dict(parameters=original, flow_mL_min=case["flow_mL_min"],
            sections=case["sections"], Pe=case["Pe"], added_resistance_s_m=resistance,
            variants=variants, baseline=baseline)
    return dict(cases=cases, primary=PRIMARY, dense_summary=path.relative_to(ROOT).as_posix(),
                dense_sha256=sha(path), prior_sha256=bypass.COLUMN_SHA)


def solve_variant(base, flow, name):
    setup = read_json(base / "plan.json")
    condition = setup["cases"][flow]
    spec = condition["variants"][name]
    feed = condition["parameters"]["inlet_concentration_kg_m3"]
    files, checks = [], []
    previous = None
    for level, (shells, dt, early) in enumerate(LEVELS):
        p = dict(condition["parameters"], film_m_s=spec["film_m_s"], radial_cells=shells,
                 dt_s=dt, startup_dt_s=early)
        s = spec["access"]
        result = kernel.solve(condition["flow_mL_min"], ROOT / spec["probe_directory"], p, s)
        transfer.verify_run(result, p, s, condition)
        np.testing.assert_allclose(result["effective_conductance_m3_s"],
            spec["film_matched_clean_bed_conductance_m3_s"], rtol=1e-12)
        filename = f"{flow}_{name}_{level}.json"
        write_json(base / filename, dict(inputs=p, scenario=s, result=result))
        files.append(filename)
        print(flow, name, level, "wall seconds", round(result["wall_s"], 1), flush=True)
        if previous is not None:
            check = full_refinement(previous, result, feed)
            checks.append(check)
            if check["passed"] and check["startup"]["passed"]:
                break
        previous = result
    write_json(base / f"{flow}_{name}_levels.json", dict(files=files, selected=files[-1],
        refinement=checks, numerically_qualified=checks[-1]["passed"] and checks[-1]["startup"]["passed"]))


def predict():
    setup = plan()
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ") + "-" + uuid4().hex[:8]
    base = OUTPUT / "runs" / stamp
    base.mkdir(parents=True)
    write_json(base / "plan.json", setup)
    (base / "solver.jl").write_text(kernel.source())
    print("Frozen fixed-mass inputs", base, sha(base / "plan.json"), flush=True)

    def process(job):
        flow, name = job
        cmd = [sys.executable, "-m", __spec__.name, str(base), "--variant", flow, name]
        timing = execute(cmd, base / f"{flow}_{name}.log", timeout=3600)
        print(flow, name, "refinement complete", timing, flush=True)

    with ThreadPoolExecutor(max_workers=3) as pool:
        list(pool.map(process, product(FLOWS, VARIANTS)))
    cases = {flow: {name: read_json(base / f"{flow}_{name}_levels.json") for name in VARIANTS}
             for flow in FLOWS}
    file = base / "frozen.json"
    write_json(file, dict(plan_sha256=sha(base / "plan.json"), solver_sha256=sha(base / "solver.jl"),
        cases=cases, curve_sha256={f: sha(base / f) for variants in cases.values()
                                  for item in variants.values() for f in item["files"]}))
    print("Predictions frozen before reading experimental curves", sha(file), flush=True)
    print(base, flush=True)


def archive_curve(spec):
    path = ROOT / spec["path"]
    if sha(path) != spec["sha256"]:
        raise ValueError("Archived comparison curve changed")
    return read_json(path)


def frozen_inputs(base):
    setup = read_json(base / "plan.json")
    frozen = read_json(base / "frozen.json")
    if (setup != plan() or frozen["plan_sha256"] != sha(base / "plan.json")
            or frozen["solver_sha256"] != sha(base / "solver.jl")
            or (base / "solver.jl").read_text() != kernel.source()):
        raise ValueError("Frozen column setup or solver source changed")
    if set(frozen["cases"]) != set(FLOWS):
        raise ValueError("Frozen prediction is missing a flow")
    for flow, variants in frozen["cases"].items():
        condition = setup["cases"][flow]
        if set(variants) != set(VARIANTS):
            raise ValueError("Frozen prediction is missing a grid comparison")
        for name, item in variants.items():
            spec = condition["variants"][name]
            previous = None
            checks = []
            for filename, (shells, dt, early) in zip(item["files"], LEVELS):
                path = base / filename
                if sha(path) != frozen["curve_sha256"][filename]:
                    raise ValueError("Frozen column curve changed")
                data = read_json(path)
                p = dict(condition["parameters"], film_m_s=spec["film_m_s"], radial_cells=shells,
                         dt_s=dt, startup_dt_s=early)
                if data["inputs"] != p or data["scenario"] != spec["access"]:
                    raise ValueError("A fixed carbon, particle or solver input changed")
                result = data["result"]
                transfer.verify_run(result, p, spec["access"], condition)
                np.testing.assert_allclose(result["effective_conductance_m3_s"],
                    spec["film_matched_clean_bed_conductance_m3_s"], rtol=1e-12)
                if previous is not None:
                    checks.append(full_refinement(previous, result, p["inlet_concentration_kg_m3"]))
                previous = result
            if not checks or len(item["files"]) != len(checks) + 1 or item["selected"] != item["files"][-1]:
                raise ValueError("Frozen prediction has incomplete numerical refinement")
            if (checks != item["refinement"] or item["numerically_qualified"] != (
                    checks[-1]["passed"] and checks[-1]["startup"]["passed"])):
                raise ValueError("Saved numerical refinement differs on replay")
    return setup, frozen


def curve_metrics(saved, basis, data, obs, common_time):
    model = checked_model(saved, basis)
    time, fraction = model["time_s"], model["outlet_fraction"]
    mass = basis["dry_mass_kg"]
    sampled = np.interp(data["time"], time, fraction)
    predictions = {f"t{i:02d}_s": crossing(time, fraction, i / 100) for i in (5, 50)}
    timing = transfer.crossing_metrics(predictions, obs)
    rmse = float(np.sqrt(np.mean((sampled - data["fraction"])**2)))
    stop = obs["end_s"]
    retained = float(np.interp(stop, time, model["retained_kg"])) / mass
    low, high = obs["reading_retention_kg_kg"]
    slow = np.asarray(saved["result"]["adsorbed_inventory_by_cohort_kg"][1])
    return dict(**timing, rmse=rmse, retained_kg_kg=retained, retained_in_reading_range=low <= retained <= high,
        common_dose_retained_kg_kg=float(np.interp(common_time, time, model["retained_kg"])) / mass,
        slow_adsorbed_kg_kg=float(np.interp(stop, time, slow)) / mass,
        outlet_after_three_liquid_residence_times=float(np.interp(
            3 * basis["liquid_volume_m3"] / basis["flow_m3_s"], time, fraction)),
        curve_screen_passed=rmse <= .1 and timing["t50_s"]["passed"] is True
            and timing["t05_s"]["passed"] is not False,
        normalized_mass_residual=saved["result"]["normalized_mass_residual"],
        normalized_interface_flux_mismatch=saved["result"]["normalized_interface_flux_mismatch"])


def sensitivity(scores):
    """Signed N45 minus N60 effects at each fixed film or transport input."""
    metrics = ("t05_s", "t50_s", "retained_kg_kg", "common_dose_retained_kg_kg", "rmse")
    values = {m: {k: (scores[k][m]["predicted_s"] if m in metrics[:2] else scores[k][m])
                  for k in VARIANTS} for m in metrics}

    def delta(a, b):
        return None if a is None or b is None else a - b

    result = {}
    for metric, v in values.items():
        complete = None not in v.values()
        result[metric] = dict(
            paired_N45_minus_N60=delta(v["film45_path45"], v["film60_path60"]),
            film_N45_minus_N60_fixed_transport={str(g): delta(v[f"film45_path{g}"], v[f"film60_path{g}"])
                                               for g in GRIDS},
            transport_N45_minus_N60_fixed_film={str(g): delta(v[f"film{g}_path45"], v[f"film{g}_path60"])
                                               for g in GRIDS},
            range_across_four_variants=[min(v.values()), max(v.values())] if complete else None)
    return result


def compare(base, setup, frozen):
    markers = source.markers()
    keys = {f"f4_{flow}_l25": flow for flow in FLOWS}
    parameters = {flow: case["parameters"] for flow, case in setup["cases"].items()}
    feeds = {flow: p["flow_m3_s"] * p["inlet_concentration_kg_m3"] for flow, p in parameters.items()}
    common_dose = min(float(markers[key]["time"][-1]) * feeds[flow] / parameters[flow]["dry_mass_kg"]
                      for key, flow in keys.items())
    keys["f5_q233_l25"] = "q233"
    comparisons = {}
    for key, flow in keys.items():
        case = setup["cases"][flow]
        p, feed = parameters[flow], feeds[flow]
        common_time = min(common_dose * p["dry_mass_kg"] / feed, markers[key]["time"][-1])
        obs = source.observation(markers[key], feed, p["dry_mass_kg"], common_time)
        curves = {name: read_json(base / item["selected"]) for name, item in frozen["cases"][flow].items()}
        curves[BASELINE] = archive_curve(case["baseline"])
        scores = {name: curve_metrics(saved, dict(p, feed_rate_kg_s=feed), markers[key], obs, common_time)
                  for name, saved in curves.items()}
        comparisons[key] = dict(flow=flow, figure=markers[key]["figure"], observed=obs,
                                variants=scores, grid_sensitivity=sensitivity(scores))
    return dict(source=read_json(source.EXTRACTION), common_dose_kg_kg=common_dose, comparisons=comparisons)


def score(base):
    path = base / "summary.json"
    if path.exists():
        raise ValueError("Use --check to verify a sealed comparison")
    setup, frozen = frozen_inputs(base)
    observed = compare(base, setup, frozen)
    files = [Path(__file__), HERE / "protocol.md", HERE / "plot.py", transfer.HERE / "kernel.py",
             transfer.HERE / "run.py", transfer.HERE / "source.py", dense.HERE / "reference_summary.json",
             kernel.access.ORIGINAL, dense.HERE / "run.py"]
    write_json(path, dict(schema="pore2process.porosity-matched-column.v1", run_id=base.name,
        sources=dict(experimental_pixels=sha(source.PIXELS), extraction=sha(source.EXTRACTION)),
        frozen_sha256=sha(base / "frozen.json"), source_sha256=source_manifest(files),
        artifacts_sha256={p.name: sha(p) for p in sorted(base.iterdir()) if p.is_file()},
        versions=runtime_versions(), figure_parameters_fitted=0, score=observed))
    print(path, flush=True)


def check(base):
    record = validate_record(base / "summary.json")
    if (record["frozen_sha256"] != sha(base / "frozen.json")
            or record["sources"] != dict(experimental_pixels=sha(source.PIXELS), extraction=sha(source.EXTRACTION))):
        raise ValueError("Sealed predictions or observation readout changed")
    if record["score"] != compare(base, *frozen_inputs(base)):
        raise ValueError("Archived scores differ on replay")
    print(dict(sources=len(record["source_sha256"]), artifacts=len(record["artifacts_sha256"])))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("base", type=Path, nargs="?")
    mode = p.add_mutually_exclusive_group()
    mode.add_argument("--variant", nargs=2)
    mode.add_argument("--score", action="store_true")
    mode.add_argument("--check", action="store_true")
    args = p.parse_args()
    if (args.variant or args.score or args.check) and args.base is None:
        p.error("Pass the column run directory")
    if args.variant:
        solve_variant(args.base.resolve(), *args.variant)
    elif args.score:
        score(args.base.resolve())
    elif args.check:
        check(args.base.resolve())
    else:
        predict()


if __name__ == "__main__":
    main()
