"""Test confined flow and distinct carbon histories at fixed total carbon."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from functools import lru_cache
from pathlib import Path
import sys
from uuid import uuid4

import numpy as np

from examples.adsorption.breakthrough_hypotheses.batch.conditional import run as equilibrium
from examples.adsorption.figure7_retention.run import integrate
from examples.adsorption.porosity_matched_column import run as column
from pore2process_demo.bridge import runtime, runtime_versions
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT

from .screen import HERE, OUTPUT, PARENT

CASES = ("q233_l10", "q233_l15", "q233_l20", "q233_l25", "q625_l25")
PERIODIC = dict(flow=[1.], liquid=[1.], carbon=[1.], mixing=0.)
HOLDOUTS = ("q233_l20", "q233_l25", "q625_l25")


def exposure_groups(flow, contact):
    """Use the existing clean tracer's contact per unit flow to group paths."""
    order = np.argsort(contact.mean(axis=1) / flow, kind="stable")
    split = np.argmin(abs(np.cumsum(flow[order])[:-1] - .5)) + 1
    labels = np.full(len(flow), 2)
    labels[order[:split]] = 1
    return labels


def exchange_per_section(transition, wall_flow, sample_length_m, section_length_m):
    """Project observed crossovers onto stationary, flux-conserving exchange."""
    transition = np.asarray(transition, dtype=float)
    if transition.shape != (2, 2) or np.any(transition < 0) or not np.isfinite(transition).all():
        raise ValueError("Expected a two-region forward-tracer transition matrix")
    np.testing.assert_allclose(transition.sum(axis=1), 1., rtol=0., atol=1e-12)
    if not 0 < wall_flow < 1 or not min(sample_length_m, section_length_m) > 0:
        raise ValueError("Invalid flow share or axial length")
    core_flow = 1 - wall_flow
    directed_flow = [wall_flow * transition[0, 1], core_flow * transition[1, 0]]
    attenuation = sum(directed_flow) / (2 * wall_flow * core_flow)
    # A first-order exchange law can attenuate a concentration difference to zero.
    # Retain the raw value in the record if the projection reaches that limit.
    mixing = -np.expm1(np.log1p(-attenuation) * section_length_m / sample_length_m) if attenuation < 1 else 1.
    return dict(mixing=float(mixing), directed_crossover_flow_fractions=directed_flow,
                projected_release_over_sample=float(attenuation))


def fixed_inputs():
    record = read_json(column.HERE / "reference_summary.json")
    parent_plan, parent_frozen = (PARENT.parent / name for name in ("plan.json", "frozen.json"))
    for path, name in ((PARENT, "run"), (parent_plan, "column_inputs"), (parent_frozen, "frozen_predictions")):
        if sha(path) != record[name]["sha256"]:
            raise ValueError(f"Archived N60 {name} changed")
    setup, frozen = read_json(parent_plan), read_json(parent_frozen)
    _, _, old_plan, _ = column.parent_inputs()
    cases = {}
    for name in CASES:
        flow, depth = name.split("_l")
        prior = setup["cases"][flow]
        spec = prior["variants"][column.PRIMARY]
        condition = old_plan["cases"][name]
        archived = None
        if depth == "25":
            archive = frozen["cases"][flow][column.PRIMARY]["selected"]
            archived = dict(path=(PARENT.parent / archive).relative_to(ROOT).as_posix(),
                            sha256=frozen["curve_sha256"][archive])
            if prior["parameters"] != condition["parameters"]:
                raise ValueError("N60 reference is on a different bed basis")
        cases[name] = dict(flow=flow, depth_cm=int(depth), sections=condition["sections"],
            flow_mL_min=condition["flow_mL_min"],
            parameters=dict(condition["parameters"], film_m_s=spec["film_m_s"]),
            probe_directory=spec["probe_directory"], access=spec["access"],
            conductance_25cm_m3_s=spec["film_matched_clean_bed_conductance_m3_s"], archive=archived)
    return cases


def region_settings(source, section_length_m):
    settings = {}
    for key, value in source.items():
        split = {name: [value[name], 1 - value[name]] for name in ("flow", "liquid", "carbon")}
        if any(min(f) <= 0 for f in split.values()):
            raise ValueError("Expected positive flows and inventories in both regions")
        exchange = exchange_per_section(value["transition"], value["flow"], value["length_m"], section_length_m)
        settings[key] = dict(split, **exchange)
    return settings


def verify_handoff(data):
    sources = [data.get("sources", {}), *(case.get("sources", {}) for case in data["cases"].values())]
    if any(sha(ROOT / ref["path"]) != ref["sha256"] for records in sources for ref in records.values()):
        raise ValueError("An input to the confined flow or tracer handoff changed")


def calibration(base):
    setup = read_json(base / "plan.json")
    chosen = read_json(base / "selection.json")
    comparison = equilibrium.compare(base, setup["qualified_profiles"], "calibration")
    if chosen != equilibrium.select(comparison):
        raise ValueError("The shared equilibrium weight differs from the frozen short-bed calibration")
    return dict(profile=setup["qualified_profiles"][chosen["candidate"]],
                sources={str((base / file).resolve().relative_to(ROOT)): sha(base / file)
                         for file in ("plan.json", "selection.json", "frozen_calibration.json")})


def plan(wall_file=None, calibration_base=None):
    wall_data = read_json(wall_file) if wall_file is not None else None
    cases = fixed_inputs()
    if wall_data is not None:
        verify_handoff(wall_data)
        cases = {name: case for name, case in cases.items() if case["flow"] in wall_data["cases"]}
        if not cases:
            raise ValueError("No matching flow rate in the confined-tracer handoff")
    choice = None
    if calibration_base is not None:
        if wall_data is None or not set(HOLDOUTS) <= cases.keys():
            raise ValueError("Joint comparison needs the calibrated profile and both confined flow rates")
        cases = {name: cases[name] for name in HOLDOUTS}
        choice = calibration(calibration_base)
    for case in cases.values():
        if wall_data is None:
            case["variants"] = dict(baseline=[PERIODIC, "shared"], histories=[PERIODIC, "exposure"])
        else:
            section_length = .01 * case["depth_cm"] / case["sections"]
            split = region_settings(wall_data["cases"], section_length)[case["flow"]]
            case["variants"] = dict(wall_flow=[split, "shared"], wall_and_histories=[split, "radial"])
            if choice is not None:
                profile = choice["profile"]
                case["parameters"]["diffusivity_m2_s"] = profile["diffusivity_m2_s"]
                case["freundlich_weight"] = profile["weight"]
                case["variants"] = dict(joint=[split, "radial"])
    return dict(cases=cases, wall_source=wall_data, parent_sha256=sha(PARENT),
                code_sha256=sha(HERE / "model.jl"), figure_parameters_fitted=int(choice is not None),
                **({"selection": choice} if choice is not None else {}))


@lru_cache(maxsize=1)
def solver(file):
    column.kernel.solver()
    jl = runtime()
    jl.include(str(file))
    return jl.BreakthroughHypotheses.solve


def simulate(file, case, variant, p):
    split, rule = case["variants"][variant]
    paths = column.kernel.transport(case["flow_mL_min"], ROOT / case["probe_directory"])
    flow = np.array(paths.flow)
    contact = np.array(paths.area)
    zones = len(split["flow"])
    groups = np.ones(flow.size * zones, dtype=int)
    if rule == "exposure":
        groups = exposure_groups(flow, contact)
    elif rule == "radial":
        groups[flow.size:] = 2
    solve = solver(file)
    particle = p
    if "freundlich_weight" in case:
        model = equilibrium.load("column")
        equilibrium.verify_inversion(model)
        law = equilibrium.law_at(case["freundlich_weight"], model)
        particle = dict(p, linear_partition_m3_kg=law)
    value = solve(paths.sections, paths.flow, paths.area,
        *[np.array(split[k]) for k in ("flow", "liquid", "carbon")], groups,
        mixing=split["mixing"], **case["access"], **particle)
    jl = runtime()
    result = {str(f): np.array(v, dtype=float).tolist() if bool(jl.isa(v, jl.AbstractArray)) else float(v)
              for f in jl.fieldnames(jl.typeof(value)) for v in (getattr(value, str(f)),)}
    verify(result, p, case)
    return result, groups.tolist()


def verify(result, p, case):
    conductance = case["conductance_25cm_m3_s"] * case["depth_cm"] / 25
    if "freundlich_weight" in case:
        spec = dict(access=case["access"], film_matched_clean_bed_conductance_m3_s=conductance)
        equilibrium.verify_column(result, p, spec, case, case["freundlich_weight"])
    else:
        column.transfer.verify_run(result, p, case["access"], case)
        np.testing.assert_allclose(result["effective_conductance_m3_s"], conductance, rtol=1e-12)


def refinement(base, name, variant, saved=None):
    case = read_json(base / "plan.json")["cases"][name]
    feed = case["parameters"]["inlet_concentration_kg_m3"]
    files, checks = [], []
    previous = None
    for level, (shells, dt, early) in enumerate(column.LEVELS):
        filename = f"{name}_{variant}_{level}.json"
        p = dict(case["parameters"], radial_cells=shells, dt_s=dt, startup_dt_s=early)
        if saved is None:
            result, groups = simulate(base / "model.jl", case, variant, p)
            write_json(base / filename, dict(inputs=p, scenario=case["access"], result=result, groups=groups))
            print(name, variant, level, round(result["wall_s"], 2), flush=True)
        else:
            if len(saved["files"]) <= level or saved["files"][level] != filename:
                raise ValueError("The frozen run is missing an expected refinement level")
            data = read_json(base / filename)
            if data["inputs"] != p or data["scenario"] != case["access"]:
                raise ValueError("A frozen bed/solver input changed")
            result = data["result"]
            verify(result, p, case)
        files.append(filename)
        if previous is not None:
            check = column.full_refinement(previous, result, feed)
            checks.append(check)
            if check["passed"] and check["startup"]["passed"]:
                break
        previous = result
    return dict(files=files, selected=files[-1], checks=checks,
                qualified=checks[-1]["passed"] and checks[-1]["startup"]["passed"])


def archived_parity(saved, case):
    archived = column.archive_curve(case["archive"])
    result, expected = saved["result"], archived["result"]
    np.testing.assert_array_equal(result["time_s"], expected["time_s"])
    for key in ("outlet_concentration_kg_m3", "liquid_inventory_kg", "adsorbed_inventory_kg"):
        np.testing.assert_allclose(result[key], expected[key], rtol=1e-9, atol=1e-13)
    return True


def metrics(saved, p, observed, obs):
    feed = p["inlet_concentration_kg_m3"] * p["flow_m3_s"]
    mass = p["dry_mass_kg"]
    basis = dict(p, feed_rate_kg_s=feed)
    score = column.curve_metrics(saved, basis, observed, obs, obs["end_s"])
    model = column.checked_model(saved, basis)
    mid = obs["t50_s"]
    retained = float(np.interp(mid, model["time_s"], model["retained_kg"])) / mass
    published = integrate(np.r_[0., observed["time"]], np.r_[0., observed["fraction"]], feed, mid) / mass
    score["post_midpoint_retained_kg_kg"] = score["retained_kg_kg"] - retained
    score["observed_post_midpoint_retained_kg_kg"] = obs["retained_kg_kg"] - published
    early, middle = (score[k]["predicted_s"] for k in ("t05_s", "t50_s"))
    score["rise_s"] = middle - early if early is not None and middle is not None else None
    return score


def score(base, setup, frozen):
    markers = column.source.markers()
    comparisons = {}
    for name, case in setup["cases"].items():
        figure = 4 if case["depth_cm"] == 25 else 5
        observed = markers[f"f{figure}_{name}"]
        p = case["parameters"]
        obs = column.source.observation(observed, p["inlet_concentration_kg_m3"] * p["flow_m3_s"],
                                        p["dry_mass_kg"], observed["time"][-1])
        curves, parity = {}, None
        for variant, level in frozen["cases"][name].items():
            saved = read_json(base / level["selected"])
            if variant == "baseline" and case["archive"] is not None:
                parity = archived_parity(saved, case)
            curves[variant] = dict(metrics(saved, p, observed, obs), qualified=level["qualified"],
                                   groups=saved["groups"])
        comparisons[name] = dict(figure=figure, depth_cm=case["depth_cm"], dry_mass_kg=p["dry_mass_kg"],
                                 observed=obs, curves=curves, archived_parity=parity)
    transfer_screen = {variant: all(
        name in comparisons and comparisons[name]["curves"][variant]["qualified"] and
        comparisons[name]["curves"][variant]["curve_screen_passed"] and
        comparisons[name]["curves"][variant]["retained_in_reading_range"] for name in HOLDOUTS)
        for variant in next(iter(setup["cases"].values()))["variants"]}
    return dict(comparisons=comparisons, all_holdouts_passed=transfer_screen)


def hashes(base, index):
    return {file: sha(base / file) for variants in index.values() for item in variants.values() for file in item["files"]}


def predict(wall, workers, calibration_base=None):
    setup = plan(wall, calibration_base)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ") + "-" + uuid4().hex[:8]
    base = OUTPUT / "column" / stamp
    base.mkdir(parents=True)
    write_json(base / "plan.json", setup)
    (base / "model.jl").write_bytes((HERE / "model.jl").read_bytes())
    cases = setup["cases"]
    print(base, flush=True)

    def process(job):
        case, variant = job
        cmd = [sys.executable, "-m", __spec__.name, str(base), "--worker", case, variant]
        print(case, variant, column.execute(cmd, base / f"{case}_{variant}.log", timeout=3600), flush=True)

    jobs = [(case, variant) for case in reversed(cases) for variant in cases[case]["variants"]]
    with ThreadPoolExecutor(max_workers=workers) as pool:
        list(pool.map(process, jobs))
    index = {case: {v: read_json(base / f"{case}_{v}_levels.json") for v in c["variants"]} for case, c in cases.items()}
    frozen = dict(cases=index, curves=hashes(base, index), plan_sha256=sha(base / "plan.json"))
    write_json(base / "frozen.json", frozen)
    print("Frozen forecasts", sha(base / "frozen.json"), flush=True)
    sources = (Path(__file__), HERE / "model.jl", HERE / "protocol.md")
    if calibration_base is not None:
        sources += (HERE / "coupled_protocol.md", Path(equilibrium.__file__), equilibrium.HERE / "equilibrium.jl")
    digest = {p.resolve().relative_to(ROOT).as_posix(): sha(p) for p in sources}
    write_json(base / "summary.json", dict(frozen_sha256=sha(base / "frozen.json"), sources=digest,
        figure_sha256=sha(column.source.PIXELS), versions=runtime_versions(),
        figure_parameters_fitted=setup["figure_parameters_fitted"], figure7_parameters_fitted=0,
        score=score(base, setup, frozen)))
    print(base / "summary.json", flush=True)


def check(base):
    setup, frozen, report = (read_json(base / file) for file in ("plan.json", "frozen.json", "summary.json"))
    if (report["frozen_sha256"] != sha(base / "frozen.json") or frozen["plan_sha256"] != sha(base / "plan.json")
            or setup["code_sha256"] != sha(base / "model.jl") or setup["parent_sha256"] != sha(PARENT)
            or hashes(base, frozen["cases"]) != frozen["curves"]
            or report["figure_sha256"] != sha(column.source.PIXELS)):
        raise ValueError("An input or curve has changed since the forecast was frozen")
    if "selection" in setup:
        selection = setup["selection"]
        if any(sha(ROOT / file) != digest for file, digest in selection["sources"].items()):
            raise ValueError("The upstream short-bed calibration changed")
    if setup["wall_source"] is not None:
        verify_handoff(setup["wall_source"])
    for name, variants in frozen["cases"].items():
        for variant, saved in variants.items():
            if refinement(base, name, variant, saved) != saved:
                raise ValueError("Numerical refinement differs on replay")
    if report["score"] != score(base, setup, frozen):
        raise ValueError("Published comparison differs on replay")
    print("Checked", base, len(frozen["curves"]), "curves", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("base", type=Path, nargs="?")
    parser.add_argument("--wall", type=Path, help="Frozen confined-flow and tracer handoff")
    parser.add_argument("--calibration", type=Path, help="Frozen two-depth equilibrium calibration for joint checks")
    parser.add_argument("--workers", type=int, default=2)
    action = parser.add_mutually_exclusive_group()
    action.add_argument("--check", action="store_true")
    action.add_argument("--worker", nargs=2)
    args = parser.parse_args()
    if args.worker:
        case, variant = args.worker
        write_json(args.base / f"{case}_{variant}_levels.json", refinement(args.base, case, variant))
    elif args.check:
        check(args.base)
    else:
        predict(args.wall, args.workers, args.calibration)
