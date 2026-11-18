"""Bound the effect of a faster annulus with zero exchange to the core."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from math import pi
from pathlib import Path
import sys
from uuid import uuid4

import numpy as np

from examples.adsorption.porosity_matched_column import run as column
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT

HERE = Path(__file__).resolve().parent
OUTPUT = ROOT / "results/adsorption/breakthrough_hypotheses"
PARENT = ROOT / read_json(column.HERE / "reference_summary.json")["run"]["path"]
WALL_VOIDAGES = (.47, .50, .55)


def annulus(flow, diameter, particle, voidage, wall_voidage, viscosity, density):
    """Conserve total flow and inventory with two parallel Ergun beds."""
    if not 0 < particle < diameter or not 0 < voidage < 1 or not 0 < wall_voidage < 1:
        raise ValueError("Invalid annulus geometry or porosity")
    area = pi * diameter**2 / 4
    wall_area = 1 - (1 - particle / diameter)**2
    sections = np.array([wall_area, 1 - wall_area])
    core_voidage = (voidage - wall_area * wall_voidage) / (1 - wall_area)
    epsilon = np.array([wall_voidage, core_voidage])
    if not 0 < core_voidage < 1 or min(flow, viscosity, density) <= 0:
        raise ValueError("Invalid core porosity or hydraulic input")
    linear = 150 * viscosity * (1 - epsilon)**2 / (epsilon**3 * particle**2)
    quadratic = 1.75 * density * (1 - epsilon) / (epsilon**3 * particle)

    def speed(gradient):
        return 2 * gradient / (linear + np.sqrt(linear**2 + 4 * gradient * quadratic))

    mean_velocity = flow / area
    lower, upper = 0., max(linear * mean_velocity + quadratic * mean_velocity**2)
    for _ in range(60):
        gradient = (lower + upper) / 2
        if sections @ speed(gradient) < mean_velocity:
            lower = gradient
        else:
            upper = gradient
    gradient = (lower + upper) / 2
    fractions = sections * speed(gradient) / mean_velocity
    liquid = sections * epsilon / voidage
    carbon = sections * (1 - epsilon) / (1 - voidage)
    for weights in (fractions, liquid, carbon):
        np.testing.assert_allclose(weights.sum(), 1., rtol=1e-12)
    return dict(width_m=particle / 2, area=sections.tolist(), voidage=epsilon.tolist(),
                flow=fractions.tolist(), liquid=liquid.tolist(), carbon=carbon.tolist(),
                gradient_Pa_m=gradient)


def combined(results, splits, p):
    """Mix branch outlet streams on total flow and sum phase inventories."""
    a, b = results
    np.testing.assert_array_equal(a["time_s"], b["time_s"])
    values = {"time_s": a["time_s"]}
    for key in ("cumulative_outlet_kg", "liquid_inventory_kg", "adsorbed_inventory_kg",
                "adsorbed_inventory_by_cohort_kg"):
        values[key] = (np.asarray(a[key]) + np.asarray(b[key])).tolist()
    q = np.asarray(splits["flow"])
    values["outlet_concentration_kg_m3"] = (q @
        np.asarray([r["outlet_concentration_kg_m3"] for r in results])).tolist()
    t = np.asarray(values["time_s"])
    balance = (np.asarray(values["liquid_inventory_kg"]) + values["adsorbed_inventory_kg"] +
               values["cumulative_outlet_kg"] - p["flow_m3_s"] * p["inlet_concentration_kg_m3"] * t)
    values["normalized_mass_residual"] = float(max(abs(balance)) /
        (p["flow_m3_s"] * p["inlet_concentration_kg_m3"] * t[-1]))
    values["normalized_interface_flux_mismatch"] = max(
        r["normalized_interface_flux_mismatch"] for r in results)
    return values


def setup():
    parents = read_json(column.HERE / "reference_summary.json")
    if sha(PARENT) != parents["run"]["sha256"]:
        raise ValueError("Archived periodic case changed")
    case = read_json(PARENT.parent / "plan.json")["cases"]["q625"]
    spec = case["variants"][column.PRIMARY]
    p = dict(case["parameters"], film_m_s=spec["film_m_s"])
    parts = {f"voidage{value:.2f}": annulus(p["flow_m3_s"], .011, p["particle_diameter_m"],
              .427, value, .0008333333333333334, 1000.) for value in (.427, *WALL_VOIDAGES)}
    return dict(input=p, access=spec["access"], splits=parts, sections=case["sections"],
                probe_directory=spec["probe_directory"], parent_sha256=sha(PARENT))


def solve(base, name):
    plan = read_json(base / "plan.json")
    p, access = plan["input"], plan["access"]
    split = plan["splits"][name]
    feed = p["inlet_concentration_kg_m3"]
    previous = None
    for level, (shells, step, early_step) in enumerate(column.LEVELS):
        physical = dict(p, radial_cells=shells, dt_s=step, startup_dt_s=early_step)
        results = []
        for i in range(2):
            q = dict(physical, flow_m3_s=p["flow_m3_s"] * split["flow"][i],
                     dry_mass_kg=p["dry_mass_kg"] * split["carbon"][i],
                     liquid_volume_m3=p["liquid_volume_m3"] * split["liquid"][i])
            r = column.kernel.solve(62.5, ROOT / plan["probe_directory"], q, access)
            column.transfer.verify_run(r, q, access, plan)
            results.append(r)
        result = combined(results, split, physical)
        if result["normalized_mass_residual"] > 1e-9:
            raise ValueError("Combined annulus and core do not conserve phenol")
        saved = dict(inputs=physical, scenario=access, result=result)
        if previous is not None:
            checks = column.full_refinement(previous, result, feed)
            if checks["passed"] and checks["startup"]["passed"]:
                write_json(base / f"{name}.json", dict(solution=saved, refinement=checks, level=level))
                return
        previous = result
    raise ValueError("Parallel-column case did not refine")


def score(base, plan, names):
    old = read_json(PARENT)["score"]
    comparison = old["comparisons"]["f4_q625_l25"]
    data = column.source.markers()["f4_q625_l25"]
    p = plan["input"]
    feed = p["flow_m3_s"] * p["inlet_concentration_kg_m3"]
    cum = min(comparison["observed"]["end_s"], old["common_dose_kg_kg"] * p["dry_mass_kg"] / feed)
    scores = {name: column.curve_metrics(read_json(base / f"{name}.json")["solution"],
              dict(p, feed_rate_kg_s=feed), data, comparison["observed"], cum) for name in names}
    prior = comparison["variants"][column.PRIMARY]
    homogeneous = scores["voidage0.43"]
    for marker in ("t05_s", "t50_s"):
        np.testing.assert_allclose(homogeneous[marker]["predicted_s"], prior[marker]["predicted_s"], rtol=1e-10)
    np.testing.assert_allclose(homogeneous["retained_kg_kg"], prior["retained_kg_kg"], rtol=1e-10)
    return dict(observed=comparison["observed"], periodic=prior, annulus=scores)


def run():
    plan = setup()
    base = OUTPUT / "screen" / (datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ") + "-" + uuid4().hex[:8])
    base.mkdir(parents=True)
    write_json(base / "plan.json", plan)
    print(base, flush=True)

    def process(name):
        command = [sys.executable, "-m", __spec__.name, str(base), "--case", name]
        wall = column.execute(command, base / f"{name}.log", timeout=1800)
        print(name, wall, flush=True)

    with ThreadPoolExecutor(max_workers=2) as pool:
        list(pool.map(process, plan["splits"]))
    freeze = dict(plan=sha(base / "plan.json"), curves={name: sha(base / f"{name}.json")
                  for name in plan["splits"]})
    write_json(base / "frozen.json", freeze)
    write_json(base / "summary.json", dict(frozen_sha256=sha(base / "frozen.json"),
        figure7_fitted_parameters=0, score=score(base, plan, plan["splits"])))
    print(base / "summary.json", flush=True)


def check(base):
    plan = read_json(base / "plan.json")
    frozen = read_json(base / "frozen.json")
    record = read_json(base / "summary.json")
    if (plan != setup() or record["frozen_sha256"] != sha(base / "frozen.json")
            or frozen["plan"] != sha(base / "plan.json")
            or frozen["curves"] != {name: sha(base / f"{name}.json") for name in plan["splits"]}):
        raise ValueError("A frozen parallel-bed input or curve changed")
    if score(base, plan, plan["splits"]) != record["score"]:
        raise ValueError("Parallel-bed scores changed on replay")
    print("Checked", base, len(frozen["curves"]), "curves", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("base", type=Path, nargs="?")
    action = parser.add_mutually_exclusive_group()
    action.add_argument("--case")
    action.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.case:
        solve(args.base.resolve(), args.case)
    elif args.check:
        check(args.base.resolve())
    else:
        run()
