"""Test batch-qualified Freundlich kinetics at the two unfitted 25 cm flows."""

import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from functools import lru_cache
from pathlib import Path
import sys
from uuid import uuid4

import numpy as np

from examples.adsorption.breakthrough_hypotheses.batch import freundlich, transfer
from examples.adsorption.figure7_retention.run import checked_model, integrate
from examples.adsorption.porosity_matched_column import run as column
from examples.adsorption.transport_accuracy.run import validate_record
from pore2process_demo.bridge import runtime, runtime_versions
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT
from pore2process_demo.provenance import source_manifest

HERE = Path(__file__).resolve().parent
SURFACE_DISPATCH = """SharedContactColumn.surface_value(rhs,beta,qmax,b,law::PublishedFreundlich.Isotherm) =
    PublishedFreundlich.surface_concentration(rhs,beta,law)"""


@lru_cache(maxsize=1)
def load_column_law():
    """Install the alternative equilibrium type in the existing kernel."""
    column.kernel.solver()
    law = freundlich.load()
    runtime().seval(SURFACE_DISPATCH)
    freundlich.verify_law()
    return law.TABLE1


def loading_at_feed(concentration):
    return freundlich.LAW["coefficient_mg_g_per_mg_L_to_exponent"] * (
        1000 * concentration) ** freundlich.LAW["exponent"] / 1000


def validation_basis(inputs):
    """Map the stated ceiling onto the historical Langmuir verification helper."""
    feed = inputs["inlet_concentration_kg_m3"]
    bc = inputs["affinity_m3_kg"] * feed
    loading = loading_at_feed(feed)
    return dict(inputs, capacity_kg_kg=loading * (1 + bc) / bc)


def verify(result, inputs, scenario, condition, spec):
    column.transfer.verify_run(result, validation_basis(inputs), scenario, condition)
    loading = loading_at_feed(inputs["inlet_concentration_kg_m3"])
    if result["loading_range_kg_kg"][1] > loading + 1e-10:
        raise ValueError("Column loading exceeds the published Freundlich loading at feed")
    np.testing.assert_allclose(result["effective_conductance_m3_s"],
        spec["film_matched_clean_bed_conductance_m3_s"], rtol=1e-12)


def plan(batch_base):
    batch = freundlich.check(batch_base)
    if not batch["qualification"]["qualifies"] or batch["figure_parameters_fitted"]:
        raise ValueError("Only batch-qualified equilibria may transfer")
    original, _, _ = transfer.archived_inputs()
    diffusivity = batch["qualification"]["diffusivity_m2_s"]
    cases = {flow: dict(case, parameters=dict(case["parameters"], diffusivity_m2_s=diffusivity),
                        variants={column.PRIMARY: case["variants"][column.PRIMARY]})
             for flow, case in original["cases"].items()}
    return dict(cases=cases, law=freundlich.LAW, batch_run=batch_base.relative_to(ROOT).as_posix(),
                batch_sha256=sha(batch_base / "summary.json"), archived_column_sha256=sha(transfer.REFERENCE))


def level_inputs(p, spec, level):
    shells, step, early = level
    return dict(p, film_m_s=spec["film_m_s"], radial_cells=shells, dt_s=step, startup_dt_s=early)


def solve(base, flow):
    law = load_column_law()
    condition = read_json(base / "plan.json")["cases"][flow]
    spec = condition["variants"][column.PRIMARY]
    p, scenario = condition["parameters"], spec["access"]
    probe = ROOT / spec["probe_directory"]
    previous = None
    for index, level in enumerate(column.LEVELS):
        inputs = level_inputs(p, spec, level)
        result = column.kernel.solve(condition["flow_mL_min"], probe,
                                     dict(inputs, linear_partition_m3_kg=law), scenario)
        verify(result, inputs, scenario, condition, spec)
        filename = f"{flow}_{index}.json"
        write_json(base / filename, dict(inputs=inputs, scenario=scenario, law=freundlich.LAW, result=result))
        print(flow, index, "wall seconds", round(result["wall_s"], 1), flush=True)
        if previous is not None:
            check = column.full_refinement(previous, result, p["inlet_concentration_kg_m3"])
            if check["passed"] and check["startup"]["passed"]:
                write_json(base / f"{flow}_levels.json", dict(files=[f"{flow}_{i}.json" for i in range(index+1)],
                           selected=filename, refinement=check))
                return
        previous = result
    raise ValueError("Freundlich transfer did not meet the numerical refinement screens")


def predict(base, batch_base):
    setup = plan(batch_base)
    write_json(base / "plan.json", setup)
    (base / "solver.jl").write_text(column.kernel.source())
    (base / "state_equilibrium.jl").write_text(freundlich.JULIA.read_text() + "\n" + SURFACE_DISPATCH + "\n")
    print("Frozen Freundlich column inputs", base, flush=True)

    def process(flow):
        command = [sys.executable, "-m", __spec__.name, str(base), "--flow", flow]
        timing = column.execute(command, base / f"{flow}.log", timeout=3600)
        print(flow, "seconds", round(timing["wall_s"], 1), flush=True)

    with ThreadPoolExecutor(max_workers=2) as pool:
        list(pool.map(process, column.FLOWS))
    write_json(base / "frozen.json", dict(
        sha256={p.name: sha(p) for p in base.iterdir() if p.is_file()}))


def curves(base):
    sealed = read_json(base / "frozen.json")
    for filename, digest in sealed["sha256"].items():
        if sha(base / filename) != digest:
            raise ValueError(f"Freundlich column artifact changed: {filename}")
    setup = read_json(base / "plan.json")
    if (setup != plan(ROOT / setup["batch_run"]) or
            (base / "solver.jl").read_text() != column.kernel.source() or
            (base / "state_equilibrium.jl").read_text() !=
            freundlich.JULIA.read_text() + "\n" + SURFACE_DISPATCH + "\n"):
        raise ValueError("Freundlich transfer setup differs on replay")
    cases = {}
    for flow, condition in setup["cases"].items():
        levels = read_json(base / f"{flow}_levels.json")
        spec = condition["variants"][column.PRIMARY]
        p = condition["parameters"]
        if (not 2 <= len(levels["files"]) <= len(column.LEVELS) or
                levels["files"][-1] != levels["selected"]):
            raise ValueError("The Freundlich column resolution check is incomplete")
        previous = None
        for filename, level in zip(levels["files"], column.LEVELS):
            current = read_json(base / filename)
            inputs = level_inputs(p, spec, level)
            if (current["inputs"] != inputs or current["scenario"] != spec["access"] or
                    current["law"] != freundlich.LAW):
                raise ValueError("A frozen Freundlich column input changed")
            result = current["result"]
            verify(result, inputs, spec["access"], condition, spec)
            check = column.full_refinement(previous, result, p["inlet_concentration_kg_m3"]) if previous else None
            previous = result
        if check != levels["refinement"] or not (check["passed"] and check["startup"]["passed"]):
            raise ValueError("The Freundlich column resolution check changed")
        cases[flow] = (current, check)
    return cases


def compare(base):
    predictions = curves(base)
    data = column.source.markers()
    keys = {"f4_q233_l25": "q233", "f4_q625_l25": "q625", "f5_q233_l25": "q233"}
    params = {flow: saved["inputs"] for flow, (saved, _) in predictions.items()}
    rates = {flow: p["flow_m3_s"] * p["inlet_concentration_kg_m3"] for flow, p in params.items()}
    common_dose = min(data[key]["time"][-1] * rates[flow] / params[flow]["dry_mass_kg"]
                      for key, flow in keys.items() if key.startswith("f4"))
    comparisons = {}
    for key, flow in keys.items():
        p, rate, markers = params[flow], rates[flow], data[key]
        mass = p["dry_mass_kg"]
        saved, check = predictions[flow]
        common_time = min(float(common_dose * mass / rate), markers["time"][-1])
        observed = column.source.observation(markers, rate, mass, common_time)
        basis = dict(p, feed_rate_kg_s=rate)
        metrics = column.curve_metrics(saved, basis, markers, observed, common_time)
        t, c = np.r_[0., markers["time"]], np.r_[0., markers["fraction"]]
        mid, stop = observed["t50_s"], observed["end_s"]
        observed["late_retained_kg_kg"] = (integrate(t, c, rate, stop) - integrate(t, c, rate, mid)) / mass
        model = checked_model(saved, basis)
        retained = np.interp([mid, stop], model["time_s"], model["retained_kg"])
        comparisons[key] = dict(observed=observed, prediction=dict(**metrics,
            late_retained_kg_kg=float(retained[1]-retained[0])/mass, radial_cells=p["radial_cells"],
            time_step_s=p["dt_s"], diffusivity_m2_s=p["diffusivity_m2_s"], numerical=check))
    return comparisons


def run(batch_base):
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ") + "-" + uuid4().hex[:8]
    base = freundlich.langmuir.OUTPUT / "freundlich_transfer" / stamp
    base.mkdir(parents=True)
    predict(base, batch_base)
    comparisons = compare(base)
    files = [Path(__file__), Path(freundlich.__file__), freundlich.JULIA, HERE / "freundlich_protocol.md",
             column.source.EXTRACTION, column.source.PIXELS, transfer.REFERENCE,
             Path(column.kernel.__file__), column.kernel.access.ORIGINAL, Path(transfer.__file__),
             ROOT / "julia/adsorption/FilmParticle.jl", ROOT / "julia/adsorption/SharedContactColumn.jl"]
    write_json(base / "summary.json", dict(schema="pore2process.breakthrough-freundlich-transfer.v1",
        frozen_sha256=sha(base / "frozen.json"), figure_parameters_fitted=0, comparisons=comparisons,
        source_sha256=source_manifest(files), versions=runtime_versions(),
        artifacts_sha256={p.name: sha(p) for p in base.iterdir() if p.is_file()}))
    print(base / "summary.json", flush=True)


def check(base):
    saved = validate_record(base / "summary.json")
    if saved["frozen_sha256"] != sha(base / "frozen.json") or saved["comparisons"] != compare(base):
        raise ValueError("Freundlich predictions or scores changed on replay")
    print(dict(sources=len(saved["source_sha256"]), artifacts=len(saved["artifacts_sha256"])))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("base", type=Path, help="Frozen Freundlich batch run, or column transfer to check")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--flow", choices=column.FLOWS)
    mode.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.flow:
        solve(args.base.resolve(), args.flow)
    elif args.check:
        check(args.base.resolve())
    else:
        run(args.base.resolve())
