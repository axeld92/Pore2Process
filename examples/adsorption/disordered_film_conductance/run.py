"""Estimate clean liquid-film conductance from the saved N160 TRT field."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from uuid import uuid4

import numpy as np

from examples.adsorption.disordered_early_bypass import run as bypass
from examples.adsorption.disordered_early_bypass.model import plateau_bound
from examples.adsorption.disordered_exchange.run import FOLDER, PARENT, PARENT_SHA, parent_inputs
from examples.adsorption.disordered_section_exchange.run import TRACKER, verify_saved_field
from examples.adsorption.dispersion_holdout.run import execute
from examples.adsorption.transport_accuracy.run import validate_record
from pore2process_demo.bridge import runtime_versions
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT
from pore2process_demo.provenance import source_manifest

from .model import effective_coefficient, required_coefficient, spherical_area_by_half, stats

HERE = Path(__file__).resolve().parent
OUTPUT = ROOT / "results/adsorption/disordered_film_conductance"
BYPASS = bypass.OUTPUT / "runs/20261117T103314Z-b73ee872/summary.json"
FLOWS = bypass.FLOWS
STEP_CHANGE_LIMIT = .03
RESIDENCE_LIMIT = .02


def inputs():
    reference = read_json(bypass.HERE / "reference_summary.json")
    if sha(BYPASS) != reference["run"]["sha256"]:
        raise ValueError("Previous clean-film result changed")
    plan, _, probes = bypass.inputs()
    prior = validate_record(BYPASS)
    n160, tracer = parent_inputs()
    verify_saved_field(n160, tracer)
    native_plan = read_json(PARENT / "plan.json")
    centres = PARENT / "disordered_centres_m.txt"
    if sha(centres) != n160["artifacts_sha256"][centres.name]:
        raise ValueError("Archived solid geometry changed")
    probe_plan = read_json(probes / "plan.json")
    return plan, prior, native_plan, probes, probe_plan


def probe_specs(base, probe_plan, step):
    specs = {}
    for flow in FLOWS:
        source = probe_plan["cases"][flow]
        command = source["command"].copy()
        command[6] = str(base / f"{flow}.f64")
        command[9] = str(step)
        command[11] = str(source["training_count"])
        specs[flow] = dict(command=command, count=int(command[11]), dt_s=step)
    return specs


def probe_level(base, probe_plan, step):
    folder = base / f"{round(step * 1e6)}us"
    folder.mkdir()
    specs = probe_specs(folder, probe_plan, step)
    write_json(folder / "plan.json", specs)

    def trace(item):
        flow, spec = item
        timing = execute(spec["command"], folder / f"{flow}.log", timeout=2400)
        print(flow, spec["dt_s"], timing, flush=True)
        return flow, timing

    with ThreadPoolExecutor(max_workers=2) as pool:
        timing = dict(pool.map(trace, specs.items()))
    write_json(folder / "execution.json", timing)


def load(base, flow, count):
    data = np.fromfile(base / f"{flow}.f64", dtype="<f8").reshape(-1, 16)
    if len(data) != count:
        raise ValueError("Unexpected marker count")
    return data


def geometry(native_plan, probe_plan):
    geometry = native_plan["geometry"]["disordered"]
    compact = read_json(PARENT / FOLDER / "compact/summary.json")
    count = np.array(compact["period_end_counts"])
    length = geometry["extent_m"][0]
    dx = float(probe_plan["cases"][FLOWS[0]]["command"][7])
    diameter = geometry["diameter_m"]
    centres = np.loadtxt(PARENT / "disordered_centres_m.txt")
    if compact["cut_index"] != 0 or centres.shape != (geometry["particle_count"], 3):
        raise ValueError("Unexpected centre coordinates or physical inlet cut")
    areas = spherical_area_by_half(centres[:, 0], diameter, length)
    np.testing.assert_allclose(areas.sum(), len(centres) * np.pi * diameter**2, rtol=1e-13)
    volumes = dx**3 * np.array([count[3], count[7] - count[3], count[7]])
    return dict(dx_m=dx, length_m=length, diameter_m=diameter, sphere_count=len(centres),
                analytical_porosity=geometry["analytical_porosity"], voxel_void_m3=volumes.tolist(),
                external_area_m2=[*areas.tolist(), float(areas.sum())],
                native_throughput_m3_s=compact["throughput_m3_s"])


def column_setup(plan, prior, spec, geo, flow):
    condition = plan["cases"][flow + "_l25"]
    p = condition["parameters"]
    scenario = condition["scenarios"]["f50_slow01"]
    baseline = prior["clean_film"][flow]
    Q, V, N = p["flow_m3_s"], p["liquid_volume_m3"], condition["sections"]
    area = 6 * p["dry_mass_kg"] / (p["particle_density_kg_m3"] * p["particle_diameter_m"])
    geometric_scale = condition["bed_cm"] / 100 / N / (geo["length_m"] / 8)
    # One section in the two descriptions has the same advective Péclet
    # number. Map archive velocities to the column's interstitial velocity.
    u_archive = geo["native_throughput_m3_s"] * geo["length_m"] / geo["voxel_void_m3"][-1]
    u_column = Q * condition["bed_cm"] / 100 / V
    velocity_scale = u_column / u_archive
    D_archive = spec["archive_coordinate_diffusivity_m2_s"]
    D_column = u_column * condition["bed_cm"] / 100 / N / condition["Pe"]
    native_pe = u_archive * (geo["length_m"] / 8) / D_archive
    np.testing.assert_allclose(native_pe, condition["Pe"], rtol=1e-5)
    np.testing.assert_allclose(D_column, D_archive * velocity_scale * geometric_scale, rtol=1e-5)
    f = scenario["inputs"]["direct_fraction"]
    R = scenario["added_resistance_s_m"]
    np.testing.assert_allclose(area * effective_coefficient(p["film_m_s"], f, R),
                               baseline["total_film_conductance_m3_s"], rtol=1e-13)
    required = Q * baseline["minimum_ntu_for_target"]
    porosity = V / (V + p["dry_mass_kg"] / p["particle_density_kg_m3"])
    return dict(flow_mL_min=condition["flow_mL_min"], flow_m3_s=Q, sections=N,
                blocks=p["blocks"], terminal_sections=p["terminal_sections"],
                dry_carbon_mass_kg=p["dry_mass_kg"], external_area_m2=area, liquid_volume_m3=V, porosity=porosity,
                length_scale=geometric_scale, velocity_scale=velocity_scale,
                physical_diffusivity_m2_s=D_column, archive_diffusivity_m2_s=D_archive, Pe=native_pe,
                direct_fraction=f, added_resistance_s_m=R, baseline_film_coefficient_m_s=p["film_m_s"],
                baseline_effective_conductance_m3_s=baseline["total_film_conductance_m3_s"],
                required_effective_conductance_m3_s=required,
                required_unrestricted_coefficient_m_s=required_coefficient(required, area, f, R))


def transfer(k, geo, column, area_based):
    """Map an archive coefficient to the fixed bed and retain access resistance."""
    area = column["external_area_m2"]
    V = column["liquid_volume_m3"]
    k = np.asarray(k) * column["velocity_scale"]
    area_per_void = geo["external_area_m2"][-1] / geo["voxel_void_m3"][-1]
    if not area_based:
        k *= area_per_void * V / (column["length_scale"] * area)
    G = area * effective_coefficient(k, column["direct_fraction"], column["added_resistance_s_m"])
    return G.tolist()


def summarize(raw, geo, column, seed):
    result = stats(raw, np.array(geo["voxel_void_m3"]), np.array(geo["external_area_m2"]), seed)
    inert_expected = geo["voxel_void_m3"][-1] / geo["native_throughput_m3_s"]
    error = abs(result["inert_transit_s"] / inert_expected - 1)
    if error > RESIDENCE_LIMIT:
        raise ValueError("Unreacted residence time differs from pore volume / flow")
    raw_k = result["archive_film_coefficient_m_s"][-1]
    interval = np.array(result["sampling_95pct_m_s"])[:, -1]
    result.update(inert_time_relative_error=error,
        film_coefficient_m_s=(np.array(result["archive_film_coefficient_m_s"]) *
                              column["velocity_scale"]).tolist(),
        direct_film_sampling_95pct_m_s=(interval * column["velocity_scale"]).tolist(),
        area_transfer_m3_s=transfer(raw_k, geo, column, True),
        area_transfer_sampling_95pct_m3_s=transfer(interval, geo, column, True),
        volume_transfer_m3_s=transfer(raw_k, geo, column, False),
        volume_transfer_sampling_95pct_m3_s=transfer(interval, geo, column, False))
    return result


def refinement(old, paired, fine, column):
    target = column["required_effective_conductance_m3_s"]
    change = abs(fine["area_transfer_m3_s"] / paired["area_transfer_m3_s"] - 1)
    spans = {}
    for key in ("area_transfer", "volume_transfer"):
        intervals = [row[key + "_sampling_95pct_m3_s"] for row in (old, fine)]
        lo = min(a for a, _ in intervals)
        hi = max(b for _, b in intervals)
        spans[key] = dict(range_m3_s=[lo, hi], crosses_target=lo <= target <= hi)
    return dict(paired_relative_change=change, spans=spans,
                needs_finer_step=change > STEP_CHANGE_LIMIT or
                any(item["crosses_target"] for item in spans.values()))


def analyse(base, inputs):
    plan, prior, native, probes, probe_plan = inputs
    geo = geometry(native, probe_plan)
    cases = {}
    for index, flow in enumerate(FLOWS):
        spec = probe_plan["cases"][flow]
        column = column_setup(plan, prior, spec, geo, flow)
        data = np.fromfile(probes / flow / "markers.f64", dtype="<f8").reshape(-1, 16)
        count = spec["training_count"]
        if len(data) != count + spec["operator_probe_count"]:
            raise ValueError("Unexpected archived ensemble size")
        seed = 1408 + index
        coarse = summarize(data, geo, column, seed)
        matched = summarize(data[:count], geo, column, seed)
        half = summarize(load(base / "40us", flow, count), geo, column, seed)
        refinements = dict(half=refinement(coarse, matched, half, column))
        levels = dict(archive_80us=coarse, matched_80us=matched, matched_40us=half)
        if (base / "20us/execution.json").exists():
            fine = summarize(load(base / "20us", flow, count), geo, column, seed)
            levels["matched_20us"] = fine
            refinements["quarter"] = refinement(half, half, fine, column)
        cases[flow] = dict(column=column, estimates=levels, refinement=refinements)
    return dict(geometry=geo, cases=cases)


def wall_survival(stages, blocks, tail):
    count = max(int(max(i.max(), j.max())) for i, j, _, _ in stages) + 1
    state = np.eye(1, count).ravel()
    remaining = [1.]
    for i, j, _, weights in stages:
        state = np.bincount(j, weights=state[i] * weights, minlength=count)
        remaining.append(float(state.sum()))
    return remaining[-1]**blocks * remaining[tail]


def clean_film(analysis, probes):
    for flow, case in analysis["cases"].items():
        column = case["column"]
        model = bypass.kernel.transport(column["flow_mL_min"], probes)
        N, Q = column["sections"], column["flow_m3_s"]
        blocks, tail = column["blocks"], column["terminal_sections"]
        finest = next(reversed(case["estimates"].values()))

        def row(G):
            ntu = G / Q
            bound = plateau_bound(ntu, N)
            wall = wall_survival(bypass.film_stages(model, ntu, N, "wall"), blocks, tail)
            uniform = wall_survival(bypass.film_stages(model, ntu, N, "uniform_per_flow"), blocks, tail)
            np.testing.assert_allclose(uniform, bound, atol=1e-11, rtol=0.)
            if wall < bound - 1e-11:
                raise ValueError("Wall-weighted survival violates fixed-conductance bound")
            return dict(effective_conductance_m3_s=G, best_escape=bound, wall_escape=wall,
                        best_below_cutoff=bound <= bypass.TARGET, wall_below_cutoff=wall <= bypass.TARGET)

        case["clean_film"] = dict(baseline=row(column["baseline_effective_conductance_m3_s"]),
            area_transfer=row(finest["area_transfer_m3_s"]),
            area_sampling_95pct_wall_escape=sorted(row(G)["wall_escape"] for G in
                                                    finest["area_transfer_sampling_95pct_m3_s"]),
            volume_transfer=row(finest["volume_transfer_m3_s"]))
    return analysis


def verify_probe_plans(base, probe_plan):
    if read_json(base / "archive.json") != archive_snapshot():
        raise ValueError("Tracer probe archive does not match its parent snapshots")
    for folder in ("40us", "20us"):
        if folder == "20us" and not (base / folder).exists():
            continue
        path = base / folder
        step = int(folder[:-2]) / 1e6
        if read_json(path / "plan.json") != probe_specs(path, probe_plan, step):
            raise ValueError("Fresh tracer settings differ from preregistered settings")
        if any(item["returncode"] != 0 for item in read_json(path / "execution.json").values()):
            raise ValueError("A fresh tracer solve did not complete")


def finish(base):
    if (base / "summary.json").exists():
        raise ValueError("Use --check to verify an already sealed run")
    setup = inputs()
    verify_probe_plans(base, setup[-1])
    analysis = analyse(base, setup)
    needs_finer = any(item["refinement"]["half"]["needs_finer_step"]
                      for item in analysis["cases"].values())
    if needs_finer and not (base / "20us").exists():
        probe_level(base, setup[-1], 2e-5)
        analysis = analyse(base, setup)
    verify_probe_plans(base, setup[-1])
    clean_film(analysis, setup[-2])
    sources = source_manifest([HERE / name for name in ("run.py", "model.py", "protocol.md")] +
        [HERE.parent / name for name in TRACKER] +
        [bypass.HERE / name for name in ("run.py", "model.py", "reference_summary.json")])
    path = base / "summary.json"
    write_json(path, dict(schema="pore2process.disordered-film-conductance.v1", run_id=base.name,
        archive=read_json(base / "archive.json"), source_sha256=sources,
        artifacts_sha256={file.relative_to(base).as_posix(): sha(file)
                          for file in sorted(base.rglob("*")) if file.is_file()},
        analysis=analysis, versions=runtime_versions(),
        new_native_solves=0, new_adsorption_runs=0, figure7_fits=0))
    print(path, flush=True)


def check(path):
    record = validate_record(path)
    setup = inputs()
    verify_probe_plans(path.parent, setup[-1])
    if record["analysis"] != clean_film(analyse(path.parent, setup), setup[-2]):
        raise ValueError("Film conductance or clean-film transfer differs on replay")
    print(dict(sources=len(record["source_sha256"]), artifacts=len(record["artifacts_sha256"])))


def archive_snapshot():
    return dict(n160_sha256=PARENT_SHA, prior_bypass_sha256=sha(BYPASS),
                probe_sha256=bypass.PROBES_SHA, column_sha256=bypass.COLUMN_SHA)


def start():
    _, _, _, _, probe_plan = inputs()
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ") + "-" + uuid4().hex[:8]
    base = OUTPUT / "runs" / stamp
    base.mkdir(parents=True)
    write_json(base / "archive.json", archive_snapshot())
    print(base, flush=True)
    probe_level(base, probe_plan, 4e-5)
    return base


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--probe", action="store_true")
    mode.add_argument("--finish", type=Path)
    mode.add_argument("--check", type=Path)
    args = parser.parse_args()
    if args.probe:
        start()
    elif args.finish:
        finish(args.finish.resolve())
    else:
        check(args.check.resolve())
