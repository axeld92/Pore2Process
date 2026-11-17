"""Compare clean-wall uptake and fixed-mass column escape at matched porosity."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from math import isclose
from pathlib import Path
from time import sleep

import numpy as np

from examples.adsorption.carbon_access_transfer import kernel
from examples.adsorption.disordered_early_bypass.model import plateau_bound, required_ntu
from examples.adsorption.disordered_film_conductance.model import (
    effective_coefficient, spherical_area_by_half, stats,
)
from examples.adsorption.disordered_film_conductance.run import transfer, wall_survival
from examples.adsorption.disordered_section_exchange.run import TRACKER
from examples.adsorption.disordered_travel_memory.run import HERE as MEMORY
from examples.adsorption.dispersion_holdout.run import execute
from examples.adsorption.transport_accuracy.run import validate_record
from pore2process_demo.bridge import runtime_versions
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT
from pore2process_demo.provenance import source_manifest
from pore2process_demo.runtime.openlb import verify_vendor

from . import native

HERE = Path(__file__).resolve().parent
TARGET = .05
FILM_LIMIT = .03
ESCAPE_LIMIT = .02


def inputs(base, grid, flow):
    config = read_json(base / "plan.json")
    if sha(base / "centres_m.txt") != config["centres_sha256"]:
        raise ValueError("Packing centers changed")
    c = config["cases"][flow]
    folder = base / f"{flow}_N{grid}"
    face = read_json(folder / "faces/summary.json")
    compact = read_json(folder / "compact/summary.json")
    if compact["face_summary_sha256"] != sha(folder / "faces/summary.json"):
        raise ValueError("Compact grid no longer matches its face record")
    dx = face["info"]["raw"]["dx_m"]
    native_flow = compact["throughput_m3_s"]
    N = c["sections"]
    length = config["extent_m"][0]
    counts = compact["period_end_counts"]
    void = dx**3 * np.array([counts[3], counts[7] - counts[3], counts[7]])
    local_velocity = native_flow * length / void[-1]
    local_diffusivity = local_velocity * (length / native.PERIODS) / c["Pe"]
    column_velocity = c["interstitial_m_s"]
    column_section = c["column_length_m"] / N
    discrepancy = abs(column_section * column_velocity / c["Pe"] / 1e-9 - 1.)
    if discrepancy > 1e-12:
        raise ValueError("Péclet reference no longer matches the frozen column diffusivity")
    centers = np.loadtxt(base / "centres_m.txt")
    shift = compact["cut_offset_m"]
    areas = spherical_area_by_half(centers[:, 0] - shift, native.DIAMETER_M, length)
    np.testing.assert_allclose(areas.sum(), len(centers) * np.pi * native.DIAMETER_M**2, rtol=1e-13)
    geo = dict(dx_m=dx, length_m=length, porosity=face["info"]["porosity"], cut_offset_m=shift,
        native_throughput_m3_s=native_flow, voxel_void_m3=void.tolist(), end_counts=counts,
        external_area_m2=[*areas.tolist(), float(areas.sum())], archive_diffusivity_m2_s=local_diffusivity)
    area = 6 * c["dry_carbon_mass_kg"] / (c["particle_density_kg_m3"] * native.DIAMETER_M)
    f = c["access"]["inputs"]["direct_fraction"]
    R = c["access"]["added_resistance_s_m"]
    column = dict(flow_mL_min=c["flow_mL_min"], flow_m3_s=c["flow_m3_s"], sections=N, Pe=c["Pe"],
        liquid_volume_m3=c["column_liquid_m3"], dry_carbon_mass_kg=c["dry_carbon_mass_kg"],
        external_area_m2=area, direct_fraction=f, added_resistance_s_m=R,
        physical_diffusivity_m2_s=1e-9, length_scale=column_section / (length / native.PERIODS),
        velocity_scale=column_velocity / local_velocity,
        baseline_effective_conductance_m3_s=area * effective_coefficient(c["baseline_film_m_s"], f, R),
        required_effective_conductance_m3_s=c["flow_m3_s"] * required_ntu(TARGET, N))
    return folder, geo, column, face


def probe(base, grid, flow, step, count, threads=None):
    if step <= 0 or count <= 0 or count % 3:
        raise ValueError("Use a positive time step and marker count divisible by three")
    folder, geo, column, _ = inputs(base, grid, flow)
    compact = folder / "compact"
    if any(sha(compact / key) != value for key, value in
           read_json(compact / "summary.json")["files"].items()):
        raise ValueError("Compact tracer payload changed")
    _, _, prior = native.bypass.inputs()
    frozen = read_json(prior / "plan.json")["cases"][flow]
    out = base / f"tracer_N{grid}_{step}us" / flow
    dt = step / 1e6
    stops = 15.
    runtime = frozen["command"][:4]
    if threads is not None:
        runtime[1] = f"--threads={threads}"
    command = [*runtime, str(native.PACKING / "track.jl"), str(folder / "compact"), str(out / "markers.f64"),
        *map(str, (geo["dx_m"], geo["archive_diffusivity_m2_s"], dt, stops, count, 4928073,
            native.grid_shape(grid)[0] // native.PERIODS, ",".join(map(str, geo["end_counts"]))))]
    request = dict(command=command, marker_count=count, grid=grid, time_step_s=dt,
                   units="archive metres and seconds", Pe=column["Pe"],
                   compact_sha256=sha(folder / "compact/summary.json"))
    try:
        out.mkdir(parents=True)
    except FileExistsError:
        for _ in range(360):
            if (out / "execution.json").is_file():
                if read_json(out / "request.json") != request:
                    raise ValueError("Existing tracer solve uses a different request")
                markers(out)
                return
            sleep(2)
        raise ValueError(f"Existing tracer solve has no execution record in {out}")
    write_json(out / "request.json", request)
    timing = execute(command, out / "tracker.log", timeout=3600)
    raw = markers(out)
    train = count // 3
    raw[:train, :2].copy(order="C").tofile(out / "training.f64")
    raw[train:].copy(order="C").tofile(out / "operator.f64")
    raw[train:, :8].copy(order="C").tofile(out / "passages.f64")
    red = [*runtime, str(MEMORY / "solve.jl"), str(out / "training.f64"),
           str(out / "passages.f64"), str(out)]
    elapsed = execute(red, out / "reduction.log", timeout=180)
    write_json(out / "execution.json", dict(tracker=timing, reduction=elapsed, reduction_command=red))
    print(f"{flow} N{grid}, {step} us tracer complete", timing, flush=True)


def markers(out):
    size = read_json(out / "request.json")["marker_count"]
    raw = np.fromfile(out / "markers.f64", dtype="<f8").reshape(size, 16)
    if not np.isfinite(raw[:, :8]).all():
        raise ValueError(f"Incomplete inert tracer recovery in {out}")
    return raw


def score_case(base, grid, flow, step):
    folder, geo, column, face = inputs(base, grid, flow)
    root = base / f"tracer_N{grid}_{step}us"
    traces = root / flow
    request = read_json(traces / "request.json")
    if (request["grid"], request["time_step_s"], request["Pe"], request["compact_sha256"]) != (
            grid, step / 1e6, column["Pe"], sha(folder / "compact/summary.json")):
        raise ValueError("Tracer request does not match the fixed carrier")
    raw = markers(traces)
    result = stats(raw, np.array(geo["voxel_void_m3"]), np.array(geo["external_area_m2"]), 427)
    expected = geo["voxel_void_m3"][-1] / geo["native_throughput_m3_s"]
    residence_error = abs(result["inert_transit_s"] / expected - 1.)
    if residence_error > .02:
        raise ValueError("Tracer residence does not match pore volume / throughput")
    coefficient = result["archive_film_coefficient_m_s"][-1]
    bootstrap = np.array(result["sampling_95pct_m_s"])[:, -1]
    model = kernel.transport(column["flow_mL_min"], root)
    N, Q = column["sections"], column["flow_m3_s"]
    blocks, tail = divmod(N, native.PERIODS)

    def wall(G, allocate="wall"):
        stages = native.bypass.film_stages(model, G / Q, N, allocate)
        return wall_survival(stages, blocks, tail)

    full = {}
    for name, by_area in (("carbon_area", True), ("liquid_volume", False)):
        G = transfer(coefficient, geo, column, by_area)
        interval = transfer(bootstrap, geo, column, by_area)
        best = plateau_bound(G / Q, N)
        np.testing.assert_allclose(wall(G, "uniform_per_flow"), best, rtol=0., atol=1e-11)
        escape = wall(G)
        if escape < best - 1e-11:
            raise ValueError("Projected escape violates the optimal-allocation bound")
        full[name] = dict(effective_conductance_m3_s=G, sampling_95pct_conductance_m3_s=interval,
            optimal_escape=best, wall_weighted_escape=escape,
            sampling_95pct_wall_weighted_escape=sorted(wall(g) for g in interval))
    return dict(column=column, geometry=geo, native=face["info"], balanced_faces=face["checks"],
        tracer_step_s=step / 1e6, markers=len(raw), captured=result["captured"][-1],
        inert_mean_time_s=result["inert_transit_s"], inert_residence_relative_error=residence_error,
        film_coefficient_m_s=[float(x * column["velocity_scale"]) for x in result["archive_film_coefficient_m_s"]],
        film_coefficient_sampling_95pct_m_s=(bootstrap * column["velocity_scale"]).tolist(),
        capture_fraction=result["captured"][-1] / len(raw),
        full_carbon_perfect_sink_escape=(1 - result["captured"][-1] / len(raw))**blocks
            * float(np.mean(~np.isfinite(raw[:, 8:8 + tail]).any(axis=1))),
        clean_film=full)


def compare(coarse, fine):
    first = coarse["clean_film"]["carbon_area"]
    second = fine["clean_film"]["carbon_area"]
    G_change = abs(second["effective_conductance_m3_s"] / first["effective_conductance_m3_s"] - 1.)
    escape_change = max(abs(fine["clean_film"][key]["wall_weighted_escape"]
                            - coarse["clean_film"][key]["wall_weighted_escape"])
                        for key in fine["clean_film"])
    return dict(relative_conductance_change=G_change, absolute_escape_change=escape_change,
                passed=G_change <= FILM_LIMIT and escape_change <= ESCAPE_LIMIT)


def analyse(base, grids, step_pairs):
    levels = {f"N{n}_{dt}us": {f: score_case(base, n, f, dt) for f in native.FLOWS}
              for n, dt in step_pairs}
    coarse, fine = grids[-2:]
    time_pair = levels[f"N{fine}_40us"], levels[f"N{fine}_20us"]
    grid_pair = levels[f"N{coarse}_20us"], levels[f"N{fine}_20us"]
    flow = native.FLOWS[-1]
    checks = {f: dict(time=compare(*(item[f] for item in time_pair)),
                      grid=compare(*(item[f] for item in grid_pair))) for f in native.FLOWS}
    high = grid_pair[-1][flow]["clean_film"]
    band = {}
    envelope = {}
    for method, row in high.items():
        lo, hi = row["sampling_95pct_wall_weighted_escape"]
        band[method] = "below" if hi < TARGET else "above" if lo > TARGET else "straddles"
        intervals = [case[flow]["clean_film"][method]["sampling_95pct_wall_weighted_escape"]
                     for case in (grid_pair[0], time_pair[0], grid_pair[1])]
        envelope[method] = [min(a for a, _ in intervals), max(b for _, b in intervals)]
    numerically_qualified = all(item["passed"] for item in checks[flow].values())
    consistent_below = all(hi < TARGET for _, hi in envelope.values())
    consistent_above = all(lo > TARGET for lo, _ in envelope.values())
    screening = dict(target=TARGET, high_flow_point_below={k: r["wall_weighted_escape"] < TARGET
                     for k, r in high.items()}, high_flow_intervals=band,
        high_flow_numerically_qualified=numerically_qualified, high_flow_comparison_envelope=envelope,
        high_flow_below_across_checked_runs=consistent_below, high_flow_above_across_checked_runs=consistent_above,
        high_flow_qualified_below=numerically_qualified and consistent_below,
        high_flow_qualified_above=numerically_qualified and consistent_above)
    return dict(levels=levels, refinement=checks, screening=screening)


def seal(base, grids):
    destination = base / "summary.json"
    if destination.exists():
        raise ValueError("Use --check to verify a sealed cohort")
    pairs = [(n, 20) for n in grids] + [(grids[-1], 40)]
    result = analyse(base, grids, pairs)
    files = [*HERE.glob("*.py"), HERE / "protocol.md", *list((native.PACKING / "openlb").glob("*")),
             *(HERE.parent / p for p in TRACKER), MEMORY / "solve.jl"]
    sources = source_manifest([p for p in files if p.is_file()])
    archived = {str(p.relative_to(base)): sha(p) for p in sorted(base.rglob("*")) if p.is_file()}
    write_json(destination, dict(schema="pore2process.porosity-matched-film.v1", run_id=base.name,
        grids=grids, pairs=pairs, analysis=result, source_sha256=sources, artifacts_sha256=archived,
        vendor=verify_vendor(), versions=runtime_versions(), figure7_fits=0))
    print(destination, flush=True)


def same_analysis(saved, current):
    """Replayed reductions may differ by floating-point roundoff."""
    if isinstance(saved, dict):
        return saved.keys() == current.keys() and all(same_analysis(v, current[k]) for k, v in saved.items())
    if isinstance(saved, list):
        return len(saved) == len(current) and all(same_analysis(a, b) for a, b in zip(saved, current))
    if isinstance(saved, float):
        return isclose(saved, current, rel_tol=1e-12)
    return saved == current


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("base", type=Path)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--trace", type=int, nargs=3, metavar=("N", "STEP_US", "MARKERS"))
    mode.add_argument("--assess", type=int, nargs=2, metavar=("N", "STEP_US"))
    mode.add_argument("--seal", type=int, nargs="+")
    mode.add_argument("--check", action="store_true")
    parser.add_argument("--flow", choices=native.FLOWS)
    parser.add_argument("--threads", type=int, help="Julia threads per tracer flow")
    args = parser.parse_args()
    if args.threads is not None and (args.threads < 1 or not args.trace):
        parser.error("--threads is a positive count for --trace")
    if args.flow and not args.trace:
        parser.error("--flow applies only to --trace")
    if args.seal and (len(args.seal) < 2 or args.seal != sorted(set(args.seal))):
        parser.error("--seal needs at least two increasing grid sizes")
    base = args.base.resolve()
    if args.trace:
        with ThreadPoolExecutor(max_workers=2) as pool:
            flows = [args.flow] if args.flow else native.FLOWS
            jobs = [pool.submit(probe, base, args.trace[0], f, *args.trace[1:], args.threads)
                    for f in flows]
            for job in jobs:
                job.result()
    elif args.seal:
        seal(base, args.seal)
    elif args.assess:
        grid, step = args.assess
        for f in native.FLOWS:
            case = score_case(base, grid, f, step)
            print(f, {key: case[key] for key in ("film_coefficient_m_s", "capture_fraction", "clean_film")})
    else:
        record = validate_record(base / "summary.json")
        if not same_analysis(record["analysis"], analyse(base, record["grids"], record["pairs"])):
            raise ValueError("Reported result differs on scientific replay")
        print(dict(sources=len(record["source_sha256"]), artifacts=len(record["artifacts_sha256"])))


if __name__ == "__main__":
    main()
