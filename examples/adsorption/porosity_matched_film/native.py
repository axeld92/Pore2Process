"""Prepare and run the separate porosity-matched stock OpenLB TRT case."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import os
from pathlib import Path
import shutil
import subprocess
from time import perf_counter
from uuid import uuid4

import numpy as np

from examples.adsorption.disordered_early_bypass import run as bypass
from examples.adsorption.disordered_packing.mesh import fluid_mask
from examples.adsorption.disordered_packing.native import HERE as PACKING
from examples.adsorption.disordered_packing.packing import DIAMETER_M, PERIODS, minimum_distance, offsets, ordered
from examples.adsorption.disordered_refinement.native import interpolate
from examples.adsorption.transport_root_cause.hydraulics.stock_omp import build_case_omp
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT

HERE = Path(__file__).resolve().parent
OUTPUT = ROOT / "results/adsorption/porosity_matched_film"
BUILD = ROOT / ".local/adsorption/porosity_matched_film"
POROSITY = .427
ASPECT = 6 / 5
SOLID_VOLUME = 4 * PERIODS * np.pi * DIAMETER_M**3 / 6
SIDE = np.cbrt(SOLID_VOLUME / (PERIODS * ASPECT * (1 - POROSITY)))
EXTENT = SIDE * np.array([PERIODS * ASPECT, 1., 1.])
FLOWS = bypass.FLOWS
VISCOSITY = 1e-3 / 1.2
SEED = 731907
SWEEPS = 2500


def packing():
    points = ordered(PERIODS) * (EXTENT / np.array([PERIODS, 1., 1.]))
    minimum = 1.01 * DIAMETER_M
    rng = np.random.default_rng(SEED)
    accepted = 0
    for _ in range(SWEEPS):
        for i in rng.permutation(len(points)):
            proposed = (points[i] + rng.uniform(-.032 * SIDE, .032 * SIDE, 3)) % EXTENT
            d = offsets(points - proposed, EXTENT)
            distances = np.einsum("ij,ij->i", d, d)
            distances[i] = np.inf
            if distances.min() >= minimum**2:
                accepted += 1
                points[i] = proposed
    if minimum_distance(points, EXTENT) < minimum or accepted < 50 * len(points):
        raise ValueError("Dense spherical packing fails its gap or mixing check")
    return points, accepted / (len(points) * SWEEPS)


def application(source):
    """Derive the rectangular project application from the accepted case."""
    changes = {
        "viscosity = .001": f"viscosity = {VISCOSITY:.17g}",
        "int n, periods;": "int n, periods;\n  int nx() const { return periods * (6*n/5); }",
        "c.periods*c.side-h": "c.nx()*h-h",
        "c.periods*c.side,c.side,c.side": "c.nx()*(c.side/c.n),c.side,c.side",
        "c.periods*c.n": "c.nx()",
        "if (c.n<16 ||": "if (c.n<16 || c.n%5!=0 ||",
    }
    for old, new in changes.items():
        if old not in source:
            raise ValueError(f"Archived application changed at {old}")
        source = source.replace(old, new)
    return source


def build():
    original = PACKING / "openlb"
    source = BUILD / "source"
    source.mkdir(parents=True, exist_ok=True)
    for file in ("Makefile", "config.mk"):
        shutil.copy2(original / file, source / file)
    (source / "capillary3d.cpp").write_text(application((original / "capillary3d.cpp").read_text()))
    return build_case_omp(source, BUILD / "build", BUILD / "build_record")


def prepare():
    old, _, _ = bypass.inputs()
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ") + "-" + uuid4().hex[:8]
    base = OUTPUT / "runs" / stamp
    base.mkdir(parents=True)
    centres = base / "centres_m.txt"
    points, acceptance = packing()
    np.savetxt(centres, points, fmt="%.17g")
    cases = {}
    for flow in FLOWS:
        case = old["cases"][flow + "_l25"]
        p = case["parameters"]
        L = case["bed_cm"] / 100
        cases[flow] = dict(flow_mL_min=case["flow_mL_min"], flow_m3_s=p["flow_m3_s"],
            interstitial_m_s=p["flow_m3_s"] * L / p["liquid_volume_m3"],
            column_length_m=L, column_liquid_m3=p["liquid_volume_m3"],
            dry_carbon_mass_kg=p["dry_mass_kg"], particle_density_kg_m3=p["particle_density_kg_m3"],
            baseline_film_m_s=p["film_m_s"], Pe=case["Pe"], sections=case["sections"],
            access=case["scenarios"]["f50_slow01"])
    write_json(base / "plan.json", dict(side_m=SIDE, extent_m=EXTENT.tolist(), porosity=POROSITY,
        solid_volume_m3=SOLID_VOLUME, centres_sha256=sha(centres), seed=SEED, packing_acceptance=acceptance,
        minimum_sphere_gap_m=minimum_distance(points, EXTENT) - DIAMETER_M,
        viscosity_Pa_s=VISCOSITY, cases=cases, column_sha256=bypass.COLUMN_SHA))
    print(base, flush=True)
    return base


def grid_shape(n):
    if n < 20 or n % 5:
        raise ValueError("Transverse resolution must be divisible by five")
    return PERIODS * (6 * n // 5), n, n


def load_field(base, flow, n, *, warm_trial=False):
    folder = base / f"{flow}_N{n}{'_trial' if warm_trial else ''}"
    raw = read_json(folder / "raw.json")
    if (raw["solver"], raw["collision"], raw["descriptor"], raw["resolution"], raw["converged"]) != (
            "OpenLB 1.9r0", "TRT", "D3Q27", n, True):
        raise ValueError("Expected converged stock OpenLB TRT")
    if not warm_trial and (raw["max_Mach"] > .04 or raw["max_density_departure"] > .002):
        raise ValueError("Native Mach or density screen failed")
    plan = read_json(base / "plan.json")
    h = raw["dx_m"]
    shape = grid_shape(n)
    centres = base / "centres_m.txt"
    if sha(centres) != plan["centres_sha256"]:
        raise ValueError("Frozen sphere centers changed")
    data = np.memmap(folder / "field.f64", mode="r", dtype="<f8").reshape(raw["fluid_nodes"], 7)
    ijk = data[:, :3].astype(np.int32)
    np.testing.assert_array_equal(ijk, np.argwhere(fluid_mask(np.loadtxt(centres), shape, h)))
    u = np.asarray(data[:, 3:6])
    target = plan["cases"][flow]["interstitial_m_s"]
    relative = np.linalg.norm(u.mean(axis=0) - [target, 0., 0.]) / target
    porosity = len(u) / np.prod(shape)
    if not np.isfinite(data).all() or relative > .002 or abs(porosity - POROSITY) > .001:
        raise ValueError("Native field, matched velocity or voxel porosity is invalid")
    info = dict(raw=raw, shape=list(shape), porosity=float(porosity),
        velocity_match_relative=float(relative), liquid_volume_m3=len(u) * h**3,
        directional_flow_m3_s=(u.sum(axis=0) * h**2 / np.array(shape)).tolist(),
        field_sha256=sha(folder / "field.f64"))
    return folder, ijk, u, info


def warm_field(base, flow, old_grid, new_grid, output, *, warm_trial):
    previous, ijk, _, info = load_field(base, flow, old_grid, warm_trial=warm_trial)
    data = np.memmap(previous / "field.f64", mode="r", dtype="<f8").reshape(-1, 7)
    old = np.zeros((*info["shape"], 4))
    old[tuple(ijk.T)] = data[:, 3:]
    old[tuple(ijk.T) + (3,)] -= 1000.
    shape = grid_shape(new_grid)
    fine = np.argwhere(fluid_mask(np.loadtxt(base / "centres_m.txt"), shape, SIDE / new_grid))
    mapped = interpolate(old, shape)[tuple(fine.T)]
    mapped[:, 3] += 1000.
    np.column_stack((fine, mapped)).astype("<f8", copy=False).tofile(output)
    return info["raw"]


def run_case(base, flow, n, threads, warm_grid=None, *, warm_trial=False, pilot=False):
    plan = read_json(base / "plan.json")
    speed = plan["cases"][flow]["interstitial_m_s"]
    grid_shape(n)
    meta = read_json(BUILD / "build_record/build.json")
    binary = Path(meta["binary"])
    if sha(binary) != meta["binary_sha256"]:
        raise ValueError("Stock TRT executable changed")
    if sha(BUILD / "source/capillary3d.cpp") != meta["source_sha256"]["capillary3d.cpp"] or (
        (BUILD / "source/capillary3d.cpp").read_text() !=
        application((PACKING / "openlb/capillary3d.cpp").read_text())
    ):
        raise ValueError("Separate OpenLB application no longer matches the verified source")
    folder = base / f"{flow}_N{n}{'_trial' if pilot else ''}"
    folder.mkdir()
    centres = base / "centres_m.txt"
    # Smaller lattice velocity at high flow keeps Mach under the existing limit.
    tau = .5 + n / 45 * (.027 if flow == FLOWS[0] else .016)
    gradient = (.53e6 * speed, 0., 0.)
    if pilot:
        tau = .552 if flow == FLOWS[0] else .521
        gradient = (1.3e6 * speed, 0., 0.)
    warm = []
    if warm_grid is not None:
        path = folder / "initial.f64"
        old = warm_field(base, flow, warm_grid, n, path, warm_trial=warm_trial)
        gradient = old["gradient_Pa_m"]
        density_scale = (((tau - .5) / n) / ((old["tau_even"] - .5) / warm_grid))**2
        warm = [path, density_scale]
    cmd = [str(binary), *map(str, (n, PERIODS, SIDE, tau, speed, *gradient, centres, *warm))]
    write_json(folder / "request.json", dict(command=cmd, threads=threads,
        plan_sha256=sha(base / "plan.json"), build_sha256=sha(BUILD / "build_record/build.json"),
        warm_field_sha256=sha(warm[0]) if warm else None, warm_grid=warm_grid, warm_trial=warm_trial))
    shutil.copy2(BUILD / "build_record/build.json", folder / "build.json")
    env = {**os.environ, "OMP_NUM_THREADS": str(threads), "OMP_DYNAMIC": "FALSE"}
    print(f"Stock TRT {flow}, N{n} started", flush=True)
    start = perf_counter()
    with (folder / "native.log").open("w") as log:
        proc = subprocess.run(cmd, cwd=folder, env=env, stdout=log, stderr=subprocess.STDOUT, timeout=18000)
    write_json(folder / "execution.json", dict(wall_s=perf_counter() - start, returncode=proc.returncode))
    if proc.returncode:
        raise RuntimeError(f"Stock TRT run failed: {folder / 'native.log'}")
    _, _, _, result = load_field(base, flow, n, warm_trial=pilot)
    write_json(folder / ("pilot_field.json" if pilot else "verified_field.json"), result)
    print(flow, n, "pilot field recorded" if pilot else "native flow verified", result, flush=True)


def run_grid(base, n, threads, warm_grid, flows, warm_trial, pilot):
    with ThreadPoolExecutor(max_workers=2) as pool:
        jobs = [pool.submit(run_case, base, flow, n, threads, warm_grid, warm_trial=warm_trial, pilot=pilot)
                for flow in flows]
        for job in jobs:
            job.result()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path)
    parser.add_argument("--grid", type=int, default=45)
    parser.add_argument("--warm-grid", type=int)
    parser.add_argument("--warm-trial", action="store_true")
    parser.add_argument("--flow", choices=FLOWS)
    parser.add_argument("--pilot", action="store_true")
    parser.add_argument("--threads", type=int, default=5, help="OpenMP threads per flow")
    args = parser.parse_args()
    if args.warm_trial and args.warm_grid is None:
        parser.error("--warm-trial requires --warm-grid")
    if args.pilot and (args.grid != 45 or args.warm_grid is not None):
        parser.error("The pilot starts from rest at N45")
    base = args.base.resolve() if args.base else prepare()
    if not (BUILD / "build_record/build.json").exists():
        build()
    run_grid(base, args.grid, args.threads, args.warm_grid,
             [args.flow] if args.flow else FLOWS, args.warm_trial, args.pilot)


if __name__ == "__main__":
    main()
