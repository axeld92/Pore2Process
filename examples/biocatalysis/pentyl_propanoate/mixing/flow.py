"""Run stock TRT water flow in a periodic sphere packing, then balance face rates."""
import argparse
import math
import os
from pathlib import Path
import shutil
import subprocess

import numpy as np

from examples.adsorption.breakthrough_hypotheses.hydrodynamics.handoff import (
    connected_pores,
    conservation,
)
from examples.adsorption.disordered_packing.mesh import FLOW_FILES, fluid_mask, shared_faces
from examples.adsorption.disordered_packing.packing import (
    DIAMETER_M as REFERENCE_DIAMETER_M,
    FCC,
    minimum_distance,
    offsets,
)
from examples.adsorption.transport_root_cause.hydraulics.stock_omp import build_case_omp
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT

HERE = Path(__file__).resolve().parent
SOURCE = ROOT / "examples/adsorption/disordered_packing"
BUILD = ROOT / ".local/biocatalysis/pentyl_propanoate/mixing"
OUTPUT = ROOT / "results/biocatalysis/pentyl_propanoate/mixing"
DIAMETER_M = 593e-6
POROSITY = 0.64
COLUMN_DIAMETER_M = 0.016
COLUMN_LENGTH_M = 0.14
WATER_FLOW_M3_S = 1.2e-6 / 60
PULSE_MIN = 8.4 / 1.2
SUPERFICIAL_M_S = WATER_FLOW_M3_S / (math.pi * (COLUMN_DIAMETER_M / 2) ** 2)
WATER_VISCOSITY_PA_S = 0.00089
WATER_DENSITY_KG_M3 = 997.05
PERIODS = 4
BEAD_COUNT = 4 * 2 * 2 * (2 * PERIODS)
SIDE_M = (BEAD_COUNT * math.pi * DIAMETER_M**3 / (6 * PERIODS * (1 - POROSITY))) ** (1 / 3)
GRIDS = (32, 48, 64)
SEED = 420619


def packing():
    """Keep the same continuous geometry during mesh refinement."""
    extent = np.array((2 * PERIODS, 2, 2), dtype=float)
    points = np.concatenate([FCC + cell for cell in np.ndindex(*(extent.astype(int)))])
    clearance = 1.025 * DIAMETER_M / (SIDE_M / 2)
    rng = np.random.default_rng(SEED)
    accepted = 0
    sweeps = 1200
    for _ in range(sweeps):
        for i in rng.permutation(len(points)):
            candidate = (points[i] + rng.uniform(-0.15, 0.15, 3)) % extent
            distance = np.linalg.norm(offsets(points - candidate, extent), axis=1)
            distance[i] = np.inf
            if distance.min() >= clearance:
                points[i] = candidate
                accepted += 1
    minimum = minimum_distance(points, extent)
    if minimum < clearance or len(points) != BEAD_COUNT:
        raise RuntimeError("The periodic bead inventory or minimum clearance changed")
    return points * (SIDE_M / 2), minimum * (SIDE_M / 2), accepted / (len(points) * sweeps)


def application():
    """Stage the existing OpenLB application with this experiment's bead and water properties."""
    case = (SOURCE / "openlb/capillary3d.cpp").read_text()
    replacements = {
        "diameter = .0008, viscosity = .001, density = 1000.": (
            f"diameter = {DIAMETER_M:.17g}, viscosity = {WATER_VISCOSITY_PA_S:.17g}, "
            f"density = {WATER_DENSITY_KG_M3:.17g}"
        ),
        "std::size_t(4*c.periods)": f"std::size_t({BEAD_COUNT})",
    }
    for old, new in replacements.items():
        if case.count(old) != 1:
            raise RuntimeError(f"OpenLB application changed at {old}")
        case = case.replace(old, new)
    return case


def build():
    source = BUILD / "source"
    source.mkdir(parents=True, exist_ok=True)
    for filename in ("Makefile", "config.mk"):
        shutil.copyfile(SOURCE / "openlb" / filename, source / filename)
    (source / "capillary3d.cpp").write_text(application())
    return build_case_omp(source, BUILD / "build", BUILD / "record")


def prepare(output):
    output.mkdir(parents=True, exist_ok=True)
    centres, minimum, acceptance = packing()
    centre_file = output / "centres_m.txt"
    np.savetxt(centre_file, centres, fmt="%.17g")
    write_json(output / "plan.json", dict(
        reference_doi="10.1021/acs.iecr.5c02780",
        reference_source="Sections 2.1–2.2 and Appendix A",
        centres_sha256=sha(centre_file),
        seed=SEED,
        beads=BEAD_COUNT,
        bead_diameter_m=DIAMETER_M,
        periodic_extent_m=[PERIODS * SIDE_M, SIDE_M, SIDE_M],
        target_void_fraction=POROSITY,
        minimum_bead_gap_m=minimum - DIAMETER_M,
        packing_move_acceptance_fraction=acceptance,
        column_diameter_m=COLUMN_DIAMETER_M,
        column_length_m=COLUMN_LENGTH_M,
        column_flow_m3_s=WATER_FLOW_M3_S,
        superficial_speed_m_s=SUPERFICIAL_M_S,
        target_interstitial_m_s=SUPERFICIAL_M_S / POROSITY,
        injection_duration_min=PULSE_MIN,
        temperature_assumption_C=25,
        water_viscosity_Pa_s=WATER_VISCOSITY_PA_S,
        water_density_kg_m3=WATER_DENSITY_KG_M3,
        particle_reynolds_superficial=(
            WATER_DENSITY_KG_M3 * SUPERFICIAL_M_S * DIAMETER_M / WATER_VISCOSITY_PA_S
        ),
    ))
    return centre_file


def solve(output, centre_file, grid, build_record, threads):
    """Run the existing native convergence and throughput checks."""
    case = output / f"N{grid}"
    case.mkdir(exist_ok=True)
    spacing = SIDE_M / grid
    shape = (PERIODS * grid, grid, grid)
    # The earlier rasterizer works in its original diameter units.
    scale = REFERENCE_DIAMETER_M / DIAMETER_M
    geometry = fluid_mask(np.loadtxt(centre_file) * scale, shape, spacing * scale)
    voxel_porosity = float(geometry.mean())
    velocity = SUPERFICIAL_M_S / voxel_porosity
    gradient = (
        150 * WATER_VISCOSITY_PA_S * SUPERFICIAL_M_S * (1 - voxel_porosity)**2
        / (DIAMETER_M**2 * voxel_porosity**3)
    )
    if abs(voxel_porosity - POROSITY) > 0.005:
        raise RuntimeError(f"Voxel void fraction {voxel_porosity:.5f} missed the packing target")
    command = [
        build_record["binary"], *map(str, (grid, PERIODS, SIDE_M, 0.72, velocity,
                                        gradient, 0.0, 0.0, centre_file)),
    ]
    env = {**os.environ, "OMP_NUM_THREADS": str(threads), "OMP_DYNAMIC": "FALSE"}
    write_json(case / "native_request.json", dict(command=command, threads=threads))
    with (case / "native.log").open("w") as log:
        print(f"Running OpenLB TRT N{grid}", flush=True)
        subprocess.run(command, cwd=case, env=env, stdout=log, stderr=subprocess.STDOUT,
                       check=True, timeout=2700)
    native = read_json(case / "raw.json")
    if (native["solver"], native["collision"], native["converged"]) != ("OpenLB 1.9r0", "TRT", True):
        raise RuntimeError("Expected a converged native OpenLB 1.9r0 TRT field")
    if native["max_Mach"] > 0.04 or native["max_density_departure"] > 0.002:
        raise RuntimeError("Native field exceeds Mach or density limits")
    field = np.fromfile(case / "field.f64", dtype="<f8").reshape(-1, 7)
    if not np.isfinite(field).all() or len(field) != native["fluid_nodes"]:
        raise RuntimeError("Native field is nonfinite or incomplete")
    cells = field[:, :3].astype(np.int32)
    np.testing.assert_array_equal(cells, field[:, :3])
    np.testing.assert_array_equal(cells, np.argwhere(geometry))
    target_flow = SUPERFICIAL_M_S * SIDE_M**2
    native_flow = field[:, 3].sum() * spacing**2 / shape[0]
    np.testing.assert_allclose(native_flow, target_flow, rtol=0.002, atol=0)
    active = connected_pores(cells, geometry, field[:, 3:6])
    return case, cells[active], field[active, 3:6], dict(
        native=native,
        grid=grid,
        shape=list(shape),
        cells_per_bead=DIAMETER_M / spacing,
        voxel_void_fraction=voxel_porosity,
        target_periodic_flow_m3_s=target_flow,
        raw_periodic_flow_m3_s=float(native_flow),
        directional_flow_m3_s=[float(native_flow), 0., 0.],
        raw_throughput_relative_error=float(native_flow / target_flow - 1),
        raw_velocity_mean_m_s=field[:, 3:6].mean(axis=0).tolist(),
        corner_only_pore_cells=int(np.count_nonzero(~active)),
        corner_only_pore_volume_fraction=float(np.mean(~active)),
        corner_only_native_axial_flow_fraction=float(field[~active, 3].sum() / field[:, 3].sum()),
        tracked_void_fraction=float(active.sum() / np.prod(shape)),
        field_sha256=sha(case / "field.f64"),
    )


def balance(case, cells, velocity, record, julia):
    shape = record["shape"]
    dx = record["native"]["dx_m"]
    folder = case / "faces"
    folder.mkdir(exist_ok=True)
    faces = shared_faces(cells, velocity, shape, dx)
    for filename, array in faces.items():
        array.tofile(folder / filename)
    volume_flow = record["raw_periodic_flow_m3_s"]
    sums = (shape[0] * volume_flow, 0.0, 0.0)
    command = [julia, "--startup-file=no", "--project=julia", str(SOURCE / "balance.jl"),
               str(folder), str(len(cells)), *map(str, sums)]
    print(f"Projecting N{record['grid']} face rates at fixed throughput", flush=True)
    with (folder / "balance.log").open("w") as log:
        subprocess.run(command, cwd=ROOT, check=True, stdout=log, stderr=subprocess.STDOUT, timeout=900)
    rates = np.fromfile(folder / "balanced.f64", dtype="<f8")
    left, right, axis = (faces[key] - 1 for key in ("left.i32", "right.i32", "axis.u8"))
    checks = conservation(cells, record, faces, rates)
    axial = axis == 0
    plane_reflux = np.bincount(
        cells[right[axial], 0], weights=np.maximum(0.0, -rates[axial]), minlength=shape[0]
    )
    # Write a periodic, one-based compact flow graph for the existing Julia tracker.
    before = np.zeros(cells.shape, dtype="<i4")
    after = np.zeros_like(before)
    incoming = np.zeros(cells.shape, dtype="<f8")
    outgoing = np.zeros_like(incoming)
    before[right, axis] = left + 1
    after[left, axis] = right + 1
    incoming[right, axis] = rates
    outgoing[left, axis] = rates
    compact = case / "compact"
    compact.mkdir(exist_ok=True)
    for filename, array in (
        ("before.i32", before), ("after.i32", after), ("incoming.f64", incoming),
        ("outgoing.f64", outgoing), ("x_index.i32", cells[:, 0].astype("<i4")),
    ):
        array.tofile(compact / filename)
    record.update(
        face_conservation=checks,
        axial_plane_reflux_fractions=(plane_reflux / volume_flow).tolist(),
        balanced_files_sha256={
            name: sha(folder / name) for name in (*FLOW_FILES, "balanced.f64")
        },
        compact_files_sha256={path.name: sha(path) for path in sorted(compact.iterdir())},
    )
    write_json(case / "field_summary.json", record)
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--julia", default="julia")
    parser.add_argument("--threads", type=int, default=6)
    parser.add_argument("--grids", nargs="+", type=int, choices=GRIDS, default=GRIDS)
    args = parser.parse_args()
    if args.threads < 1:
        parser.error("--threads must be positive")
    output = args.output.resolve()
    centres = prepare(output)
    record_path = BUILD / "record/build.json"
    if record_path.is_file():
        record = read_json(record_path)
        binary = Path(record["binary"])
        if (not binary.is_file() or record["binary_sha256"] != sha(binary)
                or (BUILD / "source/capillary3d.cpp").read_text() != application()):
            record = build()
    else:
        record = build()
    shutil.copy2(record_path, output / "build.json")
    for grid in args.grids:
        result = solve(output, centres, grid, record, args.threads)
        balance(*result, args.julia)


if __name__ == "__main__":
    main()
