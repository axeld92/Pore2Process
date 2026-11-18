"""Generate and solve a wall-bounded stock TRT packing without column fitting."""
import argparse
from datetime import datetime, timezone
import os
from pathlib import Path
import shutil
import subprocess
from time import perf_counter
from uuid import uuid4

import numpy as np

from examples.adsorption.disordered_packing.native import HERE as PERIODIC
from examples.adsorption.transport_root_cause.hydraulics.stock_omp import build_case_omp
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT
from pore2process_demo.runtime.openlb import verify_vendor

from . import packing as p

HERE = Path(__file__).resolve().parent
OUTPUT = ROOT / "results/adsorption/breakthrough_hypotheses/hydrodynamics"
BUILD = ROOT / ".local/adsorption/breakthrough_hypotheses/hydrodynamics"
APP = HERE / "openlb/capillary3d.cpp"


def build():
    stage = BUILD / "source"
    stage.mkdir(parents=True, exist_ok=True)
    for name in ("Makefile", "config.mk"):
        shutil.copy2(PERIODIC / "openlb" / name, stage / name)
    shutil.copy2(APP, stage / APP.name)
    return build_case_omp(stage, BUILD / "build", BUILD / "build_record")


def prepare():
    tag = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ") + "-" + uuid4().hex[:8]
    base = OUTPUT / "runs" / tag
    base.mkdir(parents=True)
    points, length, packing = p.create_packing()
    np.savetxt(base / "centres_m.txt", points, fmt="%.17g")
    cross_section = np.pi * p.RADIUS**2
    flows = {}
    for name, rate in p.RATES.items():
        q = rate * 1e-6 / 60
        superficial = q / cross_section
        pore_speed = superficial / p.POROSITY
        flows[name] = dict(flow_mL_min=rate, flow_m3_s=q, superficial_m_s=superficial,
            nominal_interstitial_m_s=pore_speed, Re_particle_superficial=p.DENSITY * superficial * p.DIAMETER / p.VISCOSITY,
            Pe_particle_interstitial=pore_speed * p.DIAMETER / p.DIFFUSIVITY,
            nominal_cell_residence_s=length / pore_speed)
    solid_volume = len(points) * np.pi * p.DIAMETER**3 / 6
    geometric_volume = cross_section * length
    np.testing.assert_allclose(1 - solid_volume / geometric_volume, p.POROSITY, rtol=0, atol=2e-15)
    write_json(base / "plan.json", dict(study="Alhamed, doi:10.1016/j.jhazmat.2009.05.002, Figure 5 conditions",
        geometry="short synthetic, axial periodicity only, no-slip nonadsorbing cylindrical outer wall",
        diameter_m=p.DIAMETER, radius_m=p.RADIUS, solid_spheres=len(points), length_m=length,
        full_column_length_m=.25, area_m2=cross_section, porosity=p.POROSITY,
        geometric_solid_m3=solid_volume, geometric_volume_m3=geometric_volume,
        fixed_printed_bed_density_kg_m3=740., full_column_dry_carbon_kg=740. * cross_section * .25,
        density_kg_m3=p.DENSITY, viscosity_Pa_s=p.VISCOSITY, diffusivity_m2_s=p.DIFFUSIVITY,
        flows=flows, packing=packing, shells=p.shell_geometry(points, order=256),
        centres_sha256=sha(base / "centres_m.txt")))
    print(base, flush=True)
    return base


def geometry(base):
    plan = read_json(base / "plan.json")
    centres = base / "centres_m.txt"
    if sha(centres) != plan["centres_sha256"]:
        raise ValueError("Frozen sphere centres changed")
    return plan, np.loadtxt(centres)


def field(folder):
    raw = read_json(folder / "raw.json")
    if (raw["solver"], raw["collision"], raw["descriptor"], raw["boundary"], raw["converged"]) != (
            "OpenLB 1.9r0", "TRT", "D3Q27", "stock BounceBack", True):
        raise ValueError("Expected converged stock OpenLB TRT with stock no-slip boundaries")
    data = np.memmap(folder / "field.f64", mode="r", dtype="<f8").reshape(raw["fluid_nodes"], 7)
    if not np.isfinite(data).all():
        raise ValueError("Nonfinite native field")
    return raw, data


def validate(base, folder, rate, mesh):
    plan, points = geometry(base)
    raw, data = field(folder)
    shape, h = p.grid(plan["length_m"], mesh)
    if tuple(raw["shape"]) != shape or not np.isclose(raw["dx_m"], h, rtol=1e-14, atol=0.):
        raise ValueError("Native field grid differs from frozen geometry")
    material = p.voxel_material(points, shape, h)
    ijk = data[:, :3].astype(np.int32)
    np.testing.assert_array_equal(ijk, data[:, :3])
    np.testing.assert_array_equal(ijk, np.argwhere(material == 1))
    velocity = np.asarray(data[:, 3:6])
    q = float(velocity[:, 0].sum() * h**2 / shape[0])
    expected = plan["flows"][rate]["flow_m3_s"]
    porosity = len(ijk) / np.count_nonzero(material != 3)
    relative = abs(q - expected) / expected
    if (raw["max_Mach"] > .04 or raw["max_density_departure"] > .002 or
            relative > .002 or abs(porosity - p.POROSITY) > .005):
        raise ValueError(f"Native gate failed: Mach {raw['max_Mach']:.5g}, density "
                         f"{raw['max_density_departure']:.5g}, flux error {relative:.5g}, porosity {porosity:.5g}")
    info = dict(raw=raw, shape=list(shape), porosity=porosity, liquid_volume_m3=len(ijk) * h**3,
        directional_flow_m3_s=[q, 0., 0.], velocity_match_relative=relative,
        mean_native_velocity_m_s=velocity.mean(axis=0).tolist(), field_sha256=sha(folder / "field.f64"),
        cells_per_dp=p.DIAMETER / h, voxel_solid_volume_m3=float(np.count_nonzero(material == 2) * h**3))
    return ijk, velocity, material, info


def initial_gradient(speed):
    eps = p.POROSITY
    return (150 * p.VISCOSITY * (1 - eps)**2 * speed / (p.DIAMETER**2 * eps**3) +
            1.75 * p.DENSITY * (1 - eps) * speed**2 / (p.DIAMETER * eps**3))


def warm_start(old_folder, folder, shape, h, points, speed_ratio, tau):
    raw, data = field(old_folder)
    old = np.zeros((*raw["shape"], 4))
    indexes = tuple(data[:, :3].astype(np.int32).T)
    old[indexes] = data[:, 3:]
    old[indexes + (3,)] -= p.DENSITY
    if tuple(raw["shape"]) != shape:
        # Both rasters use the same physical axis but padding changes with h.
        # Map across the real coordinates rather than across the padded box.
        ijk = np.argwhere(p.voxel_material(points, shape, h) == 1)
        units = (ijk + .5) * h / raw["dx_m"] - .5
        units[:, 1:] += (raw["shape"][1] - shape[1] * h / raw["dx_m"]) / 2
        lower = np.floor(units).astype(int)
        fraction = units - lower
        mapped = np.zeros((len(ijk), 4))
        for offset in np.ndindex(2, 2, 2):
            index = (lower + offset) % np.asarray(raw["shape"])
            weight = np.prod(np.where(offset, fraction, 1 - fraction), axis=1)
            mapped += weight[:, None] * old[tuple(index.T)]
    else:
        ijk = data[:, :3]
        mapped = old[indexes]
    mapped[:, :3] *= speed_ratio
    mapped[:, 3] = p.DENSITY + mapped[:, 3] * speed_ratio
    path = folder / "initial.f64"
    np.column_stack((ijk, mapped)).astype("<f8", copy=False).tofile(path)
    density_scale = ((tau - .5) * h / ((raw["tau_even"] - .5) * raw["dx_m"]))**2
    return [path, density_scale]


def run(base, rate, mesh, threads, tag, tau_delta, magic, warm):
    plan, points = geometry(base)
    shape, h = p.grid(plan["length_m"], mesh)
    folder = base / f"{rate}_N{mesh}{tag}"
    folder.mkdir()
    build_file = BUILD / "build_record/build.json"
    build_record = read_json(build_file)
    binary = Path(build_record["binary"])
    if sha(APP) != build_record["source_sha256"][APP.name] or sha(binary) != build_record["binary_sha256"]:
        raise ValueError("Native source or binary changed; rerun --build")
    shutil.copy2(build_file, folder / "build.json")
    centre_file = folder / "centres_native_m.txt"
    np.savetxt(centre_file, points + (0, shape[1] * h / 2, shape[1] * h / 2), fmt="%.17g")
    conditions = plan["flows"][rate]
    superficial = conditions["superficial_m_s"]
    gradient = initial_gradient(superficial)
    tau = .5 + tau_delta if tau_delta is not None else .5 + (.002 if rate == "q625" else .0048) * (mesh / 8)
    if magic is None:
        magic = .04 if rate == "q625" else .1875
    initialization = []
    if warm is not None:
        old_folder = base / warm
        old_record = read_json(old_folder / "request.json")
        if old_record["plan_sha256"] != sha(base / "plan.json"):
            raise ValueError("Warm flow belongs to a different geometry")
        ratio = conditions["flow_m3_s"] / plan["flows"][old_record["rate"]]["flow_m3_s"]
        initialization = warm_start(old_folder, folder, shape, h, points, ratio, tau)
        old_gradient = read_json(old_folder / "raw.json")["gradient_Pa_m"][0]
        gradient = old_gradient * gradient / initial_gradient(superficial / ratio)
    command = [str(binary), *map(str, (shape[0], shape[1], plan["length_m"], p.RADIUS, tau, magic,
                                     conditions["flow_m3_s"], gradient, centre_file, *initialization))]
    write_json(folder / "request.json", dict(command=command, rate=rate, mesh=mesh, threads=threads,
        plan_sha256=sha(base / "plan.json"), warm=warm, source_sha256=sha(APP), centre_sha256=sha(centre_file)))
    start = perf_counter()
    env = {**os.environ, "OMP_NUM_THREADS": str(threads), "OMP_DYNAMIC": "FALSE"}
    print(f"Launching stock OpenLB TRT {rate} {shape}, {threads} threads in {folder}", flush=True)
    with (folder / "native.log").open("w") as log:
        proc = subprocess.run(command, env=env, cwd=folder, stdout=log, stderr=subprocess.STDOUT,
                              timeout=10800, check=False)
    write_json(folder / "execution.json", dict(elapsed_s=perf_counter() - start, returncode=proc.returncode))
    if proc.returncode:
        raise RuntimeError(f"Native flow failed; see {folder / 'native.log'}")
    verify_vendor()
    _, _, _, info = validate(base, folder, rate, mesh)
    write_json(folder / "verified_field.json", info)
    print(f"Verified stock OpenLB TRT {rate} N{mesh}: {info}", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--prepare", action="store_true")
    mode.add_argument("--build", action="store_true")
    mode.add_argument("--run", choices=p.RATES)
    parser.add_argument("--base", type=Path)
    parser.add_argument("--mesh", type=int, default=8)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--tau-delta", type=float)
    parser.add_argument("--magic", type=float, help="Stock TRT odd/even product; default .1875 at q233 and .04 at q625")
    parser.add_argument("--warm", help="Earlier case folder name within the same run")
    parser.add_argument("--tag", default="", help="Append a suffix for separate convergence checks")
    args = parser.parse_args()
    if args.prepare:
        prepare()
    elif args.build:
        print(build()["binary"])
    elif args.base and args.threads > 0 and args.tag == os.path.basename(args.tag):
        run(args.base.resolve(), args.run, args.mesh, args.threads, args.tag, args.tau_delta, args.magic, args.warm)
    else:
        parser.error("--run needs --base, a positive --threads and a suffix without a path")


if __name__ == "__main__":
    main()
