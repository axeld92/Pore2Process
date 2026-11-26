"""Conserve tube throughput and summarize signed flow by distance to the wall."""
import argparse
from pathlib import Path

import numpy as np

from examples.adsorption.disordered_packing import mesh
from examples.adsorption.dispersion_holdout.run import execute
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT

from . import native
from . import packing as p


def julia(threads=1):
    config = read_json(ROOT / ".local/julia_runtime.json")
    return [config["executable"], f"--threads={threads}", "--startup-file=no", f"--project={config['project']}"]


def shell_weights(n, dx, cut_dp, order=12):
    """Fraction of each axial square in the shell, using product quadrature."""
    coordinate = (np.arange(n) + .5 - n / 2) * dx
    edge2 = (p.RADIUS - cut_dp * p.DIAMETER)**2
    weight = np.zeros((n, n))
    offsets = ((np.arange(order) + .5) / order - .5) * dx
    for dy in offsets:
        for dz in offsets:
            weight += (coordinate[:, None] + dy)**2 + (coordinate[None, :] + dz)**2 >= edge2
    return weight / order**2


def connected_pore_component(ijk, material, velocity):
    """Return the connected face graph containing the fastest axial pore.

    OpenLB also streams along diagonals. Tiny corner-only pockets consequently
    count in the native porosity but are outside the six-face tracer graph.
    """
    pore = material == 1
    component = np.zeros_like(pore)
    component[tuple(ijk[np.argmax(velocity[:, 0])])] = True
    while True:
        next_component = component.copy()
        for axis in range(3):
            next_component |= np.roll(component, 1, axis=axis) | np.roll(component, -1, axis=axis)
        next_component &= pore
        if np.array_equal(next_component, component):
            break
        component = next_component
    return component[tuple(ijk.T)]


def connected_pores(ijk, material, velocity):
    """Apply the adsorption example's exclusion limit to its tracer mesh."""
    active = connected_pore_component(ijk, material, velocity)
    if np.count_nonzero(~active) > .001 * len(ijk):
        raise ValueError("More than 0.1% of native void lacks face connectivity to the tracer mesh")
    return active


def conservation(ijk, info, faces, q):
    """Apply the project's conservation and local-correction limits separately."""
    left, right, axis = (faces[name] - 1 for name in ("left.i32", "right.i32", "axis.u8"))
    raw = faces["unbalanced.f64"]
    if q.shape != raw.shape or not np.isfinite(q).all():
        raise ValueError("Invalid balanced shared faces")
    net = np.bincount(left, weights=q, minlength=len(ijk)) - np.bincount(right, weights=q, minlength=len(ijk))
    local = float(np.max(np.abs(net)) / np.max(np.abs(q)))
    throughput = info["directional_flow_m3_s"][0]
    plane = 0.
    for dimension, width in enumerate(info["shape"]):
        use = axis == dimension
        sums = np.bincount(ijk[left[use], dimension], weights=q[use], minlength=width)
        plane = max(plane, float(np.max(np.abs(sums - (throughput if dimension == 0 else 0.))) / throughput))
    correction = float(np.linalg.norm(q - raw) / np.linalg.norm(raw))
    raw_axial = float(np.sum(raw[axis == 0]) / (info["shape"][0] * throughput))
    rescaled = raw / raw_axial
    shape_correction = float(np.linalg.norm(q - rescaled) / np.linalg.norm(rescaled))
    if local > 1e-10 or plane > 1e-9:
        raise ValueError(f"The flow handoff fails local or plane conservation: {local:.4g}, {plane:.4g}")
    return dict(max_cell_imbalance_relative=local, max_plane_error_relative=plane,
        face_rms_correction_relative=correction, raw_Q_over_native_Q=raw_axial,
        face_rms_after_global_rescale_relative=shape_correction, qualifies_local_correction=correction <= .05)


def radial_summary(plan, ijk, material, info, faces, balanced):
    dims = tuple(info["shape"])
    dx = info["raw"]["dx_m"]
    axial = faces["axis.u8"] == 1
    locations = ijk[faces["left.i32"][axial] - 1]
    sx, sy, sz = locations.T
    q = balanced[axial]
    raw_q = faces["unbalanced.f64"][axial]
    throughput = info["directional_flow_m3_s"][0]
    area_map = np.count_nonzero(material != 3, axis=0)
    fluid_map = np.count_nonzero(material == 1, axis=0)
    positive = np.bincount(sx, weights=np.maximum(q, 0), minlength=dims[0])
    back = np.bincount(sx, weights=np.maximum(-q, 0), minlength=dims[0])
    signed = positive - back
    np.testing.assert_allclose(signed, throughput, rtol=1e-9, atol=0)
    shells = []
    for geometric in plan["shells"]:
        weights = shell_weights(dims[1], dx, geometric["cut_dp"])
        axial_shell = np.bincount(sx, weights=weights[sy, sz] * q, minlength=dims[0]) / signed
        raw_wall = np.dot(weights[sy, sz], raw_q) / raw_q.sum()
        solid = geometric["wall_carbon_fraction"] * (1 - p.POROSITY)
        wall_void_fraction = (geometric["wall_area_fraction"] - solid) / p.POROSITY
        nearwall_void = float(np.sum(weights * fluid_map) / fluid_map.sum())
        nearwall_area = float(np.sum(weights * area_map) / area_map.sum())
        shells.append(dict(**geometric, wall_void_fraction=wall_void_fraction,
            wall_flow_fraction=float(axial_shell.mean()), raw_wall_flow_fraction=float(raw_wall),
            wall_flow_fraction_by_plane=axial_shell.tolist(),
            wall_flow_fraction_plane_range=[float(axial_shell.min()), float(axial_shell.max())],
            voxel_wall_void_fraction=nearwall_void, voxel_wall_geometric_fraction=nearwall_area,
            shell_continuous_porosity=p.POROSITY * wall_void_fraction / geometric["wall_area_fraction"]))
    return dict(shells=shells, radial_definition="outer shell, (tube radius minus radial position)/sphere diameter",
        flow_definition="signed axial face flow averaged over all planes; square faces weighted by shell overlap",
        throughput_m3_s=throughput, cut_backflow_fraction=(back / positive).tolist(),
        cut_forward_throughput_m3_s=positive.tolist(), porosity_voxel=info["porosity"],
        cells_per_dp=info["cells_per_dp"], period_length_m=plan["length_m"],
        field_sha256=info["field_sha256"], geometry_plan_sha256=info["geometry_plan_sha256"])


def project(base, cases):
    plan, _ = native.geometry(base)
    for case in cases:
        folder = base / case
        request = read_json(folder / "request.json")
        if request["plan_sha256"] != sha(base / "plan.json"):
            raise ValueError("Native case belongs to a different geometry")
        ijk, velocity, material, info = native.validate(base, folder, request["rate"], request["mesh"])
        info["geometry_plan_sha256"] = request["plan_sha256"]
        active = connected_pores(ijk, material, velocity)
        inactive_count = int(np.count_nonzero(~active))
        info["corner_only_void_count"] = inactive_count
        info["corner_only_void_fraction"] = inactive_count / len(ijk)
        info["corner_only_native_axial_flow_fraction"] = float(np.sum(velocity[~active, 0]) / np.sum(velocity[:, 0]))
        info["native_liquid_volume_m3"] = info["liquid_volume_m3"]
        info["liquid_volume_m3"] *= float(active.mean())
        ijk, velocity = ijk[active], velocity[active]
        path = folder / "faces"
        path.mkdir(exist_ok=True)
        if (path / "summary.json").exists():
            _, _, _, existing = mesh.load_faces(folder)
            if existing["info"] != info:
                raise ValueError("Saved face projection belongs to a different field")
            continue
        faces = mesh.shared_faces(ijk, velocity, info["shape"], info["raw"]["dx_m"])
        for filename, values in faces.items():
            values.tofile(path / filename)
        ijk.astype("<i4", copy=False).tofile(path / "coords.i32")
        total = info["shape"][0] * info["directional_flow_m3_s"][0]
        cmd = [*julia(), str(native.HERE / "balance.jl"), str(path), str(len(ijk)), str(total)]
        print(f"Projecting and checking {case} tube face flows", flush=True)
        elapsed = execute(cmd, path / "balance.log", timeout=1800)
        q = np.fromfile(path / "balanced.f64", dtype="<f8")
        checks = conservation(ijk, info, faces, q)
        save = dict(info=info, checks=checks, elapsed=elapsed, command=cmd,
            files={name: sha(path / name) for name in (*mesh.FLOW_FILES, "coords.i32", "balanced.f64")})
        radial = radial_summary(plan, ijk, material, info, faces, q)
        radial["face_checks"] = checks
        radial["corner_only_void_count"] = inactive_count
        radial["corner_only_void_fraction"] = info["corner_only_void_fraction"]
        radial["rate"] = request["rate"]
        radial["condition"] = plan["flows"][request["rate"]]
        write_json(folder / "radial.json", radial)
        write_json(path / "summary.json", save)
        print(f"{case}: half-diameter wall flow {radial['shells'][1]['wall_flow_fraction']:.2%}, "
              f"void {radial['shells'][1]['wall_void_fraction']:.2%}, "
              f"carbon {radial['shells'][1]['wall_carbon_fraction']:.2%}; local face correction "
              f"{checks['face_rms_correction_relative']:.2%}", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("base", type=Path)
    parser.add_argument("cases", nargs="+")
    args = parser.parse_args()
    project(args.base.resolve(), args.cases)


if __name__ == "__main__":
    main()
