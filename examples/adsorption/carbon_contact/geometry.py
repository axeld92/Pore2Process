"""Audit the solid surface and local flow around the archived periodic FCC cell."""
import itertools
import json

import numpy as np

from examples.adsorption.exchange_column.closure import EXCHANGE_SHA, MARKERS
from examples.adsorption.porous_face_handoff.run import archived
from pore2process_demo.io import read_json, sha

FCC_CENTRES = np.array([[0, 0, 0], [0, .5, .5], [.5, 0, .5], [.5, .5, 0]])
RADIUS_M = .0004
DIRECTIONS = tuple(itertools.product(range(3), (-1, 1)))
OCTANT_BITS = np.array([4, 2, 1])


def _sphere_distances(points, box_m):
    """Use minimum-image displacements from the four periodic grain centres."""
    offset = points[:, None, :] / box_m - FCC_CENTRES[None, :, :]
    offset -= np.rint(offset)
    return offset, np.linalg.norm(offset, axis=2) * box_m


def _patch_ids(grain, radial):
    return 8*grain + (radial >= 0) @ OCTANT_BITS


def solid_patch_counts(mask, h):
    """Count solid voxels in the 32 fixed grain octants."""
    solid_coords = np.column_stack(np.nonzero(~mask))
    box_m = len(mask) * h
    radial, distance = _sphere_distances((solid_coords+.5) * h, box_m)
    grain = distance.argmin(axis=1)
    chosen = np.arange(len(solid_coords)), grain
    if np.any(distance[chosen] > RADIUS_M+1e-14):
        raise ValueError("Archived solid mask extends beyond the FCC sphere indicators")
    counts = np.bincount(_patch_ids(grain, radial[chosen]), minlength=32)
    if counts.sum() != np.count_nonzero(~mask):
        raise ValueError("Solid patch counts do not preserve the geometry's solid inventory")
    return counts


def wall_faces(info, ijk):
    """Return solid faces keyed by zero-based archived fluid row, axis and side."""
    n, h = info["n"], info["dx_m"]
    box_m = n * h
    mask = np.zeros((n, n, n), dtype=bool)
    mask[tuple(ijk.T)] = True
    owner_parts, axis_parts, side_parts = [], [], []
    grain_parts, patch_parts, area_parts = [], [], []
    for axis, side in DIRECTIONS:
        neighbor = ijk.copy()
        neighbor[:, axis] = (neighbor[:, axis] + side) % n
        ids = np.flatnonzero(~mask[tuple(neighbor.T)])
        _, distance = _sphere_distances((neighbor[ids]+.5) * h, box_m)
        grain = distance.argmin(axis=1)
        centre = ijk[ids] + .5
        centre[:, axis] += side / 2
        radial, _ = _sphere_distances(centre * h, box_m)
        normal = radial[np.arange(len(ids)), grain]
        normal /= np.linalg.norm(normal, axis=1, keepdims=True)
        # Each Cartesian face supplies projected area. Dividing by the
        # L1 normal factor removes the orientation bias of staircase area.
        area_parts.append(h**2 / np.abs(normal).sum(axis=1))
        owner_parts.append(ids)
        axis_parts.append(np.full(len(ids), axis, dtype=np.int8))
        side_parts.append(np.full(len(ids), side, dtype=np.int8))
        grain_parts.append(grain)
        patch_parts.append(_patch_ids(grain, normal))
    _, liquid_distance = _sphere_distances((ijk+.5) * h, box_m)
    if np.any(np.min(liquid_distance, axis=1) < RADIUS_M-1e-14):
        raise ValueError("Archived liquid mask extends into an FCC sphere")
    owner, axis, side = [np.concatenate(x) for x in (owner_parts, axis_parts, side_parts)]
    face_id = 6*owner+2*axis+np.where(side < 0, 1, 2)
    if len(np.unique(face_id)) != len(face_id):
        raise ValueError("Wall face IDs must be unique")
    solid_counts = solid_patch_counts(mask, h)
    return {
        "face_id": face_id,
        "owner": owner,
        "axis": axis,
        "side": side,
        "grain": np.concatenate(grain_parts),
        "patch": np.concatenate(patch_parts),
        "smooth_area_m2": np.concatenate(area_parts),
        "solid_voxel_count_by_patch": solid_counts,
        "solid_volume_by_grain_m3": h**3 * solid_counts.reshape(4, 8).sum(axis=1),
    }


def _balanced_centres(info, ijk, faces):
    """Read the tracer's actual two-period geometry and balance at the wall."""
    metadata = MARKERS / "summary.json"
    if sha(metadata) != EXCHANGE_SHA:
        raise ValueError("The pinned path-exchange summary changed")
    record = read_json(metadata)
    names = {
        "before.i32": "<i4", "after.i32": "<i4",
        "incoming.f64": "<f8", "outgoing.f64": "<f8",
    }
    arrays = {}
    count = len(ijk)
    for name, dtype in names.items():
        relative = "compact/" + name
        path = MARKERS / relative
        if sha(path) != record["artifacts_sha256"][relative]:
            raise ValueError("The saved tracer geometry or flow changed")
        arrays[name] = np.fromfile(path, dtype=dtype).reshape(2, count, 3)
    owner, axis, side = (faces[k] for k in ("owner", "axis", "side"))
    for sign, neighbor_name, flow_name in (
        (-1, "before.i32", "incoming.f64"), (1, "after.i32", "outgoing.f64")
    ):
        selected = side == sign
        np.testing.assert_array_equal(arrays[neighbor_name][:, owner[selected], axis[selected]], 0)
        np.testing.assert_array_equal(arrays[flow_name][:, owner[selected], axis[selected]], 0.)
        # Compact topology can also omit an open fluid face at the axial ends.
        # A zero neighbor at zero rate is always an actual wall.
        actual = int(np.count_nonzero((arrays[neighbor_name] == 0) & (arrays[flow_name] == 0.)))
        if actual != 2*int(selected.sum()):
            raise ValueError("The solid mask and compact wall counts differ")
    np.testing.assert_array_equal(arrays["incoming.f64"][0], arrays["incoming.f64"][1])
    np.testing.assert_array_equal(arrays["outgoing.f64"][0], arrays["outgoing.f64"][1])
    return (arrays["incoming.f64"][0] + arrays["outgoing.f64"][0]) / (2 * info["dx_m"]**2)


def _wall_speed_stats(velocity, faces, mean_u_m_s):
    speed = np.linalg.norm(velocity, axis=1)
    order = np.argsort(speed, kind="stable")
    deciles = np.empty(len(speed), dtype=np.int8)
    for i, ids in enumerate(np.array_split(order, 10)):
        deciles[ids] = i
    owner = faces["owner"]
    areas = faces["smooth_area_m2"]
    by_speed = np.bincount(deciles[owner], weights=areas, minlength=10)
    by_speed_stair = np.bincount(deciles[owner], minlength=10)
    order_wall = np.argsort(speed[owner], kind="stable")
    cum_area = np.cumsum(areas[order_wall])
    median_speed = np.interp(cum_area[-1] / 2, cum_area, speed[owner][order_wall])
    return dict(
        local_pore_cell_speed_deciles_m_s=np.quantile(speed, np.arange(.1, 1, .1)).tolist(),
        wall_area_fraction_by_pore_cell_speed_decile=(by_speed / by_speed.sum()).tolist(),
        staircase_area_fraction_by_pore_cell_speed_decile=(by_speed_stair / by_speed_stair.sum()).tolist(),
        mean_wall_owner_speed_m_s=float(np.dot(areas, speed[owner]) / areas.sum()),
        median_wall_owner_speed_m_s=float(median_speed),
        median_wall_owner_speed_relative_to_interstitial=float(median_speed / mean_u_m_s),
        role="Pore-cell-centre speeds in volume deciles. Surface area has local speed labels only.",
    )


def _record(info, ijk, velocity, faces):
    by_grain = np.bincount(faces["grain"], weights=faces["smooth_area_m2"], minlength=4)
    area = float(by_grain.sum())
    solid = faces["solid_volume_by_grain_m3"]
    h = info["dx_m"]
    interstitial = info["carrier_directional_flow_m3_s"][0] / ((info["n"] * h)**2 * info["porosity"])
    exact = float(4 * 4 * np.pi * RADIUS_M**2)
    return dict(
        source=info,
        grain_centres_relative=FCC_CENTRES.tolist(),
        solid_face_count=int(len(faces["owner"])),
        wall_adjacent_liquid_voxels=int(len(np.unique(faces["owner"]))),
        wall_adjacent_liquid_fraction=float(len(np.unique(faces["owner"])) / len(ijk)),
        staircase_surface_area_m2=float(len(faces["owner"]) * h**2),
        normal_corrected_surface_area_m2=area,
        analytic_fcc_surface_area_m2=exact,
        normal_corrected_area_relative_error=area / exact - 1.,
        normal_corrected_area_fraction_by_grain=(by_grain / area).tolist(),
        normal_corrected_area_by_patch_m2=np.bincount(faces["patch"], weights=faces["smooth_area_m2"], minlength=32).tolist(),
        solid_voxel_count_by_patch=faces["solid_voxel_count_by_patch"].tolist(),
        solid_voxel_count=int(faces["solid_voxel_count_by_patch"].sum()),
        solid_volume_by_grain_m3=solid.tolist(),
        solid_volume_fraction_by_grain=(solid / solid.sum()).tolist(),
        mean_interstitial_velocity_m_s=interstitial,
        native_cell_speed=_wall_speed_stats(velocity, faces, interstitial),
    )


def audit():
    """Reconstruct the geometry and inspect wall-adjacent cell speeds."""
    info, ijk, velocity = archived(96)
    faces = wall_faces(info, ijk)
    n96 = _record(info, ijk, velocity, faces)
    n96["balanced_wall_and_flux_verified"] = True
    n96["balanced_cell_speed"] = _wall_speed_stats(
        _balanced_centres(info, ijk, faces), faces, n96["mean_interstitial_velocity_m_s"]
    )
    info, ijk, velocity = archived(64)
    n64 = _record(info, ijk, velocity, wall_faces(info, ijk))
    return dict(
        primary_n96=n96,
        geometry_sensitivity_n64=n64,
        compact_reference_sha256=EXCHANGE_SHA,
        patch_id_note="Grain ID times 8 plus the octant bits (4, 2, 1) for nonnegative x, y, z offsets.",
        velocity_note="Each velocity is sampled at the wall-adjacent liquid centre, not at the wall.",
        area_note="The orientation correction is based on smooth radial normals. Areas are geometric external areas.",
        interpretation="Speed classes at one location do not identify fluid histories or transit-time classes.",
    )


if __name__ == "__main__":
    print(json.dumps(audit(), indent=2, allow_nan=False))
