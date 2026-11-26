"""Deterministic, short synthetic cylinder at the Alhamed bed void fraction."""
import numpy as np

DIAMETER = .0008
RADIUS = .0055
POROSITY = .427
VISCOSITY = .001 / 1.2
DENSITY = 1000.
DIFFUSIVITY = 1e-9
MIN_GAP_DP = .004
MIN_TUBE_GAP_DP = .001
SEED = 842707
SWEEPS = 320
CUT_DP = (.25, .5, .75, 1., 1.5)
RATES = dict(q233=23.3, q625=62.5)


def minimum_separation(centres, length):
    difference = centres[:, None] - centres[None, :]
    difference[..., 0] -= length * np.rint(difference[..., 0] / length)
    distance2 = np.sum(difference**2, axis=-1)
    np.fill_diagonal(distance2, np.inf)
    return np.sqrt(distance2.min())


def initial_centres(radius):
    basis = np.array([[0, 0, 0], [0, .5, .5], [.5, 0, .5], [.5, .5, 0]])
    translations = np.stack(np.meshgrid(np.arange(3), np.arange(-5, 6),
                                         np.arange(-5, 6), indexing="ij"), axis=-1)
    candidates = (translations.reshape(-1, 1, 3) + basis).reshape(-1, 3)
    for pitch in (1.435, 1.45, 1.465, 1.48, 1.495, 1.42):
        for offset in ((.17, .29), (.2, 0.), (.07, .07), (.34, .15), (0., 0.)):
            points = candidates * pitch
            points[:, 1:] += offset
            points = points[np.linalg.norm(points[:, 1:], axis=1) <= radius - .5 - MIN_TUBE_GAP_DP]
            length = len(points) / (6 * (1 - POROSITY) * radius**2)
            points[:, 0] *= length / (3 * pitch)
            if minimum_separation(points, length) >= 1 + MIN_GAP_DP:
                return points, length, pitch, offset
    raise ValueError("Could not construct the declared nonoverlapping packing")


def create_packing():
    radius = RADIUS / DIAMETER
    points, length, pitch, offset = initial_centres(radius)
    original = points.copy()
    accepted = 0
    rng = np.random.default_rng(SEED)
    for _ in range(SWEEPS):
        for index in rng.permutation(len(points)):
            proposal = points[index] + rng.uniform(-.09, .09, 3)
            proposal[0] %= length
            if np.hypot(*proposal[1:]) > radius - .5 - MIN_TUBE_GAP_DP:
                continue
            delta = points - proposal
            delta[:, 0] -= length * np.rint(delta[:, 0] / length)
            distance = np.einsum("ij,ij->i", delta, delta)
            distance[index] = np.inf
            if distance.min() >= (1 + MIN_GAP_DP)**2:
                points[index] = proposal
                accepted += 1
    separation = minimum_separation(points, length)
    wall_gap = radius - .5 - np.max(np.linalg.norm(points[:, 1:], axis=1))
    displacement = points - original
    displacement[:, 0] -= length * np.rint(displacement[:, 0] / length)
    if separation < 1 + MIN_GAP_DP or wall_gap < MIN_TUBE_GAP_DP or accepted < 10 * len(points):
        raise ValueError("Synthetic packing fails geometric validity or displacement check")
    # Shift the arbitrary axial origin out of the original FCC plane.
    points[:, 0] = (points[:, 0] + .137 * length) % length
    return points * DIAMETER, length * DIAMETER, dict(seed=SEED, moves_per_sphere=SWEEPS,
        accepted_move_fraction=accepted / (SWEEPS * len(points)), minimum_sphere_gap_dp=separation - 1,
        minimum_wall_gap_dp=wall_gap, rms_displacement_from_template_dp=float(np.sqrt(
            np.mean(np.sum(displacement**2, axis=1)))), initial_pitch_dp=pitch, initial_offset_dp=list(offset))


def grid(length, cells_per_dp):
    if cells_per_dp < 5:
        raise ValueError("Use at least five requested cells per particle diameter")
    nx = round(length / DIAMETER * cells_per_dp)
    h = length / nx
    nt = int(np.ceil(2 * RADIUS / h)) + 4
    return (nx, nt, nt), h


def voxel_material(points, shape, h, *, radius=RADIUS, diameter=DIAMETER):
    """Independent cell-centre raster for checking the stock material IDs."""
    ny = shape[1]
    transverse = (np.arange(ny) + .5 - ny / 2) * h
    radius2 = transverse[:, None]**2 + transverse[None, :]**2
    material = np.where(radius2 <= radius**2, 1, 3).astype(np.uint8)
    material = np.broadcast_to(material, shape).copy()
    centre = ny * h / 2
    shifted = points + (0, centre, centre)
    for point in shifted:
        low = np.floor((point - diameter / 2) / h - .5).astype(int)
        high = np.ceil((point + diameter / 2) / h - .5).astype(int)
        axes = [np.arange(a, b + 1) for a, b in zip(low, high, strict=True)]
        mesh = np.meshgrid(*axes, indexing="ij", sparse=True)
        distance2 = sum(((index + .5) * h - value)**2 for index, value in zip(mesh, point, strict=True))
        selection = np.ix_(axes[0] % shape[0], *axes[1:])
        local = material[selection]
        solid = distance2 <= (diameter / 2)**2
        if np.any(local[solid] != 1):
            raise ValueError("A grain raster overlaps another grain or the solid tube wall")
        local[solid] = 2
        material[selection] = local
    return material


def shell_geometry(points, distances=CUT_DP, order=96):
    """Integrate the true sphere surface and solid volume in each radial shell.

    For each axial coordinate relative to a sphere centre, angular integration
    gives the arc share and circle-intersection geometry gives the volume share.
    """
    radius = RADIUS / DIAMETER
    offsets = np.linalg.norm(points[:, 1:], axis=1) / DIAMETER
    abscissa, weights = np.polynomial.legendre.leggauss(order)
    rings = np.sqrt(1 - abscissa**2) / 2
    centre = offsets[:, None]
    arc_product = 2 * centre * rings
    volume_total = len(points) * np.pi / 6
    results = []
    for thickness in distances:
        boundary = radius - thickness
        ratio = np.divide(boundary**2 - centre**2 - rings**2, arc_product,
                          out=np.zeros_like(arc_product), where=arc_product > 0)
        arc = np.arccos(np.clip(ratio, -1, 1)) / np.pi
        arc[centre[:, 0] == 0] = rings > boundary
        outside = centre - rings >= boundary
        inside = centre + rings <= boundary
        contains = rings - centre >= boundary
        ring_area = np.pi * rings**2
        ratio2 = np.divide(boundary**2 + centre**2 - rings**2, 2 * boundary * centre,
                           out=np.zeros_like(arc_product), where=centre > 0)
        minor = np.arccos(np.clip(ratio2, -1, 1))
        radicand = ((centre + rings)**2 - boundary**2) * (boundary**2 - (centre - rings)**2)
        overlap = rings**2 * np.pi * (1 - arc) + boundary**2 * minor - .5 * np.sqrt(np.maximum(radicand, 0))
        overlap = np.where(outside, 0., np.where(inside, ring_area,
                           np.where(contains, np.pi * boundary**2, overlap)))
        results.append(dict(cut_dp=thickness, wall_area_fraction=float(1 - (boundary / radius)**2),
            wall_carbon_fraction=float(np.sum((ring_area - overlap) @ weights) / (2 * volume_total)),
            wall_carbon_surface_fraction=float(np.mean(arc @ weights) / 2)))
    return results
