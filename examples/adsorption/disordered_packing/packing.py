"""Build a seeded hard-sphere packing and its exact-inventory FCC control."""
import numpy as np

DIAMETER_M = .0008
POROSITY = .6
PERIODS = 8
SIDE_M = DIAMETER_M*np.cbrt(2*np.pi/(3*(1-POROSITY)))
SEED = 731907
SWEEPS = 2500
MINIMUM_GAP_DIAMETERS = .025
FCC = np.array([[0., 0., 0.], [0., .5, .5], [.5, 0., .5], [.5, .5, 0.]])


def ordered(periods):
    return np.concatenate([FCC+[i, 0, 0] for i in range(periods)])


def offsets(delta, extent):
    return delta-extent*np.rint(delta/extent)


def minimum_distance(points, extent):
    delta = offsets(points[:, None, :]-points[None, :, :], extent)
    i, j = np.triu_indices(len(points), 1)
    return float(np.min(np.linalg.norm(delta[i, j], axis=1)))


def disordered():
    points = ordered(PERIODS)
    extent = np.array([PERIODS, 1., 1.])
    threshold = (DIAMETER_M*(1+MINIMUM_GAP_DIAMETERS)/SIDE_M)**2
    rng = np.random.default_rng(SEED)
    accepted = 0
    for _ in range(SWEEPS):
        for i in rng.permutation(len(points)):
            proposed = (points[i]+rng.uniform(-.16, .16, 3)) % extent
            delta = offsets(points-proposed, extent)
            distance2 = np.einsum("ij,ij->i", delta, delta)
            distance2[i] = np.inf
            if np.min(distance2) >= threshold:
                points[i] = proposed
                accepted += 1
    if minimum_distance(points, extent)**2 < threshold or accepted < len(points)*50:
        raise ValueError("Hard-sphere packing did not mix or preserve the specified gap")
    return points, accepted/(len(points)*SWEEPS)


def describe(points, periods):
    extent = np.array([periods, 1., 1.])
    if len(points) != periods*4:
        raise ValueError("Sphere count does not match the carbon inventory")
    volume = len(points)*np.pi*DIAMETER_M**3/6
    return dict(particle_count=len(points), diameter_m=DIAMETER_M, extent_m=(SIDE_M*extent).tolist(),
                solid_volume_m3=volume, analytical_porosity=1-volume/(periods*SIDE_M**3),
                minimum_center_distance_m=SIDE_M*minimum_distance(points, extent))
