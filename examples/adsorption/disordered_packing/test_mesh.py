import numpy as np

from .mesh import fluid_mask, shared_faces
from .packing import DIAMETER_M, PERIODS, SIDE_M, describe, disordered, ordered


def test_disorder_preserves_continuous_inventory():
    centres, acceptance = disordered()
    original = describe(ordered(PERIODS), PERIODS)
    changed = describe(centres, PERIODS)
    assert changed["particle_count"] == original["particle_count"] == 32
    assert changed["solid_volume_m3"] == original["solid_volume_m3"]
    assert changed["minimum_center_distance_m"] > DIAMETER_M
    assert .01 < acceptance < 1
    mask = fluid_mask(centres*SIDE_M, (PERIODS*32,32,32), SIDE_M/32)
    assert abs(mask.mean()-.6) < .003
    assert not np.array_equal(mask[:32], mask[32:64])


def test_rectangular_periodic_face_rates():
    shape = (8,3,3)
    coords = np.argwhere(np.ones(shape, dtype=bool))
    velocity = np.broadcast_to([2.,0.,0.],(len(coords),3))
    faces = shared_faces(coords,velocity,shape,.1)
    axial = faces["axis.u8"]==1
    a = coords[faces["left.i32"][axial]-1]
    b = coords[faces["right.i32"][axial]-1]
    np.testing.assert_array_equal((a[:,0]+1)%shape[0],b[:,0])
    np.testing.assert_array_equal(a[:,1:],b[:,1:])
    np.testing.assert_allclose(faces["unbalanced.f64"][axial],.02)
