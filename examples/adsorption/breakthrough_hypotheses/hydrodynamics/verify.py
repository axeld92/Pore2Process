"""Check the radial geometry integration and the tracer's two wall materials."""
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from examples.adsorption.dispersion_holdout.run import execute

from . import native
from . import packing as p
from .handoff import julia


def verify_shells():
    points, length, _ = p.create_packing()
    independent = np.random.default_rng(502)
    center = points[independent.integers(len(points), size=600_000)]
    angles = independent.normal(size=center.shape)
    angles /= np.linalg.norm(angles, axis=1, keepdims=True)
    surface = center + angles * p.DIAMETER / 2
    bulk = center + angles * (p.DIAMETER / 2 * np.cbrt(independent.random(len(center))))[:, None]
    rs = np.linalg.norm(surface[:, 1:], axis=1)
    rb = np.linalg.norm(bulk[:, 1:], axis=1)
    for shell in p.shell_geometry(points, order=256):
        cutoff = p.RADIUS - shell["cut_dp"] * p.DIAMETER
        np.testing.assert_allclose(shell["wall_carbon_surface_fraction"], np.mean(rs >= cutoff), rtol=0., atol=.002)
        np.testing.assert_allclose(shell["wall_carbon_fraction"], np.mean(rb >= cutoff), rtol=0., atol=.002)
    pore_volume = np.pi * p.RADIUS**2 * length - len(points) * np.pi * p.DIAMETER**3 / 6
    np.testing.assert_allclose(pore_volume / (np.pi * p.RADIUS**2 * length), p.POROSITY, rtol=0., atol=1e-14)


def verify_wall_materials(path):
    """Use a plug-flow voxel and relabel just one reflecting face as clean carbon."""
    dx = .0001
    speed = .01
    flow = np.array([speed * dx**2, 0., 0.], dtype="<f8")
    for n in ("before.i32", "after.i32"):
        np.zeros(3, dtype="<i4").tofile(path / n)
    for n in ("incoming.f64", "outgoing.f64"):
        flow.tofile(path / n)
    np.array([0, 1, 1], dtype="<i4").tofile(path / "coords.i32")
    np.ones(6, dtype="u1").tofile(path / "carbon_face_in_shell.u8")
    records = {}
    for test, lower_side in (("tube", 3), ("carbon", 2)):
        np.array([4, 4, lower_side, 3, 3, 3], dtype="u1").tofile(path / "material.u8")
        output = path / f"{test}.f64"
        command = [*julia(), str(native.HERE / "track.jl"), str(path), str(output),
                   *map(str, (dx, 1e-8, .0005, .04, 2048, 5654, 3, .000035))]
        execute(command, output.with_suffix(".log"), timeout=120)
        records[test] = np.fromfile(output, dtype="<f8").reshape(2048, 11)
    tube, carbon = (records[name] for name in ("tube", "carbon"))
    np.testing.assert_array_equal(tube[:, :4], carbon[:, :4])
    np.testing.assert_array_equal(tube[:, 2], 1.)
    np.testing.assert_allclose(np.mean(tube[:, 3]), dx / speed, rtol=.05, atol=0.)
    assert np.any(np.isfinite(tube[:, 6])) and not np.any(np.isfinite(tube[:, 4]))
    assert np.any(np.isfinite(carbon[:, 4])) and np.all(carbon[carbon[:, 5] != 0, 5] == 1)


def main():
    verify_shells()
    native.BUILD.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix="verification-", dir=native.BUILD) as path:
        verify_wall_materials(Path(path))
    print("Verified shell fractions, inert plug flow and carbon versus tube wall contact")


if __name__ == "__main__":
    main()
