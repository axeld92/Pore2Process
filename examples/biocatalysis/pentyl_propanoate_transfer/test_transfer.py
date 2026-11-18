import numpy as np
import pytest

from examples.biocatalysis.pentyl_propanoate_transfer.run import (
    FIGURE_S3, PROTOCOL, at_exposure, make_predictions,
)
from examples.biocatalysis.pentyl_propanoate.run import read_rows
from pore2process_demo.io import read_json


def test_exposure_matching_interpolates_and_never_extrapolates():
    series = np.array([0.0, 20.0, 40.0]), np.array([0.0, 0.5, 0.6])
    assert at_exposure(0.001, 0.006, 0.2, series) == pytest.approx((30, 0.55))
    assert at_exposure(0.001, 0.006, 0.15, series) == pytest.approx((40, 0.6))
    assert at_exposure(0.001, 0.006, 0.1, series)[1] is None


def test_source_markers_and_coverage():
    markers = [row for row in read_rows(FIGURE_S3) if row["marker_x_px_190dpi"]]
    assert len(markers) == 68
    for row in markers:
        x, y = (int(row[f"marker_{axis}_px_190dpi"]) for axis in ("x", "y"))
        assert float(row["time_min"]) == pytest.approx(600 * (x - 875.5) / 438, abs=0.05)
        assert float(row["acid_conversion"]) == pytest.approx((1376 - y) / 438, abs=0.0005)

    rows, _ = make_predictions(read_json(PROTOCOL))
    assert len(rows) == 20
    assert all(row["catalyst_mass_kg"] == 0.006 for row in rows)
    assert sum(row["predicted_from_2025_batch"] is not None for row in rows) == 18
    assert sum(row["predicted_from_2026_batch"] is not None for row in rows) == 15
