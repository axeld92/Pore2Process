import numpy as np
import pytest

from examples.biocatalysis.pentyl_propanoate_transport.phase.run import (
    PROPERTIES, acid_conversion_readout, pair_excess_volume, ternary_split,
)
from pore2process_demo.io import read_json


def test_concentration_readout_and_binary_limit():
    assert acid_conversion_readout(0.4, 300.0, 300.0) == pytest.approx(0.6)
    assert acid_conversion_readout(0.4, 240.0, 300.0) == pytest.approx(0.5)
    assert pair_excess_volume(1.0, 0.0, [-4.09]) == 0
    assert pair_excess_volume(0.5, 0.5, [-4.09]) == pytest.approx(-4.09 / 4)


def test_interpolated_tie_line_closes_phase_balances():
    rows = read_json(PROPERTIES)["ternary_liquid_liquid"]["mole_fraction_endpoints_by_temperature_K"]
    endpoints = np.array(rows["303.15"]).reshape(-1, 2, 3)
    endpoints /= endpoints.sum(axis=2, keepdims=True)
    organic, aqueous = endpoints[3:5].mean(axis=0)
    feed = 2.3 * organic + 0.7 * aqueous
    split = ternary_split(feed, rows["303.15"])
    assert split[0] == pytest.approx(2.3 * organic, abs=1e-12)
    assert split[1] == pytest.approx(0.7 * aqueous, abs=1e-12)
    assert ternary_split(np.array([0.4, 0.6, 0.0]), rows["303.15"]) is None
