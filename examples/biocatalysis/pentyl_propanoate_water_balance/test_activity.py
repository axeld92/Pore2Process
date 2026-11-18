import numpy as np
import pytest

from examples.biocatalysis.pentyl_propanoate_film_water.water_activity.run import BatchRate
from examples.biocatalysis.pentyl_propanoate_transport.phase.run import PROPERTIES
from examples.biocatalysis.pentyl_propanoate_water_balance.activity import (
    TernaryProjection, binary_activity, liquid_moles, measured_pentanol_curve,
)
from examples.biocatalysis.pentyl_propanoate_water_balance.run import (
    INPUT, conditional_batch_rate_ratio, gas_water_at_saturation, water_saturation_pressure,
)
from pore2process_demo.io import read_json


def test_stoichiometric_water_at_equal_conversion():
    for conversion in [0.0, 0.1, 0.5, 0.75]:
        reference = liquid_moles(conversion)
        assert reference.sum() == pytest.approx(3)
        assert reference[2] == reference[3]
        assert reference[2] / reference.sum() == pytest.approx(conversion / 3)
        assert liquid_moles(conversion, 0.05) - reference == pytest.approx([0, 0, 0.05, 0])
    with pytest.raises(ValueError):
        liquid_moles(0.1, -0.2)


def test_binary_reference_preserves_reported_activities_and_range():
    vle = read_json(INPUT)["pentanol_water_vle"]["observations"]
    curve = measured_pentanol_curve(vle)
    for _, pentanol, _, gamma, _ in vle:
        assert binary_activity(np.array([0, pentanol, 1 - pentanol, 0]), curve) == pytest.approx(
            (1 - pentanol) * gamma
        )
    assert binary_activity(liquid_moles(0.0), curve) == 0
    assert binary_activity(liquid_moles(0.7), curve) is None


def test_same_tie_line_gives_same_activity_for_different_phase_amounts():
    tie_lines = read_json(PROPERTIES)["ternary_liquid_liquid"]["mole_fraction_endpoints_by_temperature_K"]
    curve = measured_pentanol_curve(read_json(INPUT)["pentanol_water_vle"]["observations"])
    for rows in tie_lines.values():
        projection = TernaryProjection.from_tie_lines(rows)
        organic, aqueous = projection.endpoints[4]
        for shape in (None, curve):
            at_binodal = projection.state(np.r_[organic, 0], False, shape)
            assert at_binodal == pytest.approx((aqueous[2], False), abs=1e-10)
            for fraction in (0.1, 0.4, 0.7):
                mixture = np.r_[(1 - fraction) * organic + fraction * aqueous, 0]
                two_phase = projection.state(mixture, False, shape)
                assert two_phase == pytest.approx((aqueous[2], True), abs=1e-10)


def test_headspace_against_water_steam_table_and_volume_scaling():
    data = read_json(INPUT)["pure_water_vapor_pressure"]
    assert water_saturation_pressure(353.15, data) == pytest.approx(47414, abs=10)
    mol = gas_water_at_saturation(353.15, 10, data)
    assert mol == pytest.approx(1.615e-4, rel=1e-3)
    assert gas_water_at_saturation(353.15, 12, data) == pytest.approx(1.2 * mol)
    assert gas_water_at_saturation(353.15, 0, data) == 0


def test_analogue_rate_uses_relative_shift_inside_measured_range():
    series = [BatchRate(0.12, 1.5), BatchRate(0.33, 1.08), BatchRate(0.53, 0.5)]
    assert conditional_batch_rate_ratio(0.53, 0.33, series) == pytest.approx(0.5 / 1.08)
    assert conditional_batch_rate_ratio(0.43, 0.33, series) == pytest.approx(np.sqrt(0.5 / 1.08))
    assert conditional_batch_rate_ratio(0.45, 0.45, series) == 1
    assert conditional_batch_rate_ratio(0.65, 0.33, series) is None
