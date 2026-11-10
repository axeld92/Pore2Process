"""Actual Julia integration for the filtration boundary; analytical independent oracle."""
import pytest
from examples.filtration.layered_support.model import layer_load
from pore2process_demo.bridge import runtime
from examples.filtration.layered_support.run import SETTINGS


def test_actual_julia_layer_units_and_fixed_selective_pair():
    settings={"support_thickness_m":1e-3,"flux_m_s":100e-3/3600,
        "membrane_area_m2":2.,"pump_efficiency_fraction":.75,"other_fixed_resistance_m_inv":1e10}
    base=layer_load(1e-13,settings,2e10,1e-3)
    better=layer_load(2e-13,settings,2e10,1e-3)
    # Rs=1e10, Rl=2e10, Rfixed=1e10; better Rs=0.5e10.
    assert base["total_pressure_drop_Pa"]==pytest.approx(100e-3/3600*1e-3*4e10)
    assert base["flow_m3_s"]==pytest.approx(200e-3/3600)
    assert better["shaft_power_W"]/base["shaft_power_W"]==pytest.approx(3.5/4)
    assert better["selective_resistance_m_inv"]==base["selective_resistance_m_inv"]
    assert base["selective_dominates"] is False
    assert better["selective_dominates"] is True
    assert base["specific_shaft_energy_kWh_m3"]==pytest.approx(base["total_pressure_drop_Pa"]/.75/3.6e6)
    assert all(type(value) in (float,bool) for value in base.values())


def test_permeance_conversion_and_zero_flux():
    jl=runtime()
    # 360 LMH/bar = 1e-9 m/(Pa s), not 1e-4 or 1e-6.
    assert float(jl.permeance_LMH_bar_to_SI(360.))==pytest.approx(1e-9)
    assert float(jl.resistance_from_permeance(1e-9,1e-3))==pytest.approx(1e12)
    settings={"support_thickness_m":1e-3,"flux_m_s":0.,"membrane_area_m2":1.,
        "pump_efficiency_fraction":.75,"other_fixed_resistance_m_inv":0.}
    result=layer_load(1e-13,settings,1e12,1e-3)
    assert result["total_pressure_drop_Pa"]==result["shaft_power_W"]==result["flow_m3_s"]==0.
    assert SETTINGS["flux_m_s"]==pytest.approx(2.777777777777778e-5)
    assert SETTINGS["flux_LMH"]==100.
