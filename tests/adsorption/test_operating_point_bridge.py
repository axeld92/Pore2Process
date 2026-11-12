"""The actual Python-Julia handoff, with an independent dimensional oracle."""
import pytest
from examples.adsorption.column_transport.process import translate


def test_matched_reynolds_point_crosses_real_julia_bridge_in_si():
    p=dict(native_gradient_Pa_m=100.,native_velocity_m_s=.001,
        native_viscosity_Pa_s=.001,native_density_kg_m3=1000.,
        carrier_viscosity_Pa_s=.0016,carrier_density_kg_m3=800.,
        superficial_velocity_m_s=.002,length_m=.5,area_m2=.01,efficiency_fraction=.75)
    r=translate(**p)
    assert all(type(v) is float for v in r.values())
    assert r['pressure_drop_Pa']==pytest.approx(100*.8*2**2*.5)
    assert r['flow_m3_s']==pytest.approx(.002*.01)
    assert r['reference_shaft_power_W']==pytest.approx(160*2e-5/.75)
    with pytest.raises(Exception,match='Reynolds'):
        translate(**{**p,'superficial_velocity_m_s':.003})
