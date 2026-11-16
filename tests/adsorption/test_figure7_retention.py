import numpy as np
import pytest

from examples.adsorption.figure7_retention.run import checked_model, integrate, model_snapshot


def test_integral_uses_partial_segment_and_stops_at_observations():
    t = np.array([0.,10.,20.])
    y = np.array([0.,.5,1.])
    assert integrate(t,y,1.,15.) == 9.375
    assert integrate(t,y,1.,20.) == 10.
    with pytest.raises(ValueError,match="inside the digitized"):
        integrate(t,y,1.,21.)


def test_column_balance_uses_interval_mean_outlet_and_separate_inventories():
    p = dict(flow_m3_s=.5,inlet_concentration_kg_m3=2.,dry_mass_kg=3.,liquid_volume_m3=.2,
             capacity_kg_kg=.9,affinity_m3_kg=.2)
    result = dict(time_s=[0.,2.,3.],outlet_concentration_kg_m3=[0.,.5,1.],
                  liquid_inventory_kg=[0.,.1,.2],adsorbed_inventory_kg=[0.,1.4,1.8],
                  cumulative_outlet_kg=[0.,.5,1.])
    m = checked_model(dict(inputs=p,result=result),dict(p,feed_rate_kg_s=1.))
    s = model_snapshot(m,2.5,1.)
    assert s["retained_kg"] == 1.75
    assert s["net_retention_rate_kg_s"] == .5
    assert s["adsorbed_kg"] + s["liquid_kg"] == pytest.approx(1.75)
    with pytest.raises(AssertionError):
        checked_model(dict(inputs=p,result=dict(result,cumulative_outlet_kg=[0.,.4,1.])),
                      dict(p,feed_rate_kg_s=1.))
