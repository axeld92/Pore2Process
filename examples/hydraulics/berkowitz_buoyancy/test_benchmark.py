import numpy as np
import pytest

from .report import comparison, moments


def test_finite_pulse_linearized_and_covarying():
    entry=np.linspace(.5,299.5,300)
    transit=8000+0.12*entry+70*np.cos(np.arange(len(entry)))
    result=moments(entry,entry+transit,.2)
    var=transit.var(ddof=1)+2*np.cov(entry,transit,ddof=1)[0,1]
    assert result["alpha_mm"]==pytest.approx(100*var/transit.mean()**2)
    assert result["covariance_entry_residence_s2"]>0
    delayed=moments(entry+1000,entry+1000+transit,.2)
    assert delayed["alpha_mm"]==pytest.approx(result["alpha_mm"])


def test_matched_comparison_requires_same_field_and_inventory():
    def case(grav,alpha):
        return dict(geometry_sha256="bed",water_sha256="water",markers=20,seed=1,
                    liquid_Q_mL_min=.3,pulse_mass_mg=.75,gravity_multiplier=grav,coupling_s=60,native_sweeps=96,tracer_dt_s=4,
                    moments=[{"alpha_mm":alpha}],one_hour=dict(centroid_mm=-grav,fraction_below=.5+.1*grav))
    cases={"on":case(1,4),"off":case(0,2)}
    pair=comparison(cases,"on","off",10)
    assert pair["gap_closed_fraction"]==pytest.approx(.25)
    assert pair["centroid_shift_at_one_hour_mm"]==-1
    assert not pair["passes_published_comparison"]
    cases["on"]["water_sha256"]="other"
    with pytest.raises(ValueError,match="paired input"):
        comparison(cases,"on","off",10)
