"""Frozen input boundary, genuine Julia arrays, and independent exported balances."""
from collections import Counter

import numpy as np
import pytest

from examples.biocatalysis.pentyl_propanoate_history import aggregate,run
from pore2process_demo.io import read_json


def test_frozen_calibration_and_complete_matrix(monkeypatch):
    def forbidden():
        raise AssertionError("Reactor observations must not enter prediction preparation.")
    monkeypatch.setattr(run.eq,"read_observations",forbidden)
    plan=read_json(run.PLAN); jl=run.load_runtime(); tasks,laws=run.prepare(jl,plan)
    assert len(tasks)==560 and len(laws)==7
    assert {c["actual_mass_kg"] for c in tasks}=={.006}
    assert set(Counter((c["scenario"],c["law"],c["feed_ratio"]) for c in tasks).values())=={20}
    frozen=read_json(run.eq.FROZEN)["batch"]["all_data"]
    for law,p in laws.items():
        actual=np.array(jl.CompositionKinetics.parameter_vector(p))
        expected=[frozen["fit"][law]["parameters"][k] for k in run.eq.prior.PARAMETERS]
        assert np.array_equal(actual,expected)
    assert frozen["equilibrium_training"]["1"]["calibration_temperatures_C"]==[40.,50.,60.]
    assert all(c["alcohol0_mol_L"]==c["feed_ratio"]*c["acid0_mol_L"] for c in tasks)
    assert all(c["actual_mass_kg"]==.006 for c in tasks)
    assert sum(run.conditions_flags(c)["recovered"] for c in tasks)==80
    assert sum(run.conditions_flags(c)["surface_hazard"] for c in tasks)==120
    assert sum(run.conditions_flags(c)["equilibrium_prefill"] for c in tasks)==80


def test_real_julia_transient_refinement_and_fixed_catalyst_mass():
    plan=read_json(run.PLAN); jl=run.load_runtime(); tasks,laws=run.prepare(jl,plan)
    for feed in (1,2):
        c=next(c for c in tasks if (c["scenario"],c["law"],c["feed_ratio"],c["temperature_C"],c["flow_ml_min"])==
               ("fresh_bulk","protected_residual_decay",feed,60,3))
        r,d=run.solve_refined(jl,laws[c["law"]],c,plan)
        n=d["fine_cells"]
        assert r["activity"].shape==(n,12)
        assert r["bulk_progress"].dtype==np.float64
        assert np.array_equal(r["times_min"],plan["times_min"])
        assert np.all(np.diff(r["cumulative_hazard"],axis=1)>=0)
        assert np.all(r["activity"][:,-1]<=1.)
        assert np.isclose(r["catalyst_mass_per_cell_kg"]*n,.006,rtol=1e-15)
        assert r["initial_liquid_product_mol"]==0.
        assert np.allclose(r["generated_mol"],r["liquid_product_mol"]+r["backdiffused_product_mol"]+
                           r["effluent_product_mol"],rtol=1e-8,atol=1e-10)
        for tag in ("time","grid"):
            assert run.passes({k:d[f"{tag}_{k}"] for k in run.LIMIT_KEYS},plan["numerics"])
        ss=run.snapshots(r,c)
        assert ss[-1]["time_min"]==180 and ss[-1]["actual_mass_kg"]==.006
        assert ss[-1]["last_quarter_activity"]>=ss[-1]["first_quarter_activity"]-1e-10


def test_initial_liquid_inventory_is_not_counted_as_new_generation():
    plan=read_json(run.PLAN); jl=run.load_runtime(); tasks,laws=run.prepare(jl,plan)
    c=next(c for c in tasks if (c["scenario"],c["law"],c["feed_ratio"],c["temperature_C"],c["flow_ml_min"])==
           ("equilibrium_prefill_bulk","protected_residual_decay",1,40,3))
    r=run.one_history(jl,laws[c["law"]],c,plan,64,1.)
    assert r["initial_liquid_product_mol"]==pytest.approx(c["acid0_mol_L"]*c["void_volume_L"]*r["equilibrium_conversion"])
    assert r["generated_mol"][0]==0. and r["effluent_product_mol"][0]==0.
    assert r["liquid_product_mol"][0]==r["initial_liquid_product_mol"]
    assert np.allclose(r["generated_mol"],r["liquid_product_mol"]-r["initial_liquid_product_mol"]+
                       r["backdiffused_product_mol"]+r["effluent_product_mol"],atol=1e-10,rtol=1e-8)


def test_stability_trace_is_not_assigned_to_an_operating_condition():
    trace=aggregate.trace_reference()
    assert trace["temperature_C"] is None and trace["flow_ml_min"] is None
    old=read_json(run.eq.FROZEN)["independent_bed_time_trace"]
    assert trace["span_30_to_180_pp"]==pytest.approx(old["observed_30_to_180min_span_pp"])
    assert trace["reading_half_range_pp"]==pytest.approx(old["reading_half_range_pp"])
    assert trace["times_min"]==[10.,30.,60.,90.,120.,150.,180.]


def test_refinement_volume_averages_fine_cells():
    # Fine-cell local profile is deliberately different, but volume average
    # on the same axial support is identical to the coarse data.
    t=np.array([0.,1.,10.,30.]); a=np.array([[.8,.6,.5,.4],[.9,.8,.7,.6]])
    b=np.repeat(a,2,axis=0); b[::2]-=.01; b[1::2]+=.01
    c={"acid0_mol_L":5.3,"flow_ml_min":.7}
    x=np.array([0.,.2,.4,.6]); led={k:np.zeros(4) for k in run.LEDGERS}
    coarse={"times_min":t,"activity":a,"outlet_conversion":x,**led}
    fine={"times_min":t,"activity":b,"outlet_conversion":x,**led}
    differences=run.refinement(coarse,fine,c)
    assert max(differences.values())<1e-15
