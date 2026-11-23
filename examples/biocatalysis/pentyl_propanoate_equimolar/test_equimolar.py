"""Provenance, frozen calibration boundary and real Julia integration."""
from math import isclose

import numpy as np

from examples.biocatalysis.pentyl_propanoate_equimolar import run as audit
from examples.biocatalysis.pentyl_propanoate_equimolar import source_sign_audit as sign
from examples.biocatalysis.pentyl_propanoate.run import read_rows
from pore2process_demo.io import read_json


def test_all_twenty_new_points_trace_to_both_figures():
    meta = read_json(audit.EXTRA)
    digit = meta["digitization"]
    rows = read_rows(audit.OBSERVATIONS[1])
    assert len(rows) == 20
    assert set(range(45,65)) == {int(r["source_run_id"]) for r in rows}
    for r in rows:
        assert int(r["source_run_id"]) == 45+4*((int(r["temperature_C"])-40)//10)+(.7,1.,2.,3.).index(float(r["flow_ml_min"]))
        y0,y100 = digit["si_y0_y100_px"][r["figure"]]
        expected = (y0-float(r["marker_y_px"]))/(y0-y100)
        assert float(r["acid_conversion"]) == round(expected,3)
        y0,y100 = digit["main_figure_3b_y0_y100_px"]
        other = (y0-float(r["main_figure_3b_marker_y_px"]))/(y0-y100)
        assert float(r["main_figure_3b_conversion"]) == round(other,4)
        # One overlapped red marker in the main figure is visibly partial;
        # the more legible S1 open marker defines the held-out observation.
        assert abs(expected-other) < (.008 if int(r["source_run_id"]) == 51 else .005)
    assert len(audit.read_observations()) == 40
    assert meta["reaction_1to1"]["individual_conversion_uncertainty"] is None
    eta = meta["effectiveness"]["by_temperature_C"]
    assert [eta[str(t)]["eta_internal"] for t in (40,50,60,70,80)] == [.893,.903,.943,.990,.998]


def test_frozen_case_specs_do_not_depend_on_bed_scores():
    frozen = read_json(audit.FROZEN)
    specs = audit.case_specs(frozen)
    assert specs == audit.case_specs({"batch":frozen["batch"],"bed_comparison":{"irrelevant":123}})
    assert len(specs) == 22
    for name,r in specs.items():
        if name != audit.REFERENCE_CASE:
            assert len(r["training_assays"]) == (16 if name.startswith("all_batch") else 14)
            assert not set(r["training_assays"]) & set(range(45,85))
            assert set(r["equilibrium_training"]) == {"1","2"}
    eq = specs["all_batch_ratio_reference"]["equilibrium_training"]["1"]
    assert eq["calibration_temperatures_C"] == [40.,50.,60.]
    assert not any(r["predictive_screen_passed"] for r in frozen["batch"]["qualification"].values())


def test_real_julia_feed_transfer_and_previous_result():
    jl = audit.load_runtime()
    plan = read_json(audit.PLAN)
    source = read_json(audit.PUBLISHED.MODEL_INPUTS)
    frozen = read_json(audit.FROZEN)
    rows = audit.conditions(jl,source,read_json(audit.EXTRA),read_json(audit.PROPERTIES))
    assert len(rows) == 40 and {r["actual_mass_kg"] for r in rows} == {.006}
    a,b = [next(r for r in rows if audit.identity(r) == (feed,60,3)) for feed in (1,2)]
    assert a["acid0_mol_L"] == a["alcohol0_mol_L"] == 5.3
    assert b["acid0_mol_L"] == 3.3 and b["alcohol0_mol_L"] == 6.6
    assert a["eta_internal"] == .943 and b["eta_internal"] == .869
    assert isclose(a["catalyst_g_L"]*a["void_volume_L"],6.,rel_tol=1e-15)
    assert isclose(a["film_coefficient_m_s"],b["film_coefficient_m_s"],rel_tol=1e-14)
    assert isclose(a["capacity_mol_L_min"]/b["capacity_mol_L_min"],5.3/3.3,rel_tol=1e-14)
    assert isclose(audit.ester_g_h(1.,a)/audit.ester_g_h(1.,b),5.3/3.3,rel_tol=1e-14)
    slow = audit.transport(jl,a,"half_diffusivity")
    assert isclose(slow["film_coefficient_m_s"]/a["film_coefficient_m_s"],.5**(2/3),rel_tol=1e-14)
    assert slow["actual_mass_kg"] == a["actual_mass_kg"] and slow["eta_internal"] == a["eta_internal"]
    name = "all_batch_ratio_reference"
    spec = audit.case_specs(frozen)[name]
    for c in (a,b):
        kin,K = audit.frozen_model(jl,spec,c,None)
        assert int(kin.feed) == c["feed_ratio"]
        result = audit.solve(jl,kin,c,plan["numerical_acceptance"])
        expected = float(jl.CompositionKinetics.equilibrium_conversion(c["feed_ratio"],K))
        assert 0 < result["predicted_conversion"] < expected
        assert abs(result["integrated_source_fraction"]-result["inlet_backdiffusion_fraction"]-result["predicted_conversion"]) < 1e-9
        if c["feed_ratio"] == 2:
            old = next(r for r in read_rows(audit.PRIOR_PREDICTIONS) if r["case"] == name and audit.PUBLISHED.condition(r) == (60,3))
            assert isclose(result["predicted_conversion"],float(old["predicted_conversion"]),abs_tol=1e-10)
    assert np.isfinite(audit.ester_g_h(result["predicted_conversion"],b))


def test_sign_diagnosis_changes_only_one_printed_coefficient():
    meta=read_json(sign.PROTOCOL)["author_model_digitization"]
    for r in read_rows(sign.AUTHOR[1]):
        y0,y100=meta["conversion_y0_y100_px"]
        assert float(r["author_predicted_conversion"])==round((y0-float(r["line_y_px"]))/(y0-y100),3)
    assert set(sign.author_curves())==set(audit.read_observations())
    source=read_json(audit.PUBLISHED.MODEL_INPUTS)
    literal,alternate=[sign.source_parameters(source,case) for case in sign.LABELS]
    index=audit.PUBLISHED.PARAMETERS.index("Kd_exp")
    assert list(np.flatnonzero(literal!=alternate))==[index]
    assert literal[index] == -28091 and alternate[index] == 28091
    jl=audit.load_runtime()
    p,q=[jl.LorenzoPackedBed.PingPong(a) for a in (literal,alternate)]
    kd=lambda k,T:float(jl.LorenzoPackedBed.rate_terms(k,T,5.3,5.3,False).kd)
    assert kd(p,40)==kd(q,40)
    assert kd(p,80)<kd(p,40)<kd(q,80)
