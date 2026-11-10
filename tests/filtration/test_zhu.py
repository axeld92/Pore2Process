"""Digitized real-data fit, independent series oracle and saved-result replay."""
from pore2process_demo.paths import project_path
from copy import deepcopy
import math

import pytest

from examples.filtration.calibration.digitize_zhu import digitize, pixels_to_flux
from examples.filtration.calibration.run import ROOT, write_json
from examples.filtration.calibration.zhu import analyze, check, export, load_snapshot, validate
from pore2process_demo.bridge import runtime


def test_actual_julia_fit_oracle_and_heldout_isolation():
    snapshot=load_snapshot()
    r=analyze(snapshot,with_sensitivity=False)
    train=[c for c in r["cases"] if c["role"]=='calibration']
    expected=math.fsum(1e5/c["measured_flux_m_s"]-1/c["support_permeance_m_Pa_s"] for c in train)/len(train)
    assert r["fit"]["fixed_impedance_Pa_s_m"]==pytest.approx(expected,rel=1e-12)
    assert len(train)==11
    assert [c["point_id"] for c in r["cases"] if c["role"]=="holdout"]==list(range(2,22,2))
    assert r["metrics"]["holdout"]["max_absolute_relative_error"]<.10
    assert r["screen"]["nominal_passed"] is True
    for c in r["cases"]:
        expected_flux=1e5/(1/c["support_permeance_m_Pa_s"]+expected)
        assert c["predicted_flux_m_s"]==pytest.approx(expected_flux,rel=1e-12)
    # A synthetic held-out pixel/flux change must leave the fitted model alone.
    changed=deepcopy(snapshot)
    point=changed["digitized"]["points"][1]
    point["pixel_y"]-=2
    point["support_flux_m_h"],point["total_flux_LMH"]=pixels_to_flux(changed["digitized"]["config"],point["pixel_x"],point["pixel_y"])
    altered=analyze(changed,with_sensitivity=False)
    assert altered["fit"]==r["fit"]
    assert altered["metrics"]["holdout"]!=r["metrics"]["holdout"]
    assert snapshot["source"]["measurement_conditions"]["water_viscosity_Pa_s"] is None
    assert snapshot["source"]["measurement_conditions"]["water_temperature_K"] is None
    # Independent impedance-vs-resistance unit equivalence, using a synthetic μ.
    j=runtime()
    mu,L,Ps,Z,J,A,eta=2e-3,.003,1e-8,2e8,1e-5,2.,.5
    i=j.impedance_load(Ps,Z,J,A,eta)
    f=j.filtration_load(mu*L*Ps,L,Z/mu,mu,J,A,eta)
    assert float(i.pressure_drop_Pa)==pytest.approx(float(f.total_pressure_drop_Pa))
    assert float(i.shaft_power_W)==pytest.approx(float(f.shaft_power_W))


def test_digitization_reproduction_and_invalid_snapshot():
    snapshot=load_snapshot()
    recorded=snapshot["digitized"]
    if (project_path(ROOT/recorded["config"]["image"])).exists():
        assert digitize(recorded["config"])==recorded
        bad=deepcopy(recorded["config"])
        bad["image_sha256"]="0"*64
        with pytest.raises(ValueError,match="image changed"):
            digitize(bad)
    changed=deepcopy(snapshot)
    changed["digitized"]["points"][0]["total_flux_LMH"]*=1000
    with pytest.raises(ValueError,match="pixel/axis"):
        validate(changed)
    changed=deepcopy(snapshot)
    changed["digitized"]["points"][0]["point_id"]=2
    with pytest.raises(ValueError,match="Point IDs"):
        validate(changed)


def test_pdf_refinement_preserves_split_and_uses_actual_julia(tmp_path):
    snapshot=load_snapshot(source_reading="pdf")
    original=load_snapshot(source_reading="jpeg")
    assert snapshot["protocol"]["protocol_id"]!=original["protocol"]["protocol_id"]
    key="max_absolute_relative_flux_error"
    assert snapshot["protocol"]["primary_screen"][key]==original["protocol"]["primary_screen"][key]==.10
    assert snapshot["protocol"]["digitization"]["marker_center_bound_px"]==2
    assert snapshot["source"]["source"]["original_pdf_received"] is True
    assert [r["point_id"] for r in snapshot["digitized"]["points"]]==list(range(1,22))
    assert snapshot["digitized"]["config"]["image"]!=original["digitized"]["config"]["image"]
    if (project_path(ROOT/snapshot["digitized"]["config"]["image"])).exists():
        assert digitize(snapshot["digitized"]["config"])==snapshot["digitized"]
    result=analyze(snapshot)
    train=[c for c in result["cases"] if c["role"]=='calibration']
    expected=math.fsum(1e5/c["measured_flux_m_s"]-1/c["support_permeance_m_Pa_s"] for c in train)/len(train)
    assert result["fit"]["fixed_impedance_Pa_s_m"]==pytest.approx(expected,rel=1e-12)
    assert [c["point_id"] for c in result["cases"] if c["role"]=="holdout"]==list(range(2,22,2))
    assert result["screen"]["nominal_passed"] is True
    assert result["screen"]["all_reading_corners_passed"] is True
    assert all(r["experimental_uncertainty"] is None for r in snapshot["digitized"]["points"])
    assert snapshot["source"]["measurement_conditions"]["water_temperature_K"] is None
    assert snapshot["source"]["measurement_conditions"]["water_viscosity_Pa_s"] is None
    output=tmp_path/"pdf-run"
    export(snapshot,result,output)
    assert check(output)["integrity_and_numerical_replay"]=="passed"


def test_exports_sensitivity_energy_and_tamper(tmp_path):
    snapshot=load_snapshot()
    result=analyze(snapshot)
    output=tmp_path/"run"
    saved=export(snapshot,result,output)
    assert check(output)["integrity_and_numerical_replay"]=="passed"
    assert saved["runtime_versions"]["PythonCall"]==saved["runtime_versions"]["juliacall"]
    assert result["reported_endpoint_energy_translation"]["M1"]["conditional_fixed_flux_energy_saving_fraction"]==pytest.approx(3/7)
    low,high=result["digitization_sensitivity"]["fixed_impedance_range_Pa_s_m"]
    assert low<result["fit"]["fixed_impedance_Pa_s_m"]<high
    # 90% effective non-support share at nine times its limiting flux.
    j=runtime()
    flux=result["threshold_support_flux_at_1bar_m_h"]["non_support_90_percent"]
    value=j.series_flux(flux/3600/1e5,result["fit"]["fixed_impedance_Pa_s_m"],1e5)
    assert float(value.non_support_fraction)==pytest.approx(.9)
    loads=result["conditional_process"]
    assert loads[0]["energy_saving_vs_low_support_fraction"]==0.
    for a,b in zip(loads,loads[1:]):
        assert a["pressure_drop_Pa"]>b["pressure_drop_Pa"]
        assert a["flow_m3_s"]==b["flow_m3_s"]
    with pytest.raises(FileExistsError):
        export(snapshot,result,output)
    modified=deepcopy(saved)
    modified["fit"]["fixed_impedance_Pa_s_m"]*=2
    write_json(output/"summary.json",modified)
    with pytest.raises(ValueError,match="Different numerical value"):
        check(output)
    write_json(output/"summary.json",saved)
    with (output/"predictions.csv").open("a") as stream:
        stream.write("altered\n")
    with pytest.raises(ValueError,match="hash mismatch"):
        check(output)
