"""Real Julia fit, independent oracle, held-out isolation and artifact replay."""
from copy import deepcopy
import math

import pytest

from calibration.run import analyze, export, load_snapshot, read_json, validate_snapshot, write_json
from calibration.check import check


def test_julia_fit_against_independent_oracle_and_no_holdout_leakage():
    snapshot=load_snapshot(require_sources=False)
    result=analyze(snapshot)
    rows=snapshot["source"]["table2_rows"]
    train=[r for r in rows if r["table_row"] in (2,3,4)]
    # Independent product/geometric-mean expression, ordinary direct powers.
    expected=math.prod(r["permeability_m2"]/(400e-9**2*r["porosity_fraction"]**3/(180*(1-r["porosity_fraction"])**2)) for r in train)**(1/3)
    assert result["fit"]["multiplier"]==pytest.approx(expected,rel=1e-12)
    assert result["fit"]["calibration_count"]==3
    holdout=[c for c in result["cases"] if c["role"]=="primary_holdout"]
    assert [c["table_row"] for c in holdout]==[5,6,7]
    for c in holdout:
        e=c["porosity_fraction"]
        oracle=expected*400e-9**2*e**3/(180*(1-e)**2)
        assert c["predicted_permeability_m2"]==pytest.approx(oracle,rel=1e-12)
    assert result["metrics"]["primary_holdout"]["calibrated"]["mean_absolute_relative_error"]==pytest.approx(sum(abs(c["relative_error"]) for c in holdout)/3)
    assert result["primary_screen"]=={"threshold":.25,"passed":False,"failed_rows":[5]}
    assert all(c["permeability_uncertainty_m2"] is None for c in result["cases"])
    # Synthetic mutation of a holdout value AND its mock audit reference.
    # No real source file is changed; this specifically tests fit isolation.
    changed=deepcopy(snapshot)
    changed["source"]["table2_rows"][4]["permeability_m2"]*=10
    changed["source"]["table2_rows"][4]["reported_k_mantissa_times_1e_minus16_m2"]*=10
    changed["audit"]["visually_transcribed_table2_rows"][4][3]*=10
    altered=analyze(changed)
    assert altered["fit"]==result["fit"]
    assert altered["leave_one_calibration_out"]==result["leave_one_calibration_out"]
    assert altered["metrics"]["primary_holdout"]!=result["metrics"]["primary_holdout"]
    assert altered["cases"][4]["predicted_permeability_m2"]==result["cases"][4]["predicted_permeability_m2"]


def test_source_audit_and_split_rejection():
    snapshot=load_snapshot()
    changed=deepcopy(snapshot)
    changed["source"]["table2_rows"][1]["permeability_m2"]*=1e16
    with pytest.raises(ValueError,match="multiplier"):
        validate_snapshot(changed)
    changed=deepcopy(snapshot)
    changed["protocol"]["primary_holdout_rows"].append(2)
    with pytest.raises(ValueError,match="overlap"):
        validate_snapshot(changed)
    changed=deepcopy(snapshot)
    changed["protocol"]["nominal_primary_particle_diameter_m"]["AKP-30"]*=2
    with pytest.raises(ValueError,match="particle sizes"):
        validate_snapshot(changed)


def test_export_replay_conditional_energy_and_tamper_rejection(tmp_path):
    snapshot=load_snapshot()
    result=analyze(snapshot)
    summary=export(snapshot,result,tmp_path/"run")
    assert check(tmp_path/"run")["integrity_and_numerical_replay"]=="passed"
    with pytest.raises(ValueError,match="scientific screen failed"):
        check(tmp_path/"run",require_screen_pass=True)
    with pytest.raises(FileExistsError):
        export(snapshot,result,tmp_path/"run")
    assert read_json(tmp_path/"run/inputs.json")==snapshot
    assert summary["runtime_versions"]["julia"]
    assert summary["runtime_versions"]["PythonCall"]==summary["runtime_versions"]["juliacall"]
    measured_saving=[]
    for record in result["conditional_process"]:
        measured=record["measured_support_k"]
        alpha=record["selective_to_measured_baseline_support_ratio"]
        # Independent exact series-resistance expression; no fitted b enters.
        expected=(1-1.35/1.90)/(1+alpha)
        assert measured["fractional_shaft_energy_saving"]==pytest.approx(expected)
        measured_saving.append(expected)
        for key in ("measured_support_k","calibrated_support_k"):
            first,second=record[key]["baseline"],record[key]["alternative"]
            assert first["selective_resistance_m_inv"]==second["selective_resistance_m_inv"]==record["fixed_selective_resistance_m_inv"]
            assert first["flow_m3_s"]==pytest.approx(100e-3/3600)
            assert first["specific_shaft_energy_kWh_m3"]==pytest.approx(first["total_pressure_drop_Pa"]/(.75*3.6e6))
    assert measured_saving==sorted(measured_saving,reverse=True)
    assert result["conditional_process"][-1]["measured_support_k"]["any_total_pressure_exceeds_reported_cell_pressurization"] is True
    changed=deepcopy(summary)
    changed["fit"]["multiplier"]*=2
    write_json(tmp_path/"run/summary.json",changed)
    with pytest.raises(ValueError,match="Different numerical value"):
        check(tmp_path/"run")
    write_json(tmp_path/"run/summary.json",summary)
    with (tmp_path/"run/predictions.csv").open("a") as stream:
        stream.write("tampered\n")
    with pytest.raises(ValueError,match="hash mismatch"):
        check(tmp_path/"run")
    with pytest.raises(ValueError,match="Nonstandard JSON"):
        read_json_text=tmp_path/"bad.json"
        read_json_text.write_text('{"k":NaN}')
        read_json(read_json_text)
