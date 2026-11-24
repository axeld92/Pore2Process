import math
import numpy as np
from pore2process_demo.io import read_json
from . import audit


def test_table_volume_convention_and_source_scope():
    v=audit.volume_ledger()
    assert len(v["bajic_rows"])==12 and v["bajic_all_within_table_rounding"]
    assert any(r["shape"]=="hexagon" for r in v["bajic_rows"])
    assert all(r["shape"]=="rectangle" for r in v["strnisa_rows"] if r["reactor"] in ("A","B","C"))


def test_diameter_and_volume_units():
    f=audit.fluid_and_ports()
    np.testing.assert_allclose(f["tube_volume_per_cm_ul"],math.pi*.5**2/4*10)
    np.testing.assert_allclose(f["ports"]["Strnisa"]["equivalent_volume_ul"],30*math.pi*.5**2/4)
    np.testing.assert_allclose(f["ports"]["Panjan"]["equivalent_volume_ul"],40*math.pi*.5**2/4)
    np.testing.assert_allclose(f["detector_nominal_hydraulic_time_s"],2.4)
    assert f["diameter_peclet"]>1000 and f["reynolds"]<100


def test_equivalent_length_is_signed_mean_residual():
    f=audit.fluid_and_ports();r={v["case"]:v for v in audit.residual_rows() if v["reading"]=="bg0"}
    assert 89<r["R2"]["equivalent_same_id_length_mm"]<90
    assert r["C"]["equivalent_same_id_length_mm"]<0
    for row in r.values():
        np.testing.assert_allclose(row["equivalent_same_id_length_mm"]*f["area_mm2"],row["mean_residual_s"]*f["flow_ul_s"])


def test_analytic_slope_bound_and_shortening_sign():
    f=audit.fluid_and_ports();D=f["glucose_diffusivity_m2_s"]
    a=audit.variance_bound(dict(mean_residual_s=10.,variance_residual_s2=300.),f["radius_m"],D)
    np.testing.assert_allclose(a["possible_increment_interval_s2"][1],10*f["radius_m"]**2/(24*D))
    assert not a["necessary_variance_screen_pass"] and a["diameter_lower_bound_mm"]>.5
    b=audit.variance_bound(dict(mean_residual_s=-1.,variance_residual_s2=3.),f["radius_m"],D)
    assert b["possible_increment_interval_s2"][1]==0 and not b["necessary_variance_screen_pass"]


def test_even_optimistic_D_and_memory_sensitivity_cannot_supply_R2_variance():
    f=audit.fluid_and_ports();rows=[r for r in audit.residual_rows() if r["case"]=="R2_memory"]
    assert len(rows)==4
    for row in rows:
        bound=audit.variance_bound(row,f["radius_m"],.8*f["glucose_diffusivity_m2_s"])
        assert not bound["necessary_variance_screen_pass"]


def test_retained_actual_julia_finite_extension_and_mesh_screen():
    r=read_json(audit.REF/"summary.json")
    assert r["actual_new_Julia_tube_run"] and not r["new_porous_OpenLB_run"]
    assert all(r["numerical"][k] for k in ("mean_pass","refinement_pass","closed_form_slope_pass"))
    for row in r["residuals"]:
        for v in row["diffusivity_cases"].values():
            if row["mean_residual_s"]>0:
                assert 0<=v["independent_added_segment_variance_s2"]<=v["contiguous_extension_from_entire_blank_equivalent_s2"]+1e-7
            else:
                assert v["contiguous_extension_from_entire_blank_equivalent_s2"]<=0
    assert not r["any_same_id_tube_necessary_variance_pass"]
    middle=next(v for v in r["residuals"] if v["case"]=="R2" and v["reading"]=="bg0")
    assert not middle["diffusivity_cases"]["0.8"]["full_depth_mean_equivalent_bore_necessary_pass"]
    assert middle["diffusivity_cases"]["1.0"]["cylindrical_length_at_necessary_diameter_mm"]<20.
