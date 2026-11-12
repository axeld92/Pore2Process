"""Independent accounting fixtures for the read-only native population audit."""
import csv

import pytest

from examples.adsorption.transport_root_cause.flux_audit.analysis import summarize, summarize_uniform_stencil


def audit_file(tmp_path, *, residual=0., negative=0):
    # One fluid cell sends 0.01 native inventory to a full-way wall cell.
    # This is not a production transport model: it tests interface accounting.
    rows = []
    for i, material in enumerate((1, 2)):
        delta = -.01 if material == 1 else .01
        rows.append(dict(step=1, i=i, j=0, k=0, material=material,
                         wall_adjacent_links=1, concentration_before=.5,
                         concentration_after_collision=.5,
                         predicted_after_stream=.5+delta,
                         concentration_after_stream=.5+delta+(residual if i == 0 else 0),
                         collision_delta=0., expected_stream_delta=delta,
                         local_balance_residual=residual if i == 0 else 0,
                         net_to_fluid_per_s=-.01 if material == 2 else 0.,
                         net_to_wall_per_s=.01 if material == 1 else 0.,
                         min_unshifted_post_population=-.001 if negative else .001,
                         negative_post_populations=negative,
                         post_momentum_x=0., post_momentum_y=0., post_momentum_z=0.))
    path = tmp_path / "audit.csv"
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0])
        writer.writeheader(); writer.writerows(rows)
    return path


def test_exact_wall_storage_is_not_a_mass_leak(tmp_path):
    result = summarize(audit_file(tmp_path), dx_m=1., dt_s=1., density_scale=1000.)
    assert all(result["checks"].values())
    assert result["fluid"]["inventory_change_kg"] == pytest.approx(-10.)
    assert result["wall"]["inventory_change_kg"] == pytest.approx(10.)
    assert result["all"]["inventory_change_kg"] == pytest.approx(0.)
    assert result["fluid_to_wall_net_kg_s"] == pytest.approx(10.)
    assert result["wall_to_fluid_net_kg_s"] == pytest.approx(-10.)


def test_local_streaming_error_is_detected(tmp_path):
    result = summarize(audit_file(tmp_path, residual=1e-5), dx_m=1., dt_s=1.)
    assert not result["checks"]["streaming_matches_actual_population_links"]
    assert result["max_abs_stream_residual_native"] == pytest.approx(1e-5)


def test_negative_population_is_reported_separately(tmp_path):
    result = summarize(audit_file(tmp_path, negative=2), dx_m=1., dt_s=1.)
    assert all(result["checks"].values())
    assert result["fluid"]["negative_post_populations"] == 2
    assert result["fluid"]["min_density_after_native"] > 0


@pytest.mark.parametrize("scales", [(0., 1., 1.), (1., float("nan"), 1.), (1., 1., -1.)])
def test_invalid_units_rejected(tmp_path, scales):
    with pytest.raises(ValueError, match="positive and finite"):
        summarize(audit_file(tmp_path), dx_m=scales[0], dt_s=scales[1], density_scale=scales[2])


def stencil_file(tmp_path):
    path=tmp_path / "stencil.csv"
    path.write_text("i,j,k,material,wall_adjacent_links,uniform_before_native,predicted_uniform_after_stream_native,uniform_change_native,scalar_equilibrium_divergence_per_s\n"
                    "0,0,0,1,1,1,0.98,-0.02,0.02\n"
                    "1,0,0,2,1,1,1.02,0.02,-0.02\n")
    return path


def test_stock_stencil_comparison_uses_actual_first_step(tmp_path):
    result=summarize_uniform_stencil(stencil_file(tmp_path), first_step_audit=audit_file(tmp_path))
    assert result["max_abs_one_step_relative_constant_change"] == pytest.approx(.02)
    assert result["first_step_comparison"]["matches_to_roundoff"]


def test_stock_stencil_rejects_nonconstant_comparison(tmp_path):
    path=audit_file(tmp_path)
    contents=path.read_text().replace("0.5,0.5", "0.6,0.5")
    path.write_text(contents)
    with pytest.raises(ValueError, match="initially uniform"):
        summarize_uniform_stencil(stencil_file(tmp_path), first_step_audit=path)
