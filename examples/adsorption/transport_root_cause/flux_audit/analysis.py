"""Audit the exact discrete balance of archived native populations.

This module never reconstructs a carrier flux from cell-centred velocity.
The C++ diagnostic reads populations immediately before and after the stock
collision and immediately after its stock streaming step. Full-way BounceBack
wall populations are included explicitly as numerical inventory.
"""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np


def summarize(path: str | Path, *, dx_m: float, dt_s: float,
              density_scale: float = 1.0) -> dict:
    """Return balance diagnostics; density_scale converts native density to SI.

    For the passive scalar, density_scale=1 kg/m3 per native concentration.
    For the present water solve, density_scale=1000 kg/m3 per lattice density.
    Post-population and momentum values remain in native lattice units.
    """
    if not all(np.isfinite(x) and x > 0 for x in (dx_m, dt_s, density_scale)):
        raise ValueError("dx_m, dt_s and density_scale must be positive and finite")
    with Path(path).open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError("population audit is empty")
    arrays = {name: np.array([float(row[name]) for row in rows]) for name in rows[0]}
    if not all(np.all(np.isfinite(values)) for values in arrays.values()):
        raise ValueError("population audit contains non-finite values")
    keys = list(zip(arrays["i"], arrays["j"], arrays["k"], strict=True))
    if len(set(keys)) != len(keys):
        raise ValueError("population audit contains duplicate node IDs")
    wall = arrays["material"] == 2
    fluid = ~wall
    adjacent = fluid & (arrays["wall_adjacent_links"] > 0)
    bulk = fluid & ~adjacent
    volume = dx_m ** 3
    before = arrays["concentration_before"]
    post = arrays["concentration_after_collision"]
    predicted = arrays["predicted_after_stream"]
    after = arrays["concentration_after_stream"]
    outgoing = arrays["net_to_fluid_per_s"] + arrays["net_to_wall_per_s"]
    collision = post - before
    stream_residual = after - predicted
    link_residual = predicted - post + outgoing * dt_s / volume
    density_tolerance = 1e-11 * max(1., float(np.max(np.abs(before))))

    def maxabs(a: np.ndarray) -> float:
        return float(np.max(np.abs(a), initial=0.))

    def group(mask: np.ndarray) -> dict:
        count = int(np.count_nonzero(mask))
        return {
            "nodes": count,
            "min_density_after_native": float(np.min(after[mask])) if count else None,
            "max_density_after_native": float(np.max(after[mask])) if count else None,
            "max_abs_density_change_native": maxabs((after-before)[mask]),
            "max_abs_collision_change_native": maxabs(collision[mask]),
            "max_abs_stream_residual_native": maxabs(stream_residual[mask]),
            "max_abs_link_balance_residual_native": maxabs(link_residual[mask]),
            "negative_post_population_nodes": int(np.count_nonzero(arrays["negative_post_populations"][mask] > 0)),
            "negative_post_populations": int(np.sum(arrays["negative_post_populations"][mask])),
            "min_unshifted_post_population_native": float(np.min(arrays["min_unshifted_post_population"][mask])) if count else None,
            "inventory_change_kg": float(np.sum((after-before)[mask]) * volume * density_scale),
        }

    fluid_to_wall = float(np.sum(arrays["net_to_wall_per_s"][fluid]))
    wall_to_fluid = float(np.sum(arrays["net_to_fluid_per_s"][wall]))
    cancellation_tolerance = density_tolerance * len(rows) * volume / dt_s
    return {
        "schema": "pore2process.exact_population_audit.v1",
        "step": int(arrays["step"][0]),
        "dx_m": dx_m, "dt_s": dt_s, "density_scale_kg_m3": density_scale,
        "all": group(np.ones(len(rows), dtype=bool)), "fluid": group(fluid),
        "wall_adjacent_fluid": group(adjacent), "bulk_fluid": group(bulk), "wall": group(wall),
        "fluid_to_wall_net_kg_s": fluid_to_wall * density_scale,
        "wall_to_fluid_net_kg_s": wall_to_fluid * density_scale,
        "interface_antisymmetry_residual_kg_s": (fluid_to_wall+wall_to_fluid)*density_scale,
        "global_net_link_outflow_kg_s": float(np.sum(outgoing))*density_scale,
        "max_abs_collision_mass_change_native": maxabs(collision),
        "max_abs_stream_residual_native": maxabs(stream_residual),
        "max_abs_link_balance_residual_native": maxabs(link_residual),
        "density_tolerance_native": density_tolerance,
        "checks": {
            "collision_preserves_local_density": maxabs(collision) <= density_tolerance,
            "streaming_matches_actual_population_links": maxabs(stream_residual) <= density_tolerance,
            "local_link_balance_closes": maxabs(link_residual) <= density_tolerance,
            "wall_interface_flux_is_antisymmetric": abs(fluid_to_wall+wall_to_fluid) <= cancellation_tolerance,
            "global_periodic_link_flux_cancels": abs(float(np.sum(outgoing))) <= cancellation_tolerance,
        },
        "interpretation": (
            "These are exact population balances, not centred div(u). Local density changes and "
            "instantaneous fluid-wall exchange are allowed by full-way BounceBack numerical storage. "
            "Passing conservation does not establish concentration boundedness or transport accuracy. "
            "Negative individual populations are a diagnostic, not by themselves an acceptance criterion."
        ),
    }


def summarize_uniform_stencil(path: str | Path, *, first_step_audit: str | Path | None = None,
                              initial_concentration: float = .5) -> dict:
    """Assess a stock-equilibrium constant-preservation stencil.

    The optional comparison must be an actual initially constant scalar run's
    first step; a sine-profile run must not be passed as that comparison.
    The stencil reads stock equilibrium populations, not an invented div(u).
    """
    with Path(path).open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError("uniform stencil audit is empty")
    fluid_rows = [row for row in rows if int(row["material"]) != 2]
    if not fluid_rows:
        raise ValueError("uniform stencil audit contains no fluid")
    change = np.array([float(row["uniform_change_native"]) for row in fluid_rows])
    divergence = np.array([float(row["scalar_equilibrium_divergence_per_s"]) for row in fluid_rows])
    if not np.all(np.isfinite(change)) or not np.all(np.isfinite(divergence)):
        raise ValueError("uniform stencil audit contains non-finite values")
    result = {
        "fluid_nodes": len(fluid_rows),
        "max_abs_one_step_relative_constant_change": float(np.max(np.abs(change))),
        "rms_one_step_relative_constant_change": float(np.sqrt(np.mean(change**2))),
        "max_abs_scalar_equilibrium_divergence_per_s": float(np.max(np.abs(divergence))),
        "rms_scalar_equilibrium_divergence_per_s": float(np.sqrt(np.mean(divergence**2))),
        "first_step_comparison": None,
        "interpretation": (
            "A read-only exact first-order-equilibrium streaming stencil. This predicts the first "
            "step of an initially uniform scalar, not the evolved long-time scalar flux and not "
            "the carrier's population flux. A nonzero value demonstrates that exported cell "
            "velocities are not exactly divergence-free under this scalar equilibrium stencil."
        ),
    }
    if first_step_audit is not None:
        with Path(first_step_audit).open(newline="") as stream:
            actual = {(int(row["i"]),int(row["j"]),int(row["k"])): row for row in csv.DictReader(stream)}
        errors = []
        for row in fluid_rows:
            key = tuple(int(row[name]) for name in ("i","j","k"))
            if key not in actual:
                raise ValueError("uniform stencil and first-step audit nodes differ")
            audit = actual[key]
            if int(audit["step"]) != 1 or abs(float(audit["concentration_before"])-initial_concentration) > 1e-12:
                raise ValueError("comparison requires the first step of an initially uniform scalar")
            predicted = initial_concentration*float(row["predicted_uniform_after_stream_native"])
            errors.append(float(audit["concentration_after_stream"])-predicted)
        max_error = float(np.max(np.abs(errors)))
        result["first_step_comparison"] = {
            "initial_concentration_kg_m3": initial_concentration,
            "max_abs_prediction_error_kg_m3": max_error,
            "matches_to_roundoff": max_error <= 1e-12*max(1., abs(initial_concentration)),
        }
    return result
