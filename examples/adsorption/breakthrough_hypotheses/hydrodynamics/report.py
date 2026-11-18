"""Build the compact wall/core column input from verified or exploratory fields."""
import argparse
from pathlib import Path

import numpy as np

from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT

from . import native
from . import packing as p
from .handoff import shell_weights
from .tracer import CUT_DP


def reference(path):
    return dict(path=str(path.relative_to(ROOT)), sha256=sha(path))


def native_ring_flow(folder):
    raw, field = native.field(folder)
    weights = shell_weights(raw["shape"][1], raw["dx_m"], CUT_DP)
    square_weights = weights[field[:, 1].astype(int), field[:, 2].astype(int)]
    return float(np.dot(square_weights, field[:, 3]) / np.sum(field[:, 3]))


def candidate(base, name):
    folder = base / name
    radial = read_json(folder / "radial.json")
    tracer = read_json(folder / "tracer/summary.json")
    checked = read_json(folder / "verified_field.json")
    if (radial["geometry_plan_sha256"] != sha(base / "plan.json") or
            tracer["compact"]["radial_summary_sha256"] != sha(folder / "radial.json") or
            tracer["compact"]["face_summary_sha256"] != sha(folder / "faces/summary.json") or
            checked["field_sha256"] != radial["field_sha256"] or
            radial["field_sha256"] != sha(folder / "field.f64")):
        raise ValueError("Carrier, tracer and packing references must match")
    shell, = [entry for entry in radial["shells"] if entry["cut_dp"] == CUT_DP]
    fine = min(tracer["steps"].values(), key=lambda step: step["dt_s"])
    marker_file = folder / "tracer" / f"markers_dt_{fine['dt_s']:.5f}.f64"
    if fine["marker_sha256"] != sha(marker_file):
        raise ValueError("Frozen marker histories changed")
    marker_data = np.fromfile(marker_file, dtype="<f8").reshape(fine["marker_count"], 11)
    return dict(rate=radial["rate"], name=name, flow=shell["wall_flow_fraction"],
        liquid=shell["wall_void_fraction"], carbon=shell["wall_carbon_fraction"],
        carbon_surface=shell["wall_carbon_surface_fraction"], length_m=radial["period_length_m"],
        native_cell_averaged_wall_flow_fraction=native_ring_flow(folder),
        transition=fine["transitions"]["wall_core_conditional_on_outlet"],
        grid_cells_per_dp=radial["cells_per_dp"], local_shell_porosity=shell["shell_continuous_porosity"],
        radial_shell_flow_fractions=[row["wall_flow_fraction"] for row in radial["shells"]],
        radial_raw_shell_flow_fractions=[row["raw_wall_flow_fraction"] for row in radial["shells"]],
        face_checks=radial["face_checks"], tracer_gates=tracer["gates"],
        tracer_dt_s=fine["dt_s"], marker_count=fine["marker_count"],
        time_step_escape_delta=tracer["dt_escape_fraction_difference"],
        time_step_transition_delta=tracer["dt_transition_probability_max_difference"],
        inert=fine["inert"], tracer_transitions=fine["transitions"], clean_carbon=fine["clean_carbon"],
        integrated_wall_occupation_fraction=float(marker_data[:, 9].sum() / marker_data[:, 3].sum()),
        tracer_cut_reflux=tracer["compact"]["cut_reflux_fraction"],
        tracer_cut_fraction=tracer["compact"]["cut_fraction"],
        native=checked["raw"], carrier_condition=radial["condition"],
        sources={f: reference(folder / f) for f in ("request.json", "build.json", "verified_field.json", "faces/summary.json",
                                                   "radial.json", "tracer/summary.json")})


def compare(coarse, fine):
    flow_change = np.abs(np.asarray(fine["radial_shell_flow_fractions"]) - coarse["radial_shell_flow_fractions"])
    matrix_change = np.abs(np.asarray(fine["transition"]) - coarse["transition"])
    escape_change = abs(fine["clean_carbon"]["outlet_escape_fraction"] - coarse["clean_carbon"]["outlet_escape_fraction"])
    hold_coarse = coarse["inert"]["downstream_mean_residence_s"]
    hold_fine = fine["inert"]["downstream_mean_residence_s"]
    variance_coarse = coarse["inert"]["downstream_residence_variance_s2"]
    variance_fine = fine["inert"]["downstream_residence_variance_s2"]
    flow_max = float(flow_change.max())
    matrix_max = float(matrix_change.max())
    hold = abs(hold_coarse - hold_fine) / hold_fine
    variance = abs(variance_coarse - variance_fine) / variance_fine
    gates = dict(radial_shares=flow_max <= .02, transition=matrix_max <= .02,
                 carbon_escape=escape_change <= .02, average_residence=hold <= .05,
                 residence_variance=variance <= .10)
    return dict(coarse=coarse["name"], fine=fine["name"], radial_flow_share_absolute_differences=flow_change.tolist(),
        max_radial_flow_share_absolute_difference=flow_max, transition_max_absolute_difference=matrix_max,
        clean_carbon_escape_absolute_difference=escape_change, mean_residence_relative_difference=hold,
        variance_relative_difference=variance, gates=gates)


def summarize(base, directories):
    plan, _ = native.geometry(base)
    results = [candidate(base, directory) for directory in directories]
    cases, grids = {}, {}
    for rate in p.RATES:
        data = sorted((r for r in results if r["rate"] == rate), key=lambda x: x["grid_cells_per_dp"])
        if data:
            cases[rate] = data[-1]
        if len(data) >= 2 and data[-2]["grid_cells_per_dp"] < data[-1]["grid_cells_per_dp"]:
            grids[rate] = compare(*data[-2:])
    pass_gates = (set(grids) == set(cases) == set(p.RATES) and
                  all(all(record["gates"].values()) for record in grids.values()) and
                  all(record["face_checks"]["qualifies_local_correction"] and
                      all(record["tracer_gates"].values()) for record in cases.values()))
    return dict(cases=cases, qualification=("numerically checked synthetic scenario" if pass_gates else
        "exploratory synthetic scenario; see correction, reflux and grid gates"), passes_all_gates=pass_gates,
        mesh_sensitivity=grids, audited_runs={r["name"]: r for r in results}, wall_zone_dp=CUT_DP,
        full_column_dry_carbon_kg=plan["full_column_dry_carbon_kg"],
        sources=dict(plan=reference(base / "plan.json"), centres=reference(base / "centres_m.txt")),
        limitations="Material and flow inputs are literature matched, while sphere positions are synthetic. "
                    "A 3.68 mm axial period cannot reproduce the 25 cm specimen or its end distributors. "
                    "Tracer exchange is conditional on downstream exit. Carbon absorption is only the ideal first-hit limit.")


def compact_record(summary, base, control_name):
    keys = ("flow", "native_cell_averaged_wall_flow_fraction", "liquid", "carbon", "carbon_surface", "length_m",
            "transition", "grid_cells_per_dp", "face_checks", "tracer_gates", "integrated_wall_occupation_fraction")
    runs = {}
    for name, data in summary["audited_runs"].items():
        inert, clean, raw = (data[key] for key in ("inert", "clean_carbon", "native"))
        runs[name] = dict({key: data[key] for key in keys},
            flow_mL_min=data["carrier_condition"]["flow_mL_min"], stock_TRT_magic=raw["magic"],
            max_Mach=raw["max_Mach"], max_density_departure=raw["max_density_departure"],
            marker_count=data["marker_count"],
            inert_outlet_fraction=inert["outlet_recovered_fraction"], inert_return_fraction=inert["inlet_return_fraction"],
            inert_censored_fraction=inert["censored_fraction"], mean_transit_s=inert["downstream_mean_residence_s"],
            residence_variance_s2=inert["downstream_residence_variance_s2"],
            clean_escape_fraction=clean["outlet_escape_fraction"],
            clean_escape_by_entry_region=clean["escape_given_inlet_zone"],
            radial_summary=data["sources"]["radial.json"], tracer_summary=data["sources"]["tracer/summary.json"])
    sensitivity = None
    if control_name is not None:
        folder = base / control_name
        radial = read_json(folder / "radial.json")
        checked = read_json(folder / "verified_field.json")
        if radial["geometry_plan_sha256"] != sha(base / "plan.json") or checked["field_sha256"] != sha(folder / "field.f64"):
            raise ValueError("TRT numerical control refers to a different carrier case")
        shell, = [entry for entry in radial["shells"] if entry["cut_dp"] == CUT_DP]
        sensitivity = dict(flow_mL_min=radial["condition"]["flow_mL_min"], stock_TRT_magic=checked["raw"]["magic"],
            tau_even=checked["raw"]["tau_even"], max_Mach=checked["raw"]["max_Mach"],
            max_density_departure=checked["raw"]["max_density_departure"],
            wall_face_flow_fraction=shell["wall_flow_fraction"], wall_cell_flow_fraction=native_ring_flow(folder),
            radial_summary=reference(folder / "radial.json"), verified_field=reference(folder / "verified_field.json"))
    return dict(qualification=summary["qualification"], passes_all_gates=summary["passes_all_gates"],
        wall_zone_dp=summary["wall_zone_dp"], full_column_dry_carbon_kg=summary["full_column_dry_carbon_kg"],
        mesh_sensitivity=summary["mesh_sensitivity"], cases=runs, TRT_magic_sensitivity=sensitivity,
        sources=summary["sources"], limitations=summary["limitations"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("base", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("directories", nargs="+")
    parser.add_argument("--compact", action="store_true", help="Write a compact verification record")
    parser.add_argument("--control", help="Optional separately run TRT numerical control for the compact record")
    args = parser.parse_args()
    if args.control and not args.compact:
        parser.error("--control applies to the compact record")
    base = args.base.resolve()
    summary = summarize(base, args.directories)
    write_json(args.output, compact_record(summary, base, args.control) if args.compact else summary)
    print(summary["qualification"], sorted(summary["cases"]))


if __name__ == "__main__":
    main()
