"""Check the stock TRT water tracer against a distinct experiment, without conversion data."""
import argparse
import csv
from math import gcd, pi
import os
from pathlib import Path
import subprocess

import numpy as np

from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT

from .flow import GRIDS, HERE, OUTPUT, WATER_VISCOSITY_PA_S

SOURCE = ROOT / "literature/sources/lorenzo2025"
BOLTZMANN_J_K = 1.380649e-23
DEXTRAN_DIAMETER_M = 80e-9  # From Section 2.2, as a hydrodynamic-size scenario.
DIFFUSIVITY_M2_S = BOLTZMANN_J_K * 298.15 / (3 * pi * WATER_VISCOSITY_PA_S * DEXTRAN_DIAMETER_M)
SEED = 426620
STOP_S = 5400
MEAN_LIMIT = 0.05
VARIANCE_LIMIT = 0.10


def best_common_cut(fields):
    """Use a common physical cut with minimum worst-case upstream net flux."""
    widths = [field["shape"][0] for field in fields]
    common_width = gcd(*widths)
    phases = np.arange(common_width)
    cases = [
        np.asarray(field["axial_plane_reflux_fractions"])[phases * (width // common_width)]
        for field, width in zip(fields, widths, strict=True)
    ]
    phase = int(np.argmin(np.max(cases, axis=0)))
    return phase, common_width


def tag(grid, dt, diffusivity_factor):
    return f"N{grid}_dt{dt:.2f}_D{diffusivity_factor:g}"


def trace(output, plan, field, cut, dt, diffusivity_factor, julia, threads, count, reuse):
    case = output / f"N{field['grid']}"
    name = tag(field["grid"], dt, diffusivity_factor)
    filename = output / f"{name}.f64"
    specification = output / f"{name}.json"
    dx = field["native"]["dx_m"]
    diffusivity = DIFFUSIVITY_M2_S * diffusivity_factor
    command = [julia, "--startup-file=no", "--project=julia", f"--threads={threads}",
               str(HERE / "tracer.jl"), str(case / "compact"), str(filename),
               str(field["shape"][0]), str(cut), *map(str, (
                   dx, plan["column_length_m"], diffusivity, dt, STOP_S, count, SEED,
               ))]
    hashes = {path.name: sha(path) for path in sorted((case / "compact").iterdir())}
    if hashes != field["compact_files_sha256"]:
        raise RuntimeError(f"Saved face data changed in N{field['grid']}")
    request = dict(command=command, compact_files_sha256=hashes,
                   tracer_jl_sha256=sha(HERE / "tracer.jl"))
    if reuse and filename.is_file() and specification.is_file():
        stored = read_json(specification)
        if stored != dict(**request, sample_sha256=sha(filename)):
            raise RuntimeError(f"Saved {name} trajectory settings or output changed")
    else:
        print(f"Tracking {name} through the synthetic 14 cm bed", flush=True)
        env = {**os.environ, "OPENBLAS_NUM_THREADS": "1"}
        with (output / f"{name}.log").open("w") as log:
            subprocess.run(command, cwd=ROOT, env=env, check=True, stdout=log,
                           stderr=subprocess.STDOUT, timeout=1800)
        write_json(specification, dict(**request, sample_sha256=sha(filename)))
    samples = np.fromfile(filename, dtype="<f8").reshape(-1, 2)
    if samples.shape != (count, 2) or np.any(samples[np.isfinite(samples[:, 0]), 0] <= 0):
        raise RuntimeError(f"Invalid exit times for {name}")
    times_min = samples[:, 0] / 60
    record = dict(
        grid=field["grid"],
        release_cut=cut,
        release_cut_backflow_relative=field["axial_plane_reflux_fractions"][cut],
        returned_upstream_fraction=float(samples[:, 1].mean()),
        timestep_s=dt,
        assumed_diffusivity_m2_s=diffusivity,
        marker_count=count,
        maximum_track_min=STOP_S / 60,
        seed=SEED,
        age_sample_sha256=sha(filename),
    )
    return name, times_min, record


def moments(times, pulse_min, q=0.95):
    """Add the exact moments of an independent rectangular inlet pulse."""
    survived = times[np.isfinite(times)]
    if len(survived) < 2:
        raise RuntimeError("Insufficient recovered tracer")
    mean = float(np.mean(survived))
    variance = float(np.var(survived))
    rng = np.random.default_rng(SEED)
    draws = rng.choice(survived, size=(800, len(survived)))
    confidence = [(1 - q) / 2, 1 - (1 - q) / 2]
    return dict(
        recovery=float(len(survived) / len(times)),
        pulse_free_mean_min=mean,
        pulse_free_variance_min2=variance,
        pulse_free_relative_variance=variance / mean**2,
        rectangular_pulse_mean_min=mean + pulse_min / 2,
        rectangular_pulse_variance_min2=variance + pulse_min**2 / 12,
        pulse_free_simulation_sampling_mean_ci95_min=(
            np.quantile(draws.mean(axis=1), confidence).tolist()
        ),
        pulse_free_simulation_sampling_variance_ci95_min2=(
            np.quantile(draws.var(axis=1), confidence).tolist()
        ),
        right_censored_count=len(times) - len(survived),
        pulse_free_lower_bound_mean_min=float(np.nan_to_num(times, nan=STOP_S / 60).mean()),
    )


def pulse_response(times, ages, pulse):
    delays = times[:, None] - ages
    e = ((delays >= 0) & (delays < pulse)).mean(axis=1) / pulse
    f = np.nan_to_num(np.clip(delays / pulse, 0., 1.)).mean(axis=1)
    return e, f


def tracer_curve(output, ages_min, pulse_min, water_flow_ml_min, feed_concentration):
    """Compare the model to the digitized Figure A1a, keeping the measured time axis."""
    table = SOURCE / "water_tracer.csv"
    with table.open(newline="") as stream:
        rows = [row for row in csv.DictReader(stream) if row["tracer"] == "blue_dextran"]
    measured = sorted((float(row["time_min"]), float(row["concentration_mg_ml"])) for row in rows)
    times, concentrations = np.asarray(measured).T
    if np.any(np.diff(times) <= 0) or times[-1] > STOP_S / 60:
        raise RuntimeError("Figure A1a digitized time axis is invalid")
    integrated = np.cumsum(np.r_[0., np.diff(times) * (concentrations[:-1] + concentrations[1:]) / 2])
    area = integrated[-1]
    measured_f = integrated / area
    predicted_e, predicted_f = pulse_response(times, ages_min, pulse_min)
    columns = ("time_min", "figure_A1a_concentration_mg_ml", "model_concentration_mg_ml",
               "figure_A1a_cdf_normalized_to_digitized_area", "model_cdf")
    modeled = predicted_e * pulse_min * feed_concentration
    with (output / "tracer_curve.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(columns)
        for row in zip(times, concentrations, modeled, measured_f, predicted_f, strict=True):
            writer.writerow(f"{value:.9g}" for value in row)
    return dict(
        measured_trace_sha256=sha(table),
        measured_marker_count=len(times),
        measured_concentration_auc_mg_min_ml=float(area),
        measured_digitized_mass_mg=float(water_flow_ml_min * area),
        printed_blue_dextran_pulse_dose_mg=float(water_flow_ml_min * pulse_min * feed_concentration),
        largest_absolute_cdf_difference=float(np.max(np.abs(predicted_f - measured_f))),
        digitized_mean_min=float(np.trapezoid(times * concentrations, times) / area),
        digitized_variance_min2=float(
            np.trapezoid(times**2 * concentrations, times) / area
            - (np.trapezoid(times * concentrations, times) / area)**2
        ),
        predicted_concentration_rmse_mg_ml=float(np.sqrt(np.mean((modeled - concentrations)**2))),
    )


def check_moments(candidate, reference, pulse):
    prefix = "rectangular_pulse" if pulse else "pulse_free"
    mean_error = candidate[f"{prefix}_mean_min"] / reference["mean_min"] - 1
    variance_error = candidate[f"{prefix}_variance_min2"] / reference["variance_min2"] - 1
    return dict(mean_relative_error=mean_error,
                variance_relative_error=variance_error,
                mean_passed=abs(mean_error) <= MEAN_LIMIT,
                variance_passed=abs(variance_error) <= VARIANCE_LIMIT,
                recovery_passed=candidate["recovery"] >= .995,
                passed=(abs(mean_error) <= MEAN_LIMIT and abs(variance_error) <= VARIANCE_LIMIT
                        and candidate["recovery"] >= .995))


def refinement(candidate, reference):
    record = check_moments(candidate, dict(mean_min=reference["pulse_free_mean_min"],
                                           variance_min2=reference["pulse_free_variance_min2"]), False)
    record["comparator_recovery"] = reference["recovery"]
    record["passed"] &= reference["recovery"] >= .995
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--julia", default="julia")
    parser.add_argument("--threads", type=int, default=6)
    parser.add_argument("--markers", type=int, default=8192)
    parser.add_argument("--reuse", action="store_true")
    parser.add_argument("--trace-only", action="store_true",
                        help="Run the simulation without opening any tracer measurements")
    args = parser.parse_args()
    if args.markers < 256 or args.threads < 1:
        parser.error("Use at least 256 markers and one Julia thread")
    output = args.output.resolve()
    plan = read_json(output / "plan.json")
    fields = [read_json(output / f"N{n}/field_summary.json") for n in GRIDS]
    # The release plane was selected on N32 and N48 before any tracer result.
    # Reuse it on the refinement grid to avoid shifting the operating plane.
    phase, denominator = best_common_cut(fields[:2])
    for field in fields:
        if field["shape"][0] % denominator:
            raise RuntimeError("The chosen physical release plane is not on the refinement mesh")
    pulse_min = plan["injection_duration_min"]
    settings = (
        (0, .20, 1.), (1, .20, 1.), (1, .10, 1.), (1, .20, 2.),
        (2, .20, 1.), (2, .10, 1.),
    )
    samples, records = {}, {}
    for grid_index, dt, factor in settings:
        field = fields[grid_index]
        cut = phase * (field["shape"][0] // denominator)
        name, ages, metadata = trace(
            output, plan, field, cut, dt, factor, args.julia, args.threads, args.markers, args.reuse,
        )
        samples[name] = ages
        records[name] = dict(**metadata, **moments(ages, pulse_min))
        print(name, records[name], flush=True)
    if args.trace_only:
        return
    primary_name = tag(GRIDS[-1], .10, 1.)
    primary = records[primary_name]
    comparisons = dict(
        coarse_mesh_at_matched_step=refinement(records[tag(32, .20, 1.)], records[tag(48, .20, 1.)]),
        mesh_at_matched_step=refinement(records[tag(48, .10, 1.)], primary),
        fine_mesh_at_coarser_step=refinement(records[tag(48, .20, 1.)], records[tag(64, .20, 1.)]),
        step_at_matched_mesh=refinement(records[tag(64, .20, 1.)], primary),
        double_diffusivity_scenario=refinement(records[tag(GRIDS[1], .20, 2.)],
                                              records[tag(GRIDS[1], .20, 1.)]),
    )
    evidence = read_json(SOURCE / "extraction.json")["water_tracer"]
    np.testing.assert_allclose(evidence["water_flow_ml_min"], plan["column_flow_m3_s"] * 60e6)
    np.testing.assert_allclose(evidence["pulse_duration_min_if_injected_at_water_flow"], pulse_min)
    observed = evidence["reported_rtd_moments"]["blue_dextran"]
    comparison = check_moments(primary, observed, True)
    comparison.update(
        reported_mean_min=observed["mean_min"],
        reported_variance_min2=observed["variance_min2"],
        mean_relative_error_limit=MEAN_LIMIT,
        variance_relative_error_limit=VARIANCE_LIMIT,
        variance_upper_bound_min2_if_pulse_support_is_7min=(
            primary["pulse_free_variance_min2"] + pulse_min**2 / 4
        ),
        **tracer_curve(output, samples[primary_name], pulse_min, evidence["water_flow_ml_min"],
                       evidence["reported_pulse_concentration_mg_ml"]["blue_dextran"]),
    )
    trace_void = fields[-1]["tracked_void_fraction"]
    flow_corrected = all(field["face_conservation"]["qualifies_local_correction"] for field in fields)
    recovered_variance = primary["pulse_free_relative_variance"]
    write_json(output / "summary.json", dict(
        method="OpenLB TRT water flow, flux-conserving face field, flux-weighted Brownian first passage",
        plan=plan,
        build_record_sha256=sha(output / "build.json"),
        geometry_plan_sha256=sha(output / "plan.json"),
        tracer_evidence_sha256=sha(SOURCE / "extraction.json"),
        flow={f"N{field['grid']}": {
            key: value for key, value in field.items() if key != "axial_plane_reflux_fractions"
        } for field in fields},
        mean_expected_from_void_inventory_min=(
            trace_void * plan["column_length_m"] / plan["superficial_speed_m_s"] / 60
        ),
        water_tracer=primary,
        primary_result=primary_name,
        finite_pulse_comparison=comparison,
        instantaneous_injection_sensitivity=check_moments(primary, observed, False),
        observed_azide_moments_separate_bead_access=evidence["reported_rtd_moments"]["sodium_azide"],
        sensitivity=records,
        numerical_comparisons=comparisons,
        face_shape_qualified=flow_corrected,
        water_tracer_benchmark_passed=(
            flow_corrected and comparison["passed"] and comparisons["mesh_at_matched_step"]["passed"]
            and comparisons["step_at_matched_mesh"]["passed"]
        ),
        open_start_peclet_moment_equivalent=2 / recovered_variance,
        transfer_scope=(
            "Synthetic bed, blue-dextran water benchmark only. Relative exit-time variance from "
            "recovered markers. Transferring to reacting liquid or different flow is a scenario. "
            "The reference porosity was inferred by the study from water-tracer dispersion."
        ),
        tracer_curve_sha256=sha(output / "tracer_curve.csv"),
    ))
    print("Tracer benchmark", comparison, flush=True)


if __name__ == "__main__":
    main()
