"""Replay the published 2025 reactor model under explicit basis assumptions."""

import argparse
from datetime import datetime, timezone
from math import pi
from pathlib import Path

import numpy as np

from examples.biocatalysis.pentyl_propanoate.run import error_metrics, read_rows, write_rows
from pore2process_demo.bridge import runtime
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT


HERE = Path(__file__).resolve().parent
SOURCE = ROOT / "literature/sources/lorenzo2025"
MODEL_INPUTS = SOURCE / "author_model.json"
AUTHOR = SOURCE / "author_conversion_2to1.csv"
HOLDOUT = SOURCE / "conversion_2to1.csv"
BATCH = SOURCE / "batch_2to1_markers.csv"
MODEL = HERE / "model.jl"
OUTPUT = ROOT / "results/biocatalysis/pentyl_propanoate_transfer/published"
PARAMETERS = (
    "k1_ref", "k1_exp", "k2_ref", "k2_exp", "k3", "k4", "k5", "k6", "k7", "k8",
    "k9", "k10", "Kd_ref", "Kd1", "Kd2", "Kd_exp", "T_ref_K",
)


def condition(row):
    return int(row["temperature_C"]), float(row["flow_ml_min"])


def score(predictions, path, cases):
    measured = read_rows(path)
    observations = {condition(row): float(row["acid_conversion"]) for row in measured}
    if len(observations) != len(measured) or set(observations) != set(map(condition, predictions)):
        raise ValueError(f"{path.name} does not cover the frozen case grid.")
    observed = np.array([observations[condition(row)] for row in predictions])
    if not np.all((0 <= observed) & (observed <= 1)):
        raise ValueError(f"{path.name} contains an invalid conversion.")
    metrics = {
        name: error_metrics(np.array([row[name] for row in predictions]) - observed)
        for name in cases
    }
    return observed, metrics


def predict(model, p, source):
    bed = source["bed"]
    eps = bed["bed_porosity"]
    volume_ml = pi * bed["diameter_m"]**2 / 4 * bed["length_m"] * 1e6
    liquid_ml = eps * volume_ml
    mass_g = 1000 * bed["catalyst_mass_kg"]
    bulk_g_L = 1000 * mass_g / volume_ml
    liquid_g_L = 1000 * mass_g / liquid_ml
    acid, alcohol = bed["feed_acid_mol_per_L"], bed["feed_alcohol_mol_per_L"]
    # Eq. 3 combined with u/eps in Eq. 2 gives this Pe on z/L.
    pe = bed["length_m"] / (eps**2 * bed["bead_mean_diameter_m"])
    pe_low_re = pe * 0.20
    pe_re_one = pe * (0.20 + 0.011)
    cases = (
        "bulk_volume_eq4",
        "void_liquid_eq4",
        "void_liquid_unity_effectiveness_eq4",
        "bulk_volume_flux_inlet",
        "bulk_volume_re_one",
        "bulk_volume_temperature_ester_term",
    )
    predictions = []
    max_grid_change, max_balance_error = 0.0, 0.0
    for temperature in bed["reaction_temperature_C"]:
        eta = source["effectiveness_by_temperature"][str(temperature)]["eta_internal"]
        for flow in bed["feed_flow_ml_per_min"]:
            tau = liquid_ml / flow
            settings = (
                (bulk_g_L, eta, pe_low_re, False, False),
                (liquid_g_L, eta, pe_low_re, False, False),
                (liquid_g_L, 1.0, pe_low_re, False, False),
                (bulk_g_L, eta, pe_low_re, True, False),
                (bulk_g_L, eta, pe_re_one, False, False),
                (bulk_g_L, eta, pe_low_re, False, True),
            )
            row = {"temperature_C": temperature, "flow_ml_min": flow, "residence_min": tau}
            for name, (concentration, effectiveness, peclet, flux, dependent) in zip(
                cases, settings, strict=True
            ):
                parameters = {"acid0": acid, "alcohol0": alcohol, "inlet_flux": flux,
                              "temperature_dependent_ester_term": dependent}
                coarse = model.packed_bed(p, temperature, tau, concentration, effectiveness, peclet,
                                          cells=512, **parameters)
                fine = model.packed_bed(p, temperature, tau, concentration, effectiveness, peclet,
                                        cells=1024, **parameters)
                row[name] = float(fine.conversion)
                max_grid_change = max(max_grid_change, abs(float(coarse.conversion) - row[name]))
                max_balance_error = max(max_balance_error, float(coarse.balance_error),
                                        float(fine.balance_error))
            predictions.append(row)
    if max_grid_change > 1e-4 or max_balance_error > 1e-9:
        raise RuntimeError("The reaction grid or flux-balance check failed.")
    return predictions, cases, {
        "nominal_catalyst_mass_g": mass_g,
        "bed_volume_ml": volume_ml,
        "void_liquid_volume_ml": liquid_ml,
        "bulk_volume_catalyst_concentration_g_L": bulk_g_L,
        "void_liquid_catalyst_concentration_g_L": liquid_g_L,
        "bulk_volume_effective_mass_on_liquid_basis_g": eps * mass_g,
        "inlet_acid_disappearance_at_40C_25_3g_per_L_mol_L_min": float(
            model.disappearance_rate(p, 40, 0, 25.3, acid0=acid, alcohol0=alcohol)
        ),
        "low_re_limit_peclet": pe_low_re,
        "re_one_peclet": pe_re_one,
        "finite_volume_refinement": [512, 1024],
        "max_grid_change_absolute_conversion": max_grid_change,
        "max_normalized_stoichiometric_balance_error": max_balance_error,
    }


def batch_basis_check(model, p, source):
    # The published axes use hours while Table S3 leaves the rate units blank.
    # Compare both interpretations with the batch markers using the printed dose.
    acid = source["bed"]["feed_acid_mol_per_L"]
    alcohol = source["bed"]["feed_alcohol_mol_per_L"]
    printed_dose = 25.3
    nominal_mass_fraction = source["rate_basis"]["batch_figure_S3_catalyst_mass_fraction_relative_to_initial_acid"]
    nominal_dose = acid * nominal_mass_fraction * 74.0785
    settings = {
        "minute_basis_printed_dose": (1.0, printed_dose),
        "hour_basis_printed_dose": (1 / 60, printed_dose),
        "minute_basis_nominal_mass_fraction": (1.0, nominal_dose),
    }
    rows = []
    for marker in read_rows(BATCH):
        temperature, time, observed = (
            float(marker[key]) for key in ("temperature_C", "time_min", "acid_conversion")
        )
        if time == 0:
            continue
        predictions = {
            name: float(model.batch_prediction(p, temperature, time * scale, dose,
                                               acid0=acid, alcohol0=alcohol))
            for name, (scale, dose) in settings.items()
        }
        rows.append({**marker, **predictions})
    observed = np.array([float(row["acid_conversion"]) for row in rows])
    metrics = {
        name: {"catalyst_g_L": dose,
               **error_metrics(np.array([row[name] for row in rows]) - observed)}
        for name, (_, dose) in settings.items()
    }
    return rows, metrics


def run(output):
    output.mkdir(parents=True, exist_ok=False)
    source = read_json(MODEL_INPUTS)
    jl = runtime()
    jl.include(str(MODEL))
    model = jl.LorenzoPackedBed
    p = model.PingPong(np.array([source["rate_parameters"][name] for name in PARAMETERS]))
    rows, cases, numerics = predict(model, p, source)
    # Freeze the full prediction grid before either comparison is read.
    write_rows(output / "predictions_before_comparison.csv", rows)
    author, against_author = score(rows, AUTHOR, cases)
    observed, against_observation = score(rows, HOLDOUT, cases)
    comparison = [
        {**row, "author_figure_6c": float(fig), "observed": float(obs)}
        for row, fig, obs in zip(rows, author, observed, strict=True)
    ]
    write_rows(output / "comparison.csv", comparison)
    batch, against_batch = batch_basis_check(model, p, source)
    write_rows(output / "batch_basis_comparison.csv", batch)
    summary = {
        "scope": "Conditional acid-sink reconstruction. Eq. 2 with the Eq. S1 disappearance sign generates acid.",
        "rate_units_assumed": "mol/L/min with concentrations mol/L and catalyst g/L",
        "source_sha256": {
            str(path.relative_to(ROOT)): sha(path)
            for path in (MODEL_INPUTS, AUTHOR, HOLDOUT, BATCH, MODEL, HERE / "run.py")
        },
        "numerics": {"julia": str(jl.seval("VERSION")), **numerics},
        "batch_basis_check_on_2025_kinetic_data": against_batch,
        "comparison": {
            "scenario_against_author_figure_6c": against_author,
            "scenario_against_held_out_measurements": against_observation,
            "author_figure_6c_against_held_out_measurements": error_metrics(author - observed),
        },
        "max_shift_from_bulk_eq4_percentage_points": {
            name: 100 * max(abs(row[name] - row["bulk_volume_eq4"]) for row in rows)
            for name in cases[1:]
        },
        "comparison_sha256": sha(output / "comparison.csv"),
    }
    write_json(output / "summary.json", summary)
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    default = OUTPUT / datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    parser.add_argument("--output", type=Path, default=default)
    args = parser.parse_args()
    record = run(args.output.resolve())
    print(args.output.resolve() / "summary.json")
    print(record["comparison"]["scenario_against_author_figure_6c"]["bulk_volume_eq4"])
