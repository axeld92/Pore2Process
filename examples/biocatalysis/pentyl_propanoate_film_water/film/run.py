"""Predict the film penalty across the frozen 20-run reactor comparison."""

from argparse import ArgumentParser
from datetime import datetime, timezone
from math import pi
from pathlib import Path

import numpy as np

from examples.biocatalysis.pentyl_propanoate.run import error_metrics, read_rows, write_rows
from examples.biocatalysis.pentyl_propanoate_transfer.published import run as published
from pore2process_demo.bridge import runtime
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT


HERE = Path(__file__).resolve().parent
PROPERTIES = HERE / "properties.json"
OUTPUT = ROOT / "results/biocatalysis/pentyl_propanoate_film_water/film"
CORRELATION_CASES = {
    "wilson_geankoplis_central": ("wilson_geankoplis", "central"),
    "wilson_geankoplis_slow": ("wilson_geankoplis", "low"),
    "wilson_geankoplis_fast": ("wilson_geankoplis", "high"),
    "corrected_printed_B2_central": ("corrected_printed_B2", "central"),
    "corrected_printed_B2_slow": ("corrected_printed_B2", "low"),
    "isolated_sphere_central": ("isolated_sphere", "central"),
}
NAMES = ("no_film", *CORRELATION_CASES)


def predict(julia, source, properties):
    kinetic_model = julia.LorenzoPackedBed
    film = julia.BeadFilm
    parameters = kinetic_model.PingPong(np.array([source["rate_parameters"][name]
                                                 for name in published.PARAMETERS]))
    bed = source["bed"]
    porosity, diameter = bed["bed_porosity"], bed["bead_mean_diameter_m"]
    area = pi * bed["diameter_m"]**2 / 4
    total_volume = area * bed["length_m"]
    liquid_volume = porosity * total_volume
    catalyst_mass = bed["catalyst_mass_kg"]
    if catalyst_mass != 0.006:
        raise ValueError("The frozen film comparison uses the reported 6 g of catalyst.")
    # The numerical value in kg/m³ is also the loading in g/L.
    catalyst_g_L = catalyst_mass / liquid_volume
    bead_surface = 6 * (1 - porosity) * total_volume / diameter
    area_per_liquid_volume = bead_surface / liquid_volume
    acid, alcohol = bed["feed_acid_mol_per_L"], bed["feed_alcohol_mol_per_L"]
    # Retain the transport convention of the previous reactor audit.
    axial_pe = bed["length_m"] * 0.20 / (porosity**2 * diameter)
    lookup = {round(row["temperature_K"] - 273.15): row for row in properties["conditions"]}
    predictions, diagnostics = [], []
    max_grid_change, max_acid_balance, max_film_balance = 0.0, 0.0, 0.0
    for temperature in bed["reaction_temperature_C"]:
        eta = source["effectiveness_by_temperature"][str(temperature)]["eta_internal"]
        liquid = lookup[temperature]
        rho, visc = liquid["density_kg_m3"], liquid["viscosity_Pa_s"]
        for flow in bed["feed_flow_ml_per_min"]:
            row = {"temperature_C": temperature, "flow_ml_min": flow, "catalyst_mass_kg": catalyst_mass}
            flow_m3_s = flow / 60e6
            u = flow_m3_s / area
            residence_min = liquid_volume * 1e6 / flow

            def solve(kinetics):
                return [kinetic_model.packed_bed(
                    kinetics, temperature, residence_min, catalyst_g_L, eta, axial_pe,
                    acid0=acid, alcohol0=alcohol, cells=cells,
                ) for cells in (512, 1024)]

            coarse, fine = solve(parameters)
            row["no_film"] = float(fine.conversion)
            max_grid_change = max(max_grid_change, float(abs(fine.conversion - coarse.conversion)))
            max_acid_balance = max(max_acid_balance, float(coarse.balance_error), float(fine.balance_error))
            for name, (correlation, mobility) in CORRELATION_CASES.items():
                molecular_d = liquid["diffusivity_m2_s"][mobility]
                coeff = film.film_coefficients(u, diameter, rho, visc, molecular_d, porosity, correlation)
                # 60 converts an inverse second to an inverse minute. C stays in mol/L.
                capacity = 60 * float(coeff.coefficient_m_s) * area_per_liquid_volume * acid
                restricted = film.FilmKinetics(parameters, capacity, eta)
                local = film.inlet_balance(restricted, temperature, catalyst_g_L, acid0=acid, alcohol0=alcohol)
                coarse, fine = solve(restricted)
                row[name] = float(fine.conversion)
                diagnostics.append({
                    "temperature_C": temperature,
                    "flow_ml_min": flow,
                    "case": name,
                    "viscosity_Pa_s": visc,
                    "density_kg_m3": rho,
                    "diffusivity_m2_s": molecular_d,
                    "Re": float(coeff.reynolds),
                    "Sc": float(coeff.schmidt),
                    "Sh": float(coeff.sherwood),
                    "k_f_m_s": float(coeff.coefficient_m_s),
                    "inlet_surface_acid_relative_to_bulk": float(local.surface_acid_relative_to_bulk),
                    "inlet_rate_relative_to_no_film": float(local.rate_relative_to_no_film),
                })
                max_grid_change = max(max_grid_change, float(abs(fine.conversion - coarse.conversion)))
                max_acid_balance = max(max_acid_balance, float(coarse.balance_error), float(fine.balance_error))
                max_film_balance = max(max_film_balance, float(local.film_balance_error_mol_L_min))
            predictions.append(row)
    if max_grid_change > 1e-4 or max_acid_balance > 1e-9 or max_film_balance > 1e-10:
        raise RuntimeError("Film or axial balance and grid checks failed.")
    return predictions, diagnostics, {
        "catalyst_mass_kg": catalyst_mass,
        "total_external_area_m2": bead_surface,
        "area_per_liquid_volume_m_inv": area_per_liquid_volume,
        "catalyst_g_L_void_liquid": catalyst_g_L,
        "maximum_mesh_change_absolute_conversion": max_grid_change,
        "maximum_normalized_axial_balance_error": max_acid_balance,
        "maximum_inlet_film_balance_error_mol_L_min": max_film_balance,
    }


def compare(predictions):
    measurements = read_rows(published.HOLDOUT)
    observed = {published.condition(row): float(row["acid_conversion"]) for row in measurements}
    if (len(predictions) != 20 or len(measurements) != 20 or len(observed) != 20
            or set(observed) != set(map(published.condition, predictions))):
        raise ValueError("The comparison requires the frozen 20-run reactor grid.")
    y = np.array([observed[published.condition(row)] for row in predictions])
    if not np.all(np.isfinite(y) & (y >= 0) & (y <= 1)):
        raise ValueError("Measured acid conversion must be a finite fraction in [0, 1].")
    baseline = np.array([row["no_film"] for row in predictions])
    results = [dict(row, observed_acid_conversion=float(value))
               for row, value in zip(predictions, y, strict=True)]
    scores = {}
    for case in NAMES:
        conversion = np.array([row[case] for row in predictions])
        scores[case] = error_metrics(conversion - y) | {
            "overpredicted_count": int(np.count_nonzero(conversion > y)),
            "max_conversion_decrease_from_no_film_percentage_points": float(100 * max(baseline - conversion)),
        }
    return results, scores


def run(output):
    source = read_json(published.MODEL_INPUTS)
    properties = read_json(PROPERTIES)
    julia = runtime()
    for name, path in (("LorenzoPackedBed", published.MODEL), ("BeadFilm", HERE / "model.jl")):
        if not julia.seval(f"isdefined(@__MODULE__, :{name})"):
            julia.include(str(path))
    predictions, diagnostics, summary = predict(julia, source, properties)
    output.mkdir(parents=True, exist_ok=False)
    write_rows(output / "predictions_before_comparison.csv", predictions)
    write_rows(output / "film_diagnostics.csv", diagnostics)
    predictions_digest = sha(output / "predictions_before_comparison.csv")
    results, comparison = compare(predictions)
    write_rows(output / "bed_comparison.csv", results)
    ranges = {}
    for case in CORRELATION_CASES:
        rows = [row for row in diagnostics if row["case"] == case]
        ranges[case] = {}
        for name in ("Re", "Sc", "Sh", "k_f_m_s", "inlet_rate_relative_to_no_film"):
            values = [row[name] for row in rows]
            ranges[case][name] = [min(values), max(values)]
    summary.update({
        "comparison": comparison,
        "liquid_cases": CORRELATION_CASES,
        "film_ranges": ranges,
        "sixty_C_three_ml_min": next(row for row in results if published.condition(row) == (60, 3)),
        "julia": str(julia.seval("VERSION")),
        "source_sha256": {
            str(path.relative_to(ROOT)): sha(path)
            for path in (PROPERTIES, published.MODEL_INPUTS, published.MODEL, published.HOLDOUT,
                         HERE / "model.jl", HERE / "run.py")
        },
        "predictions_before_comparison_sha256": predictions_digest,
        "bed_comparison_sha256": sha(output / "bed_comparison.csv"),
    })
    write_json(output / "summary.json", summary)
    return summary


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    parser.add_argument("--output", type=Path, default=OUTPUT / stamp)
    args = parser.parse_args()
    report = run(args.output.resolve())
    for name, metrics in report["comparison"].items():
        print(name, f'RMSE {metrics["rmse_percentage_points"]:.2f} pp,',
              f'maximum conversion drop {metrics["max_conversion_decrease_from_no_film_percentage_points"]:.2f} pp')
