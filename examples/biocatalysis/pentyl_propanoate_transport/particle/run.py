"""Test bead transport scenarios at full catalyst mass on the frozen reactor grid."""

from argparse import ArgumentParser
from datetime import datetime, timezone
from math import log, sqrt
from pathlib import Path

import numpy as np

from examples.biocatalysis.pentyl_propanoate.run import (
    HERE as BASE,
    HOLDOUT,
    PROTOCOL,
    error_metrics,
    read_rows,
    write_rows,
)
from examples.biocatalysis.pentyl_propanoate_transfer.published.run import (
    HERE as PUBLISHED,
    MODEL as PUBLISHED_MODEL,
    MODEL_INPUTS as PUBLISHED_INPUTS,
    PARAMETERS,
)
from examples.biocatalysis.pentyl_propanoate_transport.phase import run as liquid
from pore2process_demo.bridge import runtime
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT


HERE = Path(__file__).resolve().parent
EVIDENCE = HERE / "evidence.json"
FROZEN = BASE / "reference/summary.json"
APPARENT_MODEL = BASE / "model.jl"
TRANSFER = BASE.with_name("pentyl_propanoate_transfer")
MIXING = TRANSFER / "mixing/reference/summary.json"
PUBLISHED_REFERENCE = PUBLISHED / "reference_summary.json"
OUTPUT = ROOT / "results/biocatalysis/pentyl_propanoate_transport/particle"


def bead_scenarios(evidence):
    study = evidence["lorenzo_2025"]
    porosity = study["particle_porosity"]
    # Compare the two measured-property estimates at fixed tortuosity,
    # constrictivity and molecular diffusion.
    diffusivity = porosity["dry_nitrogen_sorption"]["value"] / porosity["aqueous_tracers"]["value"]
    dry = porosity["dry_nitrogen_sorption"]
    solid_volume_cm3_g = 1000 / dry["material_density_kg_m3"]
    dry_skeletal_basis = dry["pore_volume_cm3_g"] / (dry["pore_volume_cm3_g"] + solid_volume_cm3_g)
    groups = evidence["chaussard_2024"]["diameter_distributions"]
    sieved = max((row for row in groups if row["fraction"] != "unsieved"),
                 key=lambda row: row["wet_to_dry_mean_ratio"])
    diameter = 1 / sieved["wet_to_dry_mean_ratio"]
    cases = {
        "reference": {"D_ratio": 1.0, "diameter_ratio": 1.0, "inactive_core_ratio": 0.0},
        "porosity_contrast": {"D_ratio": diffusivity, "diameter_ratio": 1.0,
                              "inactive_core_ratio": 0.0},
        "porosity_and_diameter": {"D_ratio": diffusivity, "diameter_ratio": diameter,
                                  "inactive_core_ratio": 0.0},
        "skeletal_porosity_and_diameter": {
            "D_ratio": dry_skeletal_basis / porosity["aqueous_tracers"]["value"],
            "diameter_ratio": diameter,
            "inactive_core_ratio": 0.0,
        },
    }
    for depth in evidence["mei_2003"]["shell_depth_um"]:
        cases[f"shell_{depth}_um"] = {
            "D_ratio": diffusivity,
            "diameter_ratio": diameter,
            "inactive_core_ratio": 1 - 2 * depth / study["particle_mean_diameter_um"],
        }
    return cases, sieved["fraction"], dry_skeletal_basis, sqrt(diffusivity * diameter)


def multipliers(particle, reference_eta, cases, strong_limit):
    modulus = float(particle.modulus_for(reference_eta))
    values = {
        name: float(particle.relative_effectiveness(
            modulus, case["D_ratio"], case["diameter_ratio"], core=case["inactive_core_ratio"],
        ))
        for name, case in cases.items()
    }
    # Uniform-sphere relative effectiveness tends to sqrt(D_ratio*d_ratio)
    # at high modulus. This deliberately extreme limit is independent of T.
    values["strong_diffusion_limit"] = strong_limit
    return modulus, values


def predict(julia, plan, evidence, published, frozen, mixing, published_reference):
    particle = julia.ParticleTransport
    apparent = julia.PentylPropanoate
    authors = julia.LorenzoPackedBed
    p = authors.PingPong(np.array([published["rate_parameters"][k] for k in PARAMETERS]))
    geometries, selected_group, dry_skeletal_basis, limit = bead_scenarios(evidence)
    mass = plan["bed_catalyst_mass_kg"]
    acid = plan["bed_feed_acid_mol_L"]
    alcohol = plan["bed_feed_alcohol_mol_L"]
    published_numerics = published_reference["numerics"]
    liquid_ml = published_numerics["void_liquid_volume_ml"]
    catalyst_g_L = mass * 1e6 / liquid_ml
    if mass != 0.006 or 1000 * mass != published_numerics["nominal_catalyst_mass_g"]:
        raise ValueError("A different catalyst mass would change the frozen comparison.")
    water_pe = mixing["conditional_danckwerts_transfer"]["nominal_rectangular_pulse"]["danckwerts_peclet"]
    published_pe = published_numerics["low_re_limit_peclet"]
    mesh_levels = list(zip(plan["axial_dispersion_cell_counts"],
                           published_numerics["finite_volume_refinement"], strict=True))
    exposure_convention = published["bed"]["bed_porosity"]
    eta_by_temperature = published["effectiveness_by_temperature"]
    by_temperature, rows = {}, []
    maximum_grid_change, maximum_balance = 0.0, 0.0
    for temperature in plan["bed_controlled_temperatures_C"]:
        key = str(temperature)
        reference_eta = eta_by_temperature[key]["eta_internal"]
        modulus, factors = multipliers(particle, reference_eta, geometries, limit)
        target_modulus = float(particle.modulus_for(exposure_convention * reference_eta))
        by_temperature[key] = {
            "reference_eta_from_Table_B1": reference_eta,
            "reference_modulus": modulus,
            "relative_rate_multipliers": factors,
            "D_ratio_for_volume_convention_at_fixed_diameter": (modulus / target_modulus)**2,
        }
        log_k, slope, log_K, K_slope = frozen["batch"]["primary_fit_by_temperature"][key]["coefficients"]
        for flow in plan["bed_total_feed_flow_ml_min"]:
            row = {"temperature_C": temperature, "flow_ml_min": flow, "catalyst_mass_kg": mass}
            for name, factor in factors.items():
                kinetics = apparent.Kinetics(log_k + log(factor), slope, log_K, K_slope)
                both_grids = []
                for apparent_cells, published_cells in mesh_levels:
                    first = apparent.bed_prediction(
                        kinetics, temperature, mass, acid * flow / 60_000,
                        peclet=water_pe, cells=apparent_cells,
                    )
                    second = authors.packed_bed(
                        p, temperature, liquid_ml / flow, catalyst_g_L,
                        reference_eta * factor, published_pe, acid0=acid, alcohol0=alcohol,
                        cells=published_cells,
                    )
                    maximum_balance = max(maximum_balance, float(first.balance_error), float(second.balance_error))
                    both_grids.append(np.array([float(first.conversion), float(second.conversion)]))
                coarse, fine = both_grids
                maximum_grid_change = max(maximum_grid_change, float(np.max(np.abs(coarse - fine))))
                row[f"batch_2026_{name}"], row[f"published_2025_{name}"] = map(float, fine)
            rows.append(row)
    if (maximum_grid_change > plan["max_grid_change_absolute_conversion"]
            or maximum_balance > plan["max_stoichiometric_balance_error"]):
        raise RuntimeError("Grid refinement or acid balance failed on the frozen reactor cases.")
    return rows, {
        "bead_geometry_and_mobility": geometries,
        "largest_sieved_diameter_contrast_group": selected_group,
        "dry_particle_porosity_if_density_is_skeletal": dry_skeletal_basis,
        "uniform_sphere_limit_rate_multiplier": limit,
        "volume_convention_rate_multiplier_for_scale_comparison": exposure_convention,
        "temperature_cases": by_temperature,
        "maximum_grid_change_absolute_conversion": maximum_grid_change,
        "maximum_normalized_flux_balance_error": maximum_balance,
    }


def compare(predictions):
    def key(row):
        return int(row["temperature_C"]), float(row["flow_ml_min"])

    measured = read_rows(HOLDOUT)
    observations = {key(row): float(row["acid_conversion"]) for row in measured}
    if (len(observations) != 20 or len(measured) != 20 or len(predictions) != 20
            or set(observations) != set(map(key, predictions))):
        raise ValueError("The sealed predictions and held-out runs must cover the same grid.")
    observed = np.array([observations[key(row)] for row in predictions])
    if not np.all((0 <= observed) & (observed <= 1)):
        raise ValueError("Held-out conversion fractions are invalid.")
    results = [dict(row, measured_acid_conversion=y)
               for row, y in zip(predictions, observed, strict=True)]
    metrics = {}
    for model in ("batch_2026", "published_2025"):
        baseline = np.array([row[f"{model}_reference"] for row in predictions])
        metrics[model] = {}
        for name in predictions[0]:
            if name.startswith(f"{model}_"):
                values = np.array([row[name] for row in predictions])
                errors = values - observed
                metrics[model][name.removeprefix(f"{model}_")] = error_metrics(errors) | {
                    "overpredicted_count": int(np.count_nonzero(errors > 0)),
                    "maximum_reduction_from_reference_percentage_points": float(100 * max(baseline - values)),
                }
    return results, metrics


def run(output):
    plan = read_json(PROTOCOL)
    evidence = read_json(EVIDENCE)
    published = read_json(PUBLISHED_INPUTS)
    frozen = read_json(FROZEN)
    mixing = read_json(MIXING)
    published_reference = read_json(PUBLISHED_REFERENCE)
    julia = runtime()
    for name, source in (("ParticleTransport", HERE / "model.jl"),
                         ("PentylPropanoate", APPARENT_MODEL),
                         ("LorenzoPackedBed", PUBLISHED_MODEL)):
        if not julia.seval(f"isdefined(@__MODULE__, :{name})"):
            julia.include(str(source))
    predictions, diagnostics = predict(julia, plan, evidence, published, frozen, mixing, published_reference)
    properties = read_json(liquid.PROPERTIES)
    for row in predictions:
        for branch in ("batch_2026", "published_2025"):
            for basis in ("porosity", "skeletal_porosity"):
                conversion = row[f"{branch}_{basis}_and_diameter"]
                combined = liquid.readout_scenarios(conversion, published, properties)
                row[f"{branch}_{basis}_diameter_and_omit_water"] = combined["omit_generated_water"]
    output.mkdir(parents=True, exist_ok=False)
    write_rows(output / "predictions_before_comparison.csv", predictions)
    results, metrics = compare(predictions)
    write_rows(output / "bed_comparison.csv", results)
    summary = {
        "julia": str(julia.seval("VERSION")),
        "catalyst_mass_kg": plan["bed_catalyst_mass_kg"],
        "source_sha256": {
            str(path.relative_to(ROOT)): sha(path) for path in (
                EVIDENCE, PROTOCOL, FROZEN, APPARENT_MODEL, PUBLISHED_INPUTS, PUBLISHED_MODEL,
                PUBLISHED_REFERENCE, MIXING, HOLDOUT, liquid.PROPERTIES, Path(liquid.__file__).resolve(),
                HERE / "model.jl", HERE / "run.py",
            )
        },
        **diagnostics,
        "comparison": metrics,
        "bed_comparison_sha256": sha(output / "bed_comparison.csv"),
    }
    write_json(output / "summary.json", summary)
    return summary


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    parser.add_argument("--output", type=Path, default=OUTPUT / stamp)
    args = parser.parse_args()
    record = run(args.output.resolve())
    for model, scenarios in record["comparison"].items():
        print(model, {name: round(row["rmse_percentage_points"], 3) for name, row in scenarios.items()})
