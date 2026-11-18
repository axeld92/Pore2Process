"""Qualify particle diffusivity profiles before scoring their column transfer."""

import argparse
from datetime import datetime, timezone
from pathlib import Path
from uuid import uuid4

import numpy as np

from examples.adsorption.alhamed.source import PDF
from examples.adsorption.film_particle import cases, data, fitting
from examples.adsorption.film_particle.bridge import batch, call
from examples.adsorption.film_particle.metrics import batch_metrics
from examples.adsorption.transport_accuracy.run import validate_record
from pore2process_demo.bridge import runtime_versions
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT
from pore2process_demo.provenance import source_manifest

HERE = Path(__file__).resolve().parent
OUTPUT = ROOT / "results/adsorption/breakthrough_hypotheses/batch"
ARCHIVE = ROOT / "examples/adsorption/film_particle/reference/film_particle_reference.json"
RESOLUTIONS = {"base": (24, 30.), "radial": (48, 30.), "time": (24, 15.)}
TRAINING = ("d1500", "d800")
WITHHELD = ("d450", "d225")


def inventory_and_equilibrium(source):
    """Audit the printed isotherms in their reported mg/L and mg/g units."""
    p = source["source_parameters"]
    qmax, b = p["capacity_kg_kg"] * 1000, p["affinity_m3_kg"] / 1000
    langmuir = lambda c: qmax * b * c / (1 + b * c)
    freundlich = lambda c: 24.6 * c ** (1 / 2.74)
    laws = {"Langmuir": langmuir, "Freundlich": freundlich}
    endpoints = {}
    for row in source["batch_table_3"]:
        initial, volume, carbon = (row[k] for k in ("inlet_mg_L", "solution_volume_L", "dose_g"))
        fitted = row["fitted_qe_mg_g"]
        values = {}
        for name, equilibrium in laws.items():
            lower, upper = 0., initial
            for _ in range(64):
                concentration = (lower + upper) / 2
                if volume * concentration + carbon * equilibrium(concentration) < volume * initial:
                    lower = concentration
                else:
                    upper = concentration
            concentration = (lower + upper) / 2
            loading = equilibrium(concentration)
            np.testing.assert_allclose(volume * concentration + carbon * loading, volume * initial, rtol=1e-12)
            values[name] = dict(concentration_mg_L=concentration, loading_mg_g=loading,
                                signed_relative_error=loading / fitted - 1)
        endpoints[str(row["experiment"])] = dict(initial_mg_L=initial, bath_L=volume, dry_carbon_g=carbon,
            fitted_endpoint_mg_g=fitted, models=values)
    area = source["figure_7"]["area_m2"]
    carbon = {str(height): p["bulk_density_kg_m3"] * area * height / 100 for height in (10, 15, 20, 25)}
    feed = source["figure_7"]["feed_kg_m3"] * 1000
    return dict(conditions=source["figure_2"], measured_dry_bed_mass_kg=None,
        density_basis="Table 1 packed bulk density multiplied by nominal bed volume",
        derived_mass_kg_by_bed_height_cm=carbon, column_area_m2=area,
        printed_bulk_density_kg_m3=p["bulk_density_kg_m3"], printed_voidage=p["porosity"],
        derived_particle_envelope_density_kg_m3=p["particle_density_kg_m3"],
        parameters=dict(langmuir_capacity_mg_g=qmax, langmuir_affinity_L_mg=b,
                        freundlich_coefficient=24.6, freundlich_n=2.74),
        endpoints=endpoints, maximum_endpoint_relative_error={name: max(
            abs(row["models"][name]["signed_relative_error"]) for row in endpoints.values()) for name in laws},
        column_feed_mg_L=feed, loading_at_feed_mg_g={name: fit(feed) for name, fit in laws.items()})


def verify_bath(result):
    if max(result[key] for key in ("normalized_mass_residual", "normalized_interface_flux_mismatch")) > 1e-8:
        raise ValueError("Finite-bath conservation or interface flux check failed")


def loading_at(result, times):
    return np.interp(times, result["time_s"], result["mean_loading_kg_kg"])


def training_profile(film, observations):
    def objective(log_diffusivity):
        errors = []
        for key in TRAINING:
            points = observations[key]
            inputs = cases.batch_inputs(points["diameter_m"], 10 ** log_diffusivity, film,
                                        radial_cells=12, dt_s=120.)
            result = batch(**inputs)
            verify_bath(result)
            error = loading_at(result, points["time_s"]) - points["loading_kg_kg"]
            errors.extend(error)
        return float(np.mean(np.square(errors)))

    log_diffusivity, mse = fitting.golden(objective, -14., -9., iterations=24)
    return dict(diffusivity_m2_s=10 ** log_diffusivity, batch_film_ref_m_s=film,
                training_RMSE_mg_g=mse ** .5 * 1000, origin="fixed film, training-only search")


def select_profiles(observations):
    archived = read_json(ARCHIVE)["batch_fit"]
    high_film = next(row for row in archived["kf_profiles"] if row["log10_batch_kf_ref"] == -3)
    fit = dict(diffusivity_m2_s=archived["diffusivity_m2_s"],
               batch_film_ref_m_s=archived["batch_film_ref_m_s"],
               training_RMSE_mg_g=archived["search_RMSE_mg_g"], origin="archived shared batch fit")
    slow = dict(diffusivity_m2_s=high_film["best_Ds_m2_s"], batch_film_ref_m_s=1e-3,
                training_RMSE_mg_g=high_film["RMSE_mg_g"], origin="archived high-film profile")
    return dict(primary=fit, slow=slow, fast=training_profile(4.5e-5, observations),
                edge=training_profile(4e-5, observations))


def predict(base, plan):
    artifacts = {}
    for name, profile in plan["profiles"].items():
        for diameter, size in plan["sizes"].items():
            for grid, (shells, dt) in RESOLUTIONS.items():
                inputs = cases.batch_inputs(size["diameter_m"], profile["diffusivity_m2_s"],
                                           profile["batch_film_ref_m_s"], radial_cells=shells, dt_s=dt)
                result = batch(**inputs)
                verify_bath(result)
                filename = f"{name}_{diameter}_{grid}.json"
                write_json(base / filename, dict(inputs=inputs, result=result))
                artifacts[filename] = sha(base / filename)
    write_json(base / "frozen.json", dict(plan_sha256=sha(base / "plan.json"), curve_sha256=artifacts))


def compare(base, plan, observations):
    frozen = read_json(base / "frozen.json")
    if frozen["plan_sha256"] != sha(base / "plan.json"):
        raise ValueError("Frozen batch inputs changed")
    equilibrium = plan["equilibrium"]["endpoints"]["1"]["models"]["Langmuir"]["loading_mg_g"] / 1000
    rows = {}
    for name, profile in plan["profiles"].items():
        scores = {}
        for diameter, values in observations.items():
            histories = {}
            for grid, (shells, dt) in RESOLUTIONS.items():
                filename = f"{name}_{diameter}_{grid}.json"
                if sha(base / filename) != frozen["curve_sha256"][filename]:
                    raise ValueError(f"Frozen batch history changed: {filename}")
                saved = read_json(base / filename)
                if saved["inputs"] != cases.batch_inputs(values["diameter_m"], profile["diffusivity_m2_s"],
                        profile["batch_film_ref_m_s"], radial_cells=shells, dt_s=dt):
                    raise ValueError(f"Batch inputs differ from the frozen plan: {filename}")
                verify_bath(saved["result"])
                histories[grid] = saved["result"]
            current = histories.pop("base")
            base_loading = loading_at(current, values["time_s"])
            refinements = {grid: float(np.max(abs(loading_at(r, values["time_s"]) - base_loading))) / equilibrium
                           for grid, r in histories.items()}
            metrics = batch_metrics(current, values)
            scores[diameter] = dict(**metrics, relative_refinements=refinements,
                                    numerical_passed=max(refinements.values()) <= .02)
        qualifies = all(row["RMSE_screen_passed"] and row["numerical_passed"] for row in scores.values())
        rows[name] = dict(**profile, scores=scores, qualifies=qualifies)
    return rows


def qualification(profiles):
    admissible = [name for name, row in profiles.items() if row["qualifies"]]
    if "primary" not in admissible or len(admissible) < 3:
        raise ValueError("The baseline and two particle-kinetic alternatives must pass the batch checks")
    by_diffusivity = sorted(admissible, key=lambda name: profiles[name]["diffusivity_m2_s"])
    return [by_diffusivity[0], by_diffusivity[-1]]


def run():
    start = datetime.now(timezone.utc)
    base = OUTPUT / "runs" / f"{start:%Y%m%dT%H%M%SZ}-{uuid4().hex[:8]}"
    base.mkdir(parents=True)
    print(base, flush=True)
    published = data.evidence()
    observations = data.observations()
    equilibrium = inventory_and_equilibrium(published)
    # Independent bracketed mass balance agrees with the existing Julia model.
    julia_q = call("bath_equilibrium")["loading_kg_kg"] * 1000
    np.testing.assert_allclose(equilibrium["endpoints"]["1"]["models"]["Langmuir"]["loading_mg_g"],
                               julia_q, rtol=1e-10)
    profiles = select_profiles(observations)
    sizes = {key: dict(diameter_m=row["diameter_m"], role=row["role"]) for key, row in observations.items()}
    if any(sizes[key]["role"] != "calibration" for key in TRAINING) or set(sizes) != set(TRAINING + WITHHELD):
        raise ValueError("Figure 2 training and holdout assignments changed")
    plan = dict(profiles=profiles, sizes=sizes, calibration=list(TRAINING), holdout=list(WITHHELD),
                equilibrium=equilibrium, source_pdf_sha256=sha(PDF), pixels_sha256=sha(data.PIXELS),
                archive_sha256=sha(ARCHIVE))
    write_json(base / "plan.json", plan)
    predict(base, plan)
    scores = compare(base, plan, observations)
    selected = qualification(scores)
    modules = [HERE / "protocol.md", Path(__file__), ARCHIVE, data.PIXELS, PDF,
               ROOT / "julia/adsorption/FilmParticle.jl", ROOT / "uv.lock", ROOT / "julia/Manifest.toml",
               *(Path(m.__file__) for m in (cases, data, fitting))]
    write_json(base / "summary.json", dict(schema="pore2process.breakthrough-particle-batch.v1",
        qualification=scores, selected=selected, frozen_sha256=sha(base / "frozen.json"),
        source_sha256=source_manifest(modules), versions=runtime_versions(), figure7_fitted_parameters=0,
        artifacts_sha256={p.name: sha(p) for p in base.iterdir() if p.is_file()}))
    print({name: (scores[name]["diffusivity_m2_s"], scores[name]["qualifies"]) for name in scores}, flush=True)
    print(base / "summary.json", flush=True)


def check(base):
    summary = validate_record(base / "summary.json")
    files = (base / "frozen.json", data.PIXELS, PDF, ARCHIVE)
    plan = read_json(base / "plan.json")
    digests = (summary["frozen_sha256"], plan["pixels_sha256"], plan["source_pdf_sha256"], plan["archive_sha256"])
    if any(sha(path) != digest for path, digest in zip(files, digests)):
        raise ValueError("Saved batch qualification inputs changed")
    expected = compare(base, plan, data.observations())
    if summary["qualification"] != expected or summary["selected"] != qualification(expected):
        raise ValueError("Batch scores or selected transfer profiles changed on replay")
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", type=Path, help="Replay the numerical and evidence checks on a frozen run")
    args = parser.parse_args()
    if args.check:
        checked = check(args.check)
        print(dict(sources=len(checked["source_sha256"]), artifacts=len(checked["artifacts_sha256"])))
    else:
        run()
