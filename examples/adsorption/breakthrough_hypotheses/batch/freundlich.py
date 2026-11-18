"""Qualify the published Freundlich equilibrium with the existing Julia particles."""

import argparse
from datetime import datetime, timezone
from functools import lru_cache
from pathlib import Path
from uuid import uuid4

import numpy as np

from examples.adsorption.breakthrough_hypotheses.batch import run as langmuir
from examples.adsorption.film_particle import cases, data, fitting
from examples.adsorption.film_particle.bridge import batch, module
from examples.adsorption.transport_accuracy.run import validate_record
from pore2process_demo.bridge import runtime, runtime_versions
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.provenance import source_manifest

HERE = Path(__file__).resolve().parent
JULIA = HERE / "freundlich.jl"
LAW = dict(coefficient_mg_g_per_mg_L_to_exponent=24.6, exponent=1 / 2.74,
           source="Alhamed (2009), Table 1")


@lru_cache(maxsize=1)
def load():
    jl = runtime()
    jl.include(str(JULIA))
    return jl.PublishedFreundlich


def verify_law():
    model = load()
    K, a = LAW["coefficient_mg_g_per_mg_L_to_exponent"], LAW["exponent"]
    for concentration in (0., 1e-300, 1e-150, 1e-12, 1e-6, 0.006, .036, .15, .2, .4):
        expected = K * (1000 * concentration) ** a / 1000
        np.testing.assert_allclose(model.loading(concentration, model.TABLE1), expected, rtol=1e-12)
        for beta in (0., .01, 1., 30., 1e3, 1e6):
            rhs = concentration + beta * expected
            actual = model.surface_concentration(rhs, beta, model.TABLE1)
            np.testing.assert_allclose(actual, concentration, rtol=1e-12)
    if model.surface_concentration(1.4803e-319, 419.945, model.TABLE1) != 0:
        raise ValueError("A physically unresolvable concentration did not underflow to zero")


def predict_bath(inputs):
    module()
    result = batch(**dict(inputs, linear_partition_m3_kg=load().TABLE1))
    langmuir.verify_bath(result)
    return result


def fit(observations):
    @lru_cache(maxsize=None)
    def profile(log_film):
        def objective(log_D):
            errors = []
            for size in langmuir.TRAINING:
                points = observations[size]
                inputs = cases.batch_inputs(points["diameter_m"], 10 ** log_D, 10 ** log_film,
                                            radial_cells=12, dt_s=120.)
                result = predict_bath(inputs)
                errors.extend(langmuir.loading_at(result, points["time_s"]) - points["loading_kg_kg"])
            return float(np.mean(np.square(errors)))

        return fitting.golden(objective, -14., -9., iterations=24)

    grid = [-5., -4.5, -4., -3.5, -3.]
    best = min(range(len(grid)), key=lambda i: profile(grid[i])[1])
    film, _ = fitting.golden(lambda x: profile(x)[1], grid[max(0, best-1)], grid[min(len(grid)-1, best+1)],
                             iterations=16)
    film = min((film, *grid), key=lambda x: profile(x)[1])
    log_D, mse = profile(film)
    return dict(diffusivity_m2_s=10 ** log_D, batch_film_ref_m_s=10 ** film,
                training_RMSE_mg_g=mse ** .5 * 1000, film_log10_at_bound=film in (grid[0], grid[-1]),
                diffusivity_log10_at_bound=log_D in (-14., -9.))


def run():
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ") + "-" + uuid4().hex[:8]
    base = langmuir.OUTPUT / "freundlich_batch" / stamp
    base.mkdir(parents=True)
    print(base, flush=True)
    module()
    verify_law()
    observations = data.observations()
    equilibrium = langmuir.inventory_and_equilibrium(data.evidence())
    if equilibrium["maximum_endpoint_relative_error"]["Freundlich"] > .05:
        raise ValueError("Published Freundlich equilibrium failed the batch endpoint screen")
    parameters = fit(observations)
    print("Frozen Freundlich kinetics", parameters, flush=True)
    sizes = {key: dict(diameter_m=row["diameter_m"], role=row["role"]) for key, row in observations.items()}
    plan = dict(profiles={"freundlich": parameters}, law=LAW, calibration=list(langmuir.TRAINING),
                holdout=list(langmuir.WITHHELD), sizes=sizes, equilibrium=equilibrium)
    write_json(base / "plan.json", plan)
    artifacts = {}
    for name, size in sizes.items():
        for grid, (shells, step) in langmuir.RESOLUTIONS.items():
            inputs = cases.batch_inputs(size["diameter_m"], parameters["diffusivity_m2_s"],
                parameters["batch_film_ref_m_s"], radial_cells=shells, dt_s=step)
            path = base / f"freundlich_{name}_{grid}.json"
            write_json(path, dict(inputs=inputs, equilibrium=LAW, result=predict_bath(inputs)))
            artifacts[path.name] = sha(path)
    write_json(base / "frozen.json", dict(plan_sha256=sha(base / "plan.json"), curve_sha256=artifacts))
    scores = langmuir.compare(base, plan, observations)["freundlich"]
    files = [Path(__file__), JULIA, HERE / "freundlich_protocol.md", Path(langmuir.__file__), langmuir.PDF,
             data.PIXELS, langmuir.ROOT / "julia/adsorption/FilmParticle.jl"]
    write_json(base / "summary.json", dict(schema="pore2process.breakthrough-freundlich-batch.v1",
        qualification=scores, law=LAW, frozen_sha256=sha(base / "frozen.json"),
        source_sha256=source_manifest(files), versions=runtime_versions(), figure_parameters_fitted=0,
        artifacts_sha256={p.name: sha(p) for p in base.iterdir() if p.is_file()}))
    print("Qualification", {name: row["RMSE_mg_g"] for name, row in scores["scores"].items()},
          "passed", scores["qualifies"], flush=True)
    print(base / "summary.json", flush=True)


def check(base):
    summary = validate_record(base / "summary.json")
    plan = read_json(base / "plan.json")
    if (plan["law"] != LAW or summary["frozen_sha256"] != sha(base / "frozen.json") or
            summary["qualification"] != langmuir.compare(base, plan, data.observations())["freundlich"]):
        raise ValueError("Freundlich fit or held-out scores differ on replay")
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", type=Path)
    args = parser.parse_args()
    if args.check:
        print(check(args.check)["qualification"]["qualifies"])
    else:
        run()
