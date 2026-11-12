"""Transfer the qualified survivor model to both shared-carbon columns."""
import argparse
from datetime import datetime, timezone
from pathlib import Path
import shutil
from uuid import uuid4

import numpy as np

from examples.adsorption.alhamed.study import compare
from examples.adsorption.contact_column.inputs import CONTACT, CONTACT_SHA, EVENT, shares_from_events
from examples.adsorption.contact_column.run import SHORT_S, coordinates, crossings, scoring_curve, startup_refinement
from examples.adsorption.carbon_contact.geometry import wall_faces
from examples.adsorption.exchange_column.run import CASES, PIXELS, source_inputs
from examples.adsorption.film_particle.metrics import column_passed
from examples.adsorption.porous_face_handoff.run import archived
from examples.adsorption.shared_contact.run import LEVELS, STARTUP, full_refinement, review
from examples.adsorption.transport_accuracy.run import validate_record
from examples.adsorption.wall_uptake import kernel
from examples.adsorption.wall_uptake.estimate import coefficients
from examples.adsorption.wall_uptake.run import conservation
from pore2process_demo.bridge import runtime_versions
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT
from pore2process_demo.provenance import source_manifest

from .predict import states
from .run import ARCHIVE, HERE, OUTPUT, fields, inputs

HOLDOUT = OUTPUT/"runs/20261112T150521Z-cc8dc4c0/summary.json"
VARIANTS = dict(quartiles=("primary", "twin_quartiles"),
                parent_only=("primary", "parent_only"),
                half_step=("half_step", "twin_quartiles"))


def prepare(holdout):
    qualified = validate_record(holdout)
    if not qualified["scorecard"]["wall_holdout_passed"]:
        raise ValueError("The saved wall holdout failed")
    frozen_path = holdout.parent/"frozen_prediction.json"
    if sha(frozen_path) != qualified["frozen_prediction_sha256"]:
        raise ValueError("Frozen holdout changed")
    frozen = read_json(frozen_path)
    flow, twins, mass, _ = inputs()
    old = read_json(ARCHIVE/"model.json")
    if frozen["dry_carbon_mass_kg"] != mass or old["dry_carbon_mass_kg"] != mass:
        raise ValueError("Carbon inventory changed")
    if sha(CONTACT/"summary.json") != CONTACT_SHA:
        raise ValueError("Pinned wall contact archive changed")
    validate_record(CONTACT/"summary.json")
    info, ijk, _ = archived(96)
    wall = wall_faces(info, ijk)
    delays = np.array(old["closure"]["normalized_transit_times"])
    cases = {}
    for size, p in old["sizes"].items():
        cases[size] = {}
        for name, (step, rule) in VARIANTS.items():
            key = f"{size}_{step}"
            data = fields(key, len(twins))
            model = frozen["cases"][key]["models"][rule]
            groups = model["subgroups"]
            labels = states(data[:, :2], twins, flow, groups)
            fine_flow = np.array(model["flow"])
            joint = np.array(model["carrier_flux"])
            mix = joint.T/joint.sum(axis=0)[:, None]
            fine_delays = np.repeat(delays, groups)
            parent_flux = (mix*fine_flow[:, None]).T.reshape(len(flow), groups, len(flow), groups)
            np.testing.assert_allclose(parent_flux.sum(axis=(1, 3)),
                                       np.array(p["mix"]).T*flow[None, :], rtol=2e-11, atol=0.)
            np.testing.assert_allclose((fine_flow*fine_delays).reshape(len(flow), groups).sum(axis=1),
                                       flow*delays, rtol=1e-13, atol=0.)
            np.testing.assert_allclose(mix.T@fine_flow, fine_flow, rtol=2e-11, atol=0.)
            # Only chronological first-period reflections contribute to sector area.
            events = np.fromfile(CONTACT/(key+".f64.faces"), dtype=EVENT)
            first = events[events["marker_phase"] > 0]
            marker = first["marker_phase"].astype(np.int64)-1
            np.testing.assert_array_equal(np.bincount(marker, minlength=len(twins)) > 0, np.isfinite(data[:, 2]))
            shares = shares_from_events(first, wall, marker, 0, labels, len(fine_flow))
            s = np.array(model["fresh_survival"])
            all_hit = s == 0
            count = np.bincount(labels[:, 0], minlength=len(s))
            finite_s = np.where(all_hit, .5/count, s)
            sink_two = float(finite_s @ joint @ finite_s)
            film = coefficients(-np.log(finite_s), shares, fine_flow, p["parameters"], p["trace"])
            cases[size][name] = dict(closure=dict(flow_fractions=fine_flow.tolist(),
                                                  normalized_transit_times=fine_delays.tolist()),
                                     mix=mix.tolist(), weights=shares.tolist(), trace=p["trace"],
                                     parameters=p["parameters"], coefficients=film,
                                     parent_mix=p["mix"], all_hit_count=int(all_hit.sum()),
                                     all_hit_flow_share=float(fine_flow[all_hit].sum()),
                                     survival_floor_change=sink_two-model["two_period_survival"])
    return dict(holdout_sha256=sha(holdout), dry_carbon_mass_kg=mass, cases=cases)


def result_at(base, files):
    return read_json(base/files[-1])["result"] if files else None


def refine_column(base, generated, size, name, setup, parent_resolution=None):
    p0 = setup["parameters"]
    feed = p0["inlet_concentration_kg_m3"]
    files = dict(full=[], short=[])

    def save(p, phase):
        r = kernel.solve(generated, setup["closure"], setup["mix"], setup["weights"],
                         setup["coefficients"]["path_film_m_s"], p)
        conservation(r, p, setup)
        filename = f"{size}_{name}_{phase}{len(files[phase])}.json"
        write_json(base/filename, dict(inputs=p, result=r))
        files[phase].append(filename)
        print(size, name, filename, crossings(r, feed), "runtime s", round(r["total_s"], 2), flush=True)
        return r

    if parent_resolution is not None:
        for phase in ("full", "short"):
            selected = parent_resolution[phase]
            if selected:
                p = read_json(base/selected[-1])["inputs"]
                save(p, phase)
    else:
        previous = None
        for subcells, radial_cells, dt_s, startup_dt_s in LEVELS:
            p = {**p0, **dict(subcells=subcells, radial_cells=radial_cells, dt_s=dt_s, startup_dt_s=startup_dt_s)}
            full = save(p, "full")
            if previous is not None:
                changes = full_refinement(previous, full, feed)
                if changes["passed"]:
                    break
            previous = full
        if not changes["startup"]["passed"]:
            previous = full
            for subcells, radial_cells, startup_dt_s in STARTUP:
                short = save({**p, **dict(subcells=subcells, radial_cells=radial_cells,
                             startup_dt_s=startup_dt_s, startup_s=SHORT_S, t_end_s=SHORT_S)}, "short")
                if startup_refinement(previous, short, feed)["passed"]:
                    break
                previous = short
        files["review"] = review(base, files)
    return files


def curve(base, files):
    return scoring_curve(result_at(base, files["full"]), result_at(base, files["short"]))


def score(base, frozen, prepared):
    _, source, _, _ = source_inputs()
    if sha(PIXELS) != source["pixel_csv_sha256"]:
        raise ValueError("Figure 7 markers changed")
    cases = {}
    for size, (series, _) in CASES.items():
        setup = prepared["cases"][size]
        files = frozen["cases"][size]
        feed = setup["quartiles"]["parameters"]["inlet_concentration_kg_m3"]
        cases[size] = {}
        for name, item in files.items():
            r = curve(base, item)
            metrics = compare(r, series)
            cases[size][name] = dict(metrics=metrics, physical_screen_passed=column_passed(metrics),
                                     crossing_s=crossings(r, feed),
                                     normalized_mass_residual=r["normalized_mass_residual"],
                                     normalized_interface_flux_mismatch=r["normalized_interface_flux_mismatch"])
        t, c = coordinates(curve(base, files["quartiles"]), feed)
        tf, cf = coordinates(curve(base, files["half_step"]), feed)
        cases[size]["wall_step_sensitivity_max_fraction"] = float(max(abs(c-np.interp(t, tf, cf))))
    return cases


def run(holdout):
    prepared = prepare(holdout)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")+"-"+uuid4().hex[:8]
    base = OUTPUT/"column_runs"/stamp
    base.mkdir(parents=True)
    sources = source_manifest([HERE/"column.py", HERE/"column_protocol.md", HERE/"predict.py",
                               HERE.parent/"wall_uptake/kernel.py", kernel.ORIGINAL,
                               kernel.JULIA/"ExchangeColumn.jl", kernel.JULIA/"FilmParticle.jl",
                               HERE.parent/"alhamed/study.py"])
    for name in sources:
        copy = base/"source_snapshot"/name
        copy.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT/name, copy)
    write_json(base/"model.json", prepared)
    generated = base/"FilmByPathColumn.jl"
    generated.write_text(kernel.source())
    checks = kernel.verify(generated)
    frozen = dict(model_sha256=sha(base/"model.json"), cases={})
    for size, cases in prepared["cases"].items():
        item = {}
        for name, setup in cases.items():
            resolution = item["quartiles"] if name == "half_step" else None
            item[name] = refine_column(base, generated, size, name, setup, resolution)
        frozen["cases"][size] = item
    write_json(base/"frozen_predictions.json", frozen)
    digest = sha(base/"frozen_predictions.json")
    print("Column curves frozen before Figure 7 scoring:", digest, flush=True)
    scored = score(base, frozen, prepared)
    summary = base/"summary.json"
    write_json(summary, dict(schema="pore2process.survivor-exposure-column.v1", run_id=stamp,
                             holdout=str(holdout.relative_to(ROOT)), source_sha256=sources,
                             artifacts_sha256={str(p.relative_to(base)): sha(p)
                                               for p in base.rglob("*") if p.is_file()},
                             frozen_predictions_sha256=digest, figure7=scored, solver_checks=checks,
                             figure7_fitted_parameters=0, versions=runtime_versions()))
    write_json(OUTPUT/"latest_column.json", dict(summary=str(summary.relative_to(ROOT))))
    print("Completed Figure 7 comparisons:", summary, flush=True)


def check(summary):
    record = validate_record(summary)
    base = summary.parent
    prepared = prepare(ROOT/record["holdout"])
    frozen = read_json(base/"frozen_predictions.json")
    if (sha(base/"frozen_predictions.json") != record["frozen_predictions_sha256"]
            or read_json(base/"model.json") != prepared
            or sha(base/"model.json") != frozen["model_sha256"]
            or (base/"FilmByPathColumn.jl").read_text() != kernel.source()):
        raise ValueError("Frozen model or predictions changed")
    for size, cases in frozen["cases"].items():
        for name, item in cases.items():
            setup = prepared["cases"][size][name]
            for filename in item["full"]+item["short"]:
                run = read_json(base/filename)
                conservation(run["result"], run["inputs"], setup)
            if "review" in item and item["review"] != review(base, item):
                raise ValueError("Column refinement differs on replay")
    if score(base, frozen, prepared) != record["figure7"]:
        raise ValueError("Figure 7 comparison differs on replay")
    print(dict(sources=len(record["source_sha256"]), artifacts=len(record["artifacts_sha256"])))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--holdout", type=Path, default=HOLDOUT)
    parser.add_argument("--check", type=Path)
    args = parser.parse_args()
    check(args.check.resolve()) if args.check else run(args.holdout.resolve())
