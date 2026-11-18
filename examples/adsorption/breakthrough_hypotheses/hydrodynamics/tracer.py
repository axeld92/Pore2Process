"""Check one confined period with reflecting-tube and clean-carbon tracer outcomes."""
import argparse
from pathlib import Path

import numpy as np

from examples.adsorption.disordered_packing import mesh
from examples.adsorption.dispersion_holdout.run import execute
from pore2process_demo.io import read_json, sha, write_json

from . import native
from . import packing as p
from .handoff import julia

CUT_DP = .5
STEPS = (.001, .0005)
SEED = 8401
MARKERS = 16384
MAX_REFLUX = .02


def select_cut(records, fraction=None):
    """Minimize worst measured reflux, before any tracer/adsorption calculation."""
    width = max(len(r["cut_backflow_fraction"]) for r in records)
    candidates = np.arange(width) / width if fraction is None else np.array([fraction])
    reflux = []
    for r in records:
        fraction_array = np.array(r["cut_backflow_fraction"])
        cuts = np.rint(candidates * len(fraction_array)).astype(int)
        reflux.append(fraction_array[(cuts - 1) % len(fraction_array)])
    cost = np.max(reflux, axis=0)
    return float(candidates[np.argmin(cost)])


def compact(base, folder, fraction):
    ijk, faces, q, record = mesh.load_faces(folder)
    info = record["info"]
    shape = tuple(info["shape"])
    nx, width, _ = shape
    h = info["raw"]["dx_m"]
    cut = round(fraction * nx) % nx
    a, b, axis = (faces[name]-1 for name in ("left.i32", "right.i32", "axis.u8"))
    external = (axis == 0) & (ijk[b, 0] == cut)
    before = np.zeros((len(ijk), 3), dtype="<i4")
    after = np.zeros_like(before)
    low = np.zeros(before.shape, dtype="<f8")
    high = np.zeros_like(low)
    low[b, axis] = q
    high[a, axis] = q
    internal = ~external
    before[b[internal], axis[internal]] = a[internal] + 1
    after[a[internal], axis[internal]] = b[internal] + 1
    _, spheres = native.geometry(base)
    raster = p.voxel_material(spheres, shape, h)
    types = np.zeros((len(ijk), 6), dtype="u1")
    wall = np.zeros_like(types)
    for direction in (-1, 1):
        for dimension in range(3):
            column = 2 * dimension + (direction > 0)
            adjacent = ijk.copy()
            adjacent[:, dimension] = (adjacent[:, dimension] + direction) % shape[dimension]
            types[:, column] = raster[tuple(adjacent.T)]
            transverse = (ijk[:, 1:] + .5 - width / 2) * h
            if dimension:
                transverse[:, dimension-1] += direction * h / 2
            wall[:, column] = np.sum(transverse**2, axis=1) >= (p.RADIUS - CUT_DP * p.DIAMETER)**2
    types[a[external], 1] = 4
    types[b[external], 0] = 4
    np.testing.assert_array_equal(types[:, 0::2] == 1, before > 0)
    np.testing.assert_array_equal(types[:, 1::2] == 1, after > 0)
    qf = float(np.sum(np.maximum(0, q[external])))
    qr = float(np.sum(np.maximum(0, -q[external])))
    np.testing.assert_allclose(qf - qr, info["directional_flow_m3_s"][0], rtol=1e-9, atol=0.)
    imbalance = np.sum(high - low, axis=1)
    # A uniform unit concentration has exactly the same discrete divergence as the carrier.
    uniform_error = float(np.max(np.abs(imbalance)) / np.max(np.abs(q)))
    if uniform_error > 1e-10:
        raise ValueError("Uniform-concentration shared flux is not balanced")
    path = folder / "tracer"
    path.mkdir(exist_ok=True)
    data = {"before.i32": before, "after.i32": after, "incoming.f64": low, "outgoing.f64": high,
            "coords.i32": ijk.astype("<i4", copy=False), "material.u8": types, "carbon_face_in_shell.u8": wall}
    for filename, value in data.items():
        value.tofile(path / filename)
    carbon = types == 2
    summary = dict(cut_index=cut, cut_fraction=float(cut / nx), selected_cut_fraction=fraction,
        section_length_m=nx * h, cut_dp=CUT_DP, throughflow_m3_s=qf-qr,
        inlet_forward_flow_m3_s=qf, cut_reflux_m3_s=qr, cut_reflux_fraction=qr / qf,
        discretized_wall_carbon_face_fraction=float(np.sum(carbon * wall) / np.sum(carbon)),
        connected_liquid_volume_m3=info["liquid_volume_m3"], uniform_tracer_flux_imbalance_relative=uniform_error,
        dx_m=h, transverse_cells=width, files={name: sha(path / name) for name in data},
        face_summary_sha256=sha(folder / "faces/summary.json"), radial_summary_sha256=sha(folder / "radial.json"))
    write_json(path / "compact.json", summary)
    return path, summary


def statistics(data, holding):
    start, end, direction, endtime, first_carbon, carbon_zone, first_tube, local_wall, local_core, wall_s, surviving_wall_s = data.T
    if not (np.all(np.isin(start, (1, 2))) and np.all(np.isin(end, (1, 2))) and
            np.all(np.isin(direction, (-1, 0, 1))) and np.all(np.isin(carbon_zone, (0, 1, 2)))):
        raise ValueError("Unknown tracer region or fate")
    if not (np.all((wall_s >= 0.) & (wall_s <= endtime + 1e-10)) and np.all((surviving_wall_s >= 0.) &
              (surviving_wall_s <= np.where(carbon_zone, first_carbon, endtime) + 1e-10))):
        raise ValueError("Invalid wall residence times")
    outlet = direction == 1
    returns = direction == -1
    incomplete = direction == 0
    assert int(outlet.sum()+returns.sum()+incomplete.sum()) == len(data)
    tau = endtime[outlet]
    if not len(tau):
        raise ValueError("No downstream tracer arrivals")
    clean_exit = outlet & (carbon_zone == 0)
    unresolved = incomplete & (carbon_zone == 0)
    hit = carbon_zone > 0
    assert np.array_equal(hit, np.isfinite(first_carbon))
    assert np.all(first_carbon[hit] <= endtime[hit])
    survival_s = np.where(hit, first_carbon, endtime)
    hit_by_zone = np.array([np.sum(carbon_zone == zone) for zone in (1, 2)])
    exposure_s = np.sum([surviving_wall_s, survival_s - surviving_wall_s], axis=1)
    return dict(marker_count=len(data), region_order=["wall", "core"],
        inert=dict(outlet_recovered_fraction=float(outlet.mean()), inlet_return_fraction=float(returns.mean()),
            mass_recovered_fraction=float(1 - incomplete.mean()), censored_fraction=float(incomplete.mean()),
            downstream_mean_residence_s=float(tau.mean()), downstream_residence_variance_s2=float(tau.var()),
            expected_holdup_over_net_flow_s=holding, downstream_relative_mean_error=float(abs(tau.mean() - holding) / holding),
            transit_quantiles_s=dict(zip(("t05", "t10", "t50", "t90", "t95"), np.quantile(tau, (.05,.1,.5,.9,.95)).tolist(), strict=True)),
            mean_wall_residence_fraction=float(np.mean(wall_s[outlet] / tau)),
            downstream_mean_residence_by_inlet_zone_s=[float(endtime[outlet & (start == zone)].mean()) for zone in (1, 2)],
            downstream_mean_residence_by_outlet_zone_s=[float(endtime[outlet & (end == zone)].mean()) for zone in (1, 2)],
            conditional_wall_residence_fraction=[float(np.mean(wall_s[outlet & (start == zone)] / endtime[outlet & (start == zone)]))
                                                for zone in (1, 2)]),
        transitions=dict(start_fraction=[float(np.mean(start == zone)) for zone in (1, 2)],
            outlet_fraction_given_outlet=[float(np.mean(end[outlet] == zone)) for zone in (1, 2)],
            wall_core_conditional_on_outlet=[[float(np.mean(end[outlet & (start == source)] == dest))
                                            for dest in (1, 2)] for source in (1, 2)],
            outlet_from_source=[[float(np.mean(outlet[start == source] & (end[start == source] == dest)))
                                 for dest in (1, 2)] for source in (1, 2)],
            return_from_source=[float(returns[start == zone].mean()) for zone in (1, 2)],
            unresolved_from_source=[float(incomplete[start == zone].mean()) for zone in (1, 2)]),
        clean_carbon=dict(outlet_escape_fraction=float(clean_exit.mean()),
            escape_given_inlet_zone=[float(clean_exit[start == zone].mean()) for zone in (1, 2)],
            first_hit_count_by_wall_core=hit_by_zone.tolist(),
            first_hit_probability_by_wall_core=(hit_by_zone / len(data)).tolist(),
            first_tube_touch_fraction=float(np.isfinite(first_tube).mean()),
            inlet_return_fraction=float(np.mean(returns & ~hit)), unresolved_fraction=float(unresolved.mean()),
            surviving_marker_seconds_by_zone=exposure_s.tolist(),
            effective_first_contact_rate_per_s=(hit_by_zone / exposure_s).tolist(),
            mean_wall_and_core_reflection_length_m=[float(local_wall.mean()), float(local_core.mean())]))


def particle_checks(base, cases, fraction, threads, markers):
    records = [read_json(base / case / "radial.json") for case in cases]
    selected = select_cut(records, fraction)
    for case in cases:
        folder = base / case
        path, comp = compact(base, folder, selected)
        holding = comp["connected_liquid_volume_m3"] / comp["throughflow_m3_s"]
        refinements = {}
        for dt in STEPS:
            out = path / f"markers_dt_{dt:.5f}.f64"
            # Use the same physical stop with both step sizes.
            stop = np.ceil(16 * holding / STEPS[0]) * STEPS[0]
            cmd = [*julia(threads), str(native.HERE / "track.jl"), str(path), str(out),
                   *map(str, (comp["dx_m"], p.DIFFUSIVITY, dt, stop, markers, SEED, comp["transverse_cells"],
                              p.RADIUS - CUT_DP * p.DIAMETER))]
            print(f"Tracing {markers} markers for {case} at dt={dt:g} s", flush=True)
            elapsed = execute(cmd, out.with_suffix(".log"), timeout=1800)
            sample = np.fromfile(out, dtype="<f8").reshape(markers, 11)
            entry = statistics(sample, holding)
            entry.update(dt_s=dt, stop_s=stop, marker_sha256=sha(out), command=cmd, elapsed=elapsed)
            refinements[str(dt)] = entry
        coarse, fine = [refinements[str(dt)] for dt in STEPS]
        exposure_change = abs(coarse["clean_carbon"]["outlet_escape_fraction"] - fine["clean_carbon"]["outlet_escape_fraction"])
        transition_change = float(np.max(abs(np.array(coarse["transitions"]["wall_core_conditional_on_outlet"]) -
                                            np.array(fine["transitions"]["wall_core_conditional_on_outlet"]))))
        mean_residence_change = abs(coarse["inert"]["downstream_mean_residence_s"] - fine["inert"]["downstream_mean_residence_s"]) / holding
        radial = read_json(folder / "radial.json")
        gates = dict(reflux=comp["cut_reflux_fraction"] <= MAX_REFLUX,
            inert_recovery=fine["inert"]["mass_recovered_fraction"] >= .99,
            mean_residence=fine["inert"]["downstream_relative_mean_error"] <= .05,
            time_step_escape=exposure_change <= .02, time_step_transition=transition_change <= .02,
            time_step_mean_residence=mean_residence_change <= .05,
            local_face_correction=radial["face_checks"]["qualifies_local_correction"])
        path_result = dict(compact=comp, diffusion_m2_s=p.DIFFUSIVITY, seed=SEED, steps=refinements,
            gates=gates, qualified=all(gates.values()), dt_escape_fraction_difference=exposure_change,
            dt_transition_probability_max_difference=transition_change, dt_mean_residence_relative_change=mean_residence_change,
            limitations="One short synthetic period. Forward-flux inlet pulse with counted upstream returns. "
                        "Molecular diffusion reflected at artificial inlet and outlet cuts. First carbon hit is an ideal sink; "
                        "neither adsorption nor empirical curve fitting is used.")
        write_json(path / "summary.json", path_result)
        print(f"{case}: outcomes {fine['inert']}, {fine['transitions']}, {fine['clean_carbon']}; gates {gates}", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("base", type=Path)
    parser.add_argument("cases", nargs="+")
    parser.add_argument("--cut-fraction", type=float)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--markers", type=int, default=MARKERS)
    args = parser.parse_args()
    if (args.cut_fraction is not None and not 0 <= args.cut_fraction < 1) or args.threads < 1 or args.markers < 1:
        parser.error("Expected a cut in [0,1), positive threads and marker count")
    particle_checks(args.base.resolve(), args.cases, args.cut_fraction, args.threads, args.markers)


if __name__ == "__main__":
    main()
