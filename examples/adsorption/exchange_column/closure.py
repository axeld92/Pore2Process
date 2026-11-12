"""Coarsen the qualified marker rank pairs without reading Figure 7."""
import numpy as np

from examples.adsorption.path_exchange.predict import ranks, summary
from examples.adsorption.path_exchange.run import ARCHIVE, DT, PATHS_SHA
from pore2process_demo.io import read_json, sha
from pore2process_demo.paths import ROOT

MARKERS = ROOT / "results/adsorption/path_exchange/runs/20261112T073700Z-e27daca0"
EXCHANGE_SHA = "54dc188e143b6fe9431eae7f2e3bcb46f1ecf82eee5fa1bb82425113ed9690b5"
PRIMARY_SHA = "3c63cc7edb1dbf736350cedc657225a81a7389ce594d3892d7c3bdda5e37b746"
GRID_OPTIONS = (4, 8, 16, 32)
LIMITS = dict(max_cdf=.005, t99=.05, variance=.05)


def group(trajectories, delays, weights, edges):
    count = len(trajectories)
    boundaries = np.rint(count * edges).astype(int)
    shares = np.diff(boundaries) / count
    labels = [np.searchsorted(boundaries[1:], ranks(t), side="right") for t in trajectories.T]
    counts = np.zeros((len(shares), len(shares)), dtype=int)
    np.add.at(counts, tuple(labels), 1)
    np.testing.assert_array_equal(counts.sum(axis=0), np.diff(boundaries))
    np.testing.assert_array_equal(counts.sum(axis=1), np.diff(boundaries))
    probability = np.r_[0., np.cumsum(weights)]
    probability[-1] = 1.
    integral = np.r_[0., np.cumsum(delays * weights)]
    q = boundaries / count
    # The integral of the quantile function gives the exact bin mean.
    j = np.minimum(np.searchsorted(probability[1:], q, side="right"), len(delays)-1)
    conditional = np.diff(integral[j] + (q-probability[j])*delays[j]) / shares
    np.testing.assert_allclose(np.dot(shares, conditional), np.dot(weights, delays), rtol=1e-12)
    return shares, conditional, counts


def closure():
    marker_summary = MARKERS / "summary.json"
    if sha(marker_summary) != EXCHANGE_SHA or not read_json(marker_summary)["analysis"]["qualified"]:
        raise ValueError("Qualified path exchange reference changed")
    paths = ARCHIVE / "paths.csv"
    pairs = MARKERS / "primary.f64"
    if sha(paths) != PATHS_SHA or sha(pairs) != PRIMARY_SHA:
        raise ValueError("Frozen one-cell calibration or marker archive changed")
    parent = read_json(marker_summary)
    frozen = read_json(MARKERS / "frozen_prediction.json")
    archive_curve = MARKERS / "prediction.csv"
    if sha(MARKERS / "frozen_prediction.json") != parent["frozen_prediction_sha256"] or sha(archive_curve) != frozen["curve_sha256"]:
        raise ValueError("Frozen tracer prediction changed")
    times, _, base_cdf, *_ = np.loadtxt(archive_curve, delimiter=",", skiprows=1, unpack=True)
    reference = frozen["cases"]["primary"]
    delays, weights = np.loadtxt(paths, delimiter=",", skiprows=1, unpack=True)
    data = np.fromfile(pairs, dtype="<f8").reshape(-1, 2)
    data = data[np.isfinite(data).all(axis=1)]
    attempts = []
    for resolution in GRID_OPTIONS:
        edges = np.r_[np.linspace(0, .9, 3*resolution, endpoint=False),
                      np.linspace(.9, .99, resolution, endpoint=False),
                      np.linspace(.99, 1, resolution+1)]
        flow, delay, counts = group(data, delays, weights, edges)
        delay2 = (delay[:,None] + delay).ravel()
        sort = np.argsort(delay2)
        cdf, stats = summary(times, delay2[sort], (counts/data.shape[0]).ravel()[sort])
        errors = dict(max_cdf=float(max(abs(cdf-base_cdf))), t99=abs(stats["t99"]/reference["t99"]-1),
                      variance=abs(stats["variance"]/reference["variance"]-1))
        passed = all(errors[k] <= LIMITS[k] for k in LIMITS)
        attempts.append(dict(classes=len(flow), errors=errors, passed=passed))
        if passed:
            tau = float(np.dot(flow, delay))
            return dict(flow_fractions=flow.tolist(), normalized_transit_times=(delay/tau).tolist(),
                        transition_counts=counts.tolist(), reference_period_length_m=parent["mesh"]["length_m"]/2,
                        reference_liquid_transit_s=tau, attempts=attempts, limits=LIMITS,
                        primary_trajectory_sha256=PRIMARY_SHA, one_cell_delays_sha256=PATHS_SHA,
                        exchange_summary_sha256=EXCHANGE_SHA)
    raise ValueError(f"No predeclared path resolution passed the tracer coarsening check: {attempts}")
