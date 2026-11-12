"""Transfer residence-time rank dependence with fixed one-cell marginals."""
import numpy as np

from examples.adsorption.porous_outlet_tracer.run import PULSE_S
from examples.adsorption.residence_time_holdout.paths import moments, pulse_cdf
from examples.adsorption.residence_time_holdout.run import LIMITS, metrics


def ranks(values):
    order = np.argsort(values, kind="stable")
    rank = np.empty(len(values), dtype=int)
    rank[order] = np.arange(len(values))
    return rank


def equiprobable_delays(delays, weights, count):
    """Use the conditional mean delay in each equal-flow quantile bin."""
    probability = np.r_[0., np.cumsum(weights)]
    probability[-1] = 1.
    integral = np.r_[0., np.cumsum(delays * weights)]
    edges = np.linspace(0., 1., count + 1)
    j = np.minimum(np.searchsorted(probability[1:], edges, side="right"), len(delays) - 1)
    integrated = integral[j] + (edges - probability[j]) * delays[j]
    return count * np.diff(integrated)


def summary(times, delay, weight):
    cumulative = pulse_cdf(times, delay, weight, PULSE_S)
    mean, variance = moments(delay, weight, PULSE_S)
    return cumulative, dict(**metrics(times, cumulative), mean=mean, variance=variance)


def pair_prediction(data, delays, weights, times):
    finite = np.isfinite(data)
    completed = data[np.all(finite, axis=1)]
    if len(completed) / len(data) < .999 or np.any(completed <= 0.):
        raise ValueError("Marker recovery below 99.9% or nonpositive completed travel time")
    count = len(completed)
    r1, r2 = (ranks(x) for x in completed.T)
    equal = equiprobable_delays(delays, weights, count)
    np.testing.assert_allclose(equal.mean(), np.dot(delays, weights), rtol=1e-12, atol=0)
    paired = np.sort(equal[r1] + equal[r2])
    cumulative, result = summary(times, paired, np.full(count, 1 / count))
    d1, d2 = (np.minimum(r * 10 // count, 9) for r in (r1, r2))
    counts = np.zeros((10, 10), dtype=int)
    np.add.at(counts, (d1, d2), 1)
    conditional = counts / counts.sum(axis=1, keepdims=True)
    result.update(pair_recovery=count / len(data), first_crossing_recovery=float(finite[:, 0].mean()),
                  rank_correlation=float(np.corrcoef(r1, r2)[0, 1]),
                  calibrated_time_correlation=float(np.corrcoef(equal[r1], equal[r2])[0, 1]),
                  fast_retention=float(conditional[0, 0]), slow_retention=float(conditional[-1, -1]),
                  decile_counts=counts.tolist())
    return cumulative, result


def compare(actual, measured, prediction, estimated):
    error = {k: float(abs(estimated[k] / measured[k] - 1)) for k in LIMITS if k != "profile_rmse"}
    error["profile_rmse"] = float(np.sqrt(np.mean((actual - prediction)**2)))
    return dict(errors=error, passes={k: value <= LIMITS[k] for k, value in error.items()})
