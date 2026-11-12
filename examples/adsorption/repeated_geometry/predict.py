"""Propagate the fixed tracer exchange with and without within-class memory."""
import numpy as np

PERIODS = 8
FAST = .9
LIMITS = dict(survival=.01, conditional_capture=.01, fast_survival=.01)


def profiles(survival, flow, joint, groups, periods=PERIODS, reset_every=None):
    survival, flow, joint = map(np.asarray, (survival, flow, joint))
    transition = joint/joint.sum(axis=1, keepdims=True)
    np.testing.assert_allclose(flow@transition, flow, rtol=2e-11, atol=0.)
    edges = np.r_[0., np.cumsum(flow)]
    initially_fast = np.diff(np.minimum(edges, FAST))
    mass = np.stack((flow, initially_fast))
    start = mass.sum(axis=1)
    cumulative = []
    conditional = []
    for period in range(1, periods+1):
        before = mass.sum(axis=1)
        mass = mass*survival
        after = mass.sum(axis=1)
        cumulative.append(after/start)
        conditional.append((before-after)/before)
        mass = mass@transition
        if reset_every is not None and period % reset_every == 0:
            state = mass.reshape(2, -1, groups)
            shares = flow.reshape(-1, groups)
            shares = shares/shares.sum(axis=1, keepdims=True)
            mass = (state.sum(axis=2, keepdims=True)*shares).reshape(2, -1)
    return dict(survival=np.array(cumulative)[:, 0].tolist(),
                conditional_capture=np.array(conditional)[:, 0].tolist(),
                fast_survival=np.array(cumulative)[:, 1].tolist(),
                fast_conditional_capture=np.array(conditional)[:, 1].tolist())


def observe(data):
    arrivals, times = np.split(data, 2, axis=1)
    reached = np.isfinite(arrivals)
    if (np.isinf(data).any() or np.any(reached[:, 1:] & ~reached[:, :-1])
            or np.any(arrivals[:, 1:] < arrivals[:, :-1])
            or np.any(times > arrivals+1e-12)
            or np.any(times[:, 1:] < arrivals[:, :-1]-1e-12)):
        raise ValueError("Invalid first-passage sequence")
    contacted = np.logical_or.accumulate(np.isfinite(times), axis=1)
    uncontacted = ~contacted
    resolved = reached | contacted
    exits = reached & uncontacted
    fastest = np.argsort(arrivals[:, 0], kind="stable")[:round(FAST*len(data))]
    fast = exits[fastest]

    def conditional(survival):
        incoming = np.r_[1., survival[:-1]]
        return (1.-survival/incoming).tolist()

    total_survival = exits.mean(axis=0)
    fast_survival = fast.mean(axis=0)
    return dict(survival=total_survival.tolist(), conditional_capture=conditional(total_survival),
                fast_survival=fast_survival.tolist(),
                fast_conditional_capture=conditional(fast_survival),
                maximum_unresolved_fraction=float(max((~resolved).mean(axis=0))),
                recovery=reached.mean(axis=0).tolist(),
                survival_upper=(exits | ~resolved).mean(axis=0).tolist(),
                fast_survival_upper=(fast | ~resolved[fastest]).mean(axis=0).tolist(),
                median_inert_exit_s=float(np.median(arrivals[reached[:, -1], -1])))


def comparison(reference, prediction):
    errors = {key: float(max(abs(np.array(reference[key])[2:]-np.array(prediction[key])[2:])))
              for key in LIMITS}
    return dict(errors=errors, passed=all(value <= LIMITS[key] for key, value in errors.items()))
