"""Reduce paired inert trajectories to an exchange with wall-exposure memory."""
import numpy as np

from examples.adsorption.contact_column.inputs import rank_bins
from examples.adsorption.path_exchange.predict import ranks

QUARTILES = 4
LIMITS = dict(two_period_survival_error=.01,
              fresh_advection_class_capture_error=.01,
              survivor_conditional_capture_error=.01)


def states(transit, twins, flow, subgroups):
    """Preserve the archived parent ranks, then resolve the local twin ranks."""
    parent = rank_bins(transit, np.r_[0., np.cumsum(flow)])
    if subgroups == 1:
        return parent
    result = np.empty_like(parent)
    for period in range(2):
        for i in range(len(flow)):
            members = np.flatnonzero(parent[:, period] == i)
            local = ranks(twins[members, period])
            if len(members) < subgroups:
                raise ValueError("Too few trajectories to resolve wall-exposure subgroups")
            result[members, period] = i*subgroups + local*subgroups//len(members)
    return result


def carrier_flux(labels, flow, receiving, subgroups):
    """Return joint inlet-to-outlet flux, with the original parent-pair fluxes."""
    n = len(flow)
    target = receiving.T * flow[None, :]
    for marginal in (target.sum(axis=0), target.sum(axis=1)):
        np.testing.assert_allclose(marginal, flow, rtol=0., atol=5e-15)
    if subgroups == 1:
        return target
    fine_flow = np.repeat(flow/subgroups, subgroups)
    size = len(fine_flow)
    joint = np.zeros((size, size))
    np.add.at(joint, tuple(labels.T), 1)
    block = joint.reshape(n, subgroups, n, subgroups)
    block += (target[:, None, :, None] > 0)/(subgroups**2)
    block *= target[:, None, :, None] > 0
    for _ in range(2000):
        joint *= (fine_flow/joint.sum(axis=1))[:, None]
        joint *= (fine_flow/joint.sum(axis=0))[None, :]
        total = block.sum(axis=(1, 3))
        block *= np.divide(target, total, out=np.zeros_like(target), where=total > 0)[:, None, :, None]
        error = max(np.max(abs(joint.sum(axis=axis)-fine_flow)/fine_flow) for axis in (0, 1))
        if error <= 1e-11:
            np.testing.assert_allclose(block.sum(axis=(1, 3)), target, rtol=2e-14, atol=0.)
            return joint
    raise ValueError("Could not balance subgroup fluxes to the fixed Péclet exchange")


def predict(labels, first_survives, flow, receiving, subgroups):
    joint = carrier_flux(labels, flow, receiving, subgroups)
    fine_flow = np.repeat(flow/subgroups, subgroups)
    counts = np.bincount(labels[:, 0], minlength=len(fine_flow))
    survival = np.bincount(labels[:, 0], weights=first_survives, minlength=len(fine_flow))/counts
    incoming = survival @ joint
    first_survival = float(fine_flow @ survival)
    survived_twice = float(incoming @ survival)
    np.testing.assert_allclose(incoming.sum(), first_survival, rtol=1.1e-11, atol=0.)
    return dict(subgroups=subgroups, carrier_flux=joint.tolist(), flow=fine_flow.tolist(),
                fresh_survival=survival.tolist(),
                seam_survival_fraction=(incoming/fine_flow).tolist(),
                first_survival=first_survival, two_period_survival=survived_twice,
                survivor_conditional_capture=1.-survived_twice/first_survival)


def class_average(values, labels, count):
    return np.bincount(labels, weights=values, minlength=count)/np.bincount(labels, minlength=count)


def score(frozen, labels, advection, first_survives, second_hit, flow):
    groups = frozen["subgroups"]
    parent = labels[:, 1]//groups
    fresh = 1.-np.array(frozen["fresh_survival"])[labels[:, 1]]
    local_if_observed_seam = float(fresh[first_survives].mean())
    observed_survivor_capture = float(second_hit[first_survives].mean())
    two_period_survival = float(np.mean(first_survives & ~second_hit))
    count = len(flow)
    errors = {}
    classes = {}
    for name, local_labels in (("advection", advection), ("diffusive", parent)):
        measured = class_average(second_hit, local_labels, count)
        predicted = class_average(fresh, local_labels, count)
        errors[f"fresh_{name}_class_capture_error"] = float(flow @ abs(predicted-measured))
        classes[name] = dict(measured_capture=measured.tolist(), predicted_capture=predicted.tolist())
    errors.update(two_period_survival_error=abs(frozen["two_period_survival"]-two_period_survival),
                  survivor_conditional_capture_error=abs(frozen["survivor_conditional_capture"]
                                                        -observed_survivor_capture))
    return dict(errors=errors, passed=all(errors[k] <= v for k, v in LIMITS.items()),
                two_period_survival=two_period_survival,
                observed_survivor_conditional_capture=observed_survivor_capture,
                predicted_survivor_conditional_capture=frozen["survivor_conditional_capture"],
                observed_seam_predicted_conditional_capture=local_if_observed_seam,
                observed_reset_total_capture=float(second_hit.mean()),
                predicted_reset_total_capture=float(fresh.mean()), classes=classes)
