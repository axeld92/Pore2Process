"""Audit the fixed Julia prediction using exact increments of inert variance."""
import argparse
from pathlib import Path

import numpy as np

from examples.adsorption.path_exchange.predict import ranks
from examples.adsorption.transport_accuracy.run import validate_record
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT

from .run import HERE, OUTPUT, archive_case, parent_inputs

HOLDOUT = OUTPUT/"runs/20261113T170400Z-e462ad77/summary.json"
HOLDOUT_SHA = "e5c81163a8224473cb61c96ac661529730cdddc32e0933463fc44e6695808484"
TERMS = ("exchange", "local_delay", "memory", "propagation")


def gain(joint, delay, delay2, donor_age):
    """Variance increment with elapsed age independent of the next delay given donor."""
    mean_delay = np.sum(joint*delay)
    mean_age = joint.sum(axis=1)@donor_age
    covariance = np.sum(joint*(donor_age[:,None]-mean_age)*delay)
    return float(np.sum(joint*delay2)-mean_delay**2+2*covariance)


def conditional_moments(counts, first_sum, second_sum, nearest=False):
    """Fill unsampled pairs within the same receiver class."""
    receiver_count = counts.sum(axis=0)
    if np.any(receiver_count==0):
        raise ValueError("A transit-time class has no markers")
    moments = [np.divide(x,counts,out=np.broadcast_to(x.sum(axis=0)/receiver_count,
                       counts.shape).copy(),where=counts>0) for x in (first_sum,second_sum)]
    if nearest:
        for i,j in zip(*np.nonzero(counts==0),strict=True):
            donors = np.flatnonzero(counts[:,j])
            index = donors[np.argmin(abs(donors-i))]
            for moment in moments:
                moment[i,j] = moment[index,j]
    return moments


def pair_table(labels, values=None):
    return np.bincount(labels,weights=values,minlength=100).reshape(10,10)


def seam(age, delay, donor, receiver, frozen_counts, frozen_delay, frozen_delay2, frozen_gain, *, nearest=False):
    labels = donor*10+receiver
    counts = pair_table(labels)
    donor_count = counts.sum(axis=1)
    if not np.array_equal(donor_count,frozen_counts.sum(axis=1)):
        raise ValueError("The donor flow fractions changed")
    observed = counts/len(age)
    frozen = frozen_counts/len(age)
    sums = pair_table(labels,delay)
    sums2 = pair_table(labels,delay**2)
    mean,second = conditional_moments(counts,sums,sums2,nearest)
    first,first2 = conditional_moments(frozen_counts,frozen_counts*frozen_delay,
                                      frozen_counts*frozen_delay2,nearest)
    donor_age = np.bincount(donor,weights=age,minlength=10)/donor_count
    f00 = gain(frozen,first,first2,donor_age)
    f10 = gain(observed,first,first2,donor_age)
    f01 = gain(frozen,mean,second,donor_age)
    f11 = gain(observed,mean,second,donor_age)
    exchange = .5*((f10-f00)+(f11-f01))
    local = .5*((f01-f00)+(f11-f10))
    joint_age = pair_table(labels,age)
    memory_sort = 2*np.sum((joint_age-counts*donor_age[:,None])*mean)/len(age)
    memory_pair = 2*np.sum(pair_table(labels,age*delay)-joint_age*mean)/len(age)
    actual_gain = float(np.var(age+delay)-np.var(age))
    result = dict(exchange=exchange,local_delay=local,memory=float(memory_sort+memory_pair),
                  propagation=f00-frozen_gain)
    np.testing.assert_allclose(sum(result.values()),actual_gain-frozen_gain,rtol=0.,atol=2e-15)
    return dict(**result,memory_sort=float(memory_sort),memory_within_pair=float(memory_pair),
                measured_gain=actual_gain,frozen_gain=float(frozen_gain),
                training_gap_probability=float(observed[frozen_counts==0].sum()),
                observed_gap_probability=float(frozen[counts==0].sum()))


def audit_case(arrivals, baseline, coefficients, *, nearest=False):
    duration = np.diff(arrivals,axis=1,prepend=0.)
    if not np.isfinite(duration).all() or np.any(duration<=0.):
        raise ValueError("The variance audit requires complete inert trajectories")
    limits = np.rint(len(arrivals)*np.arange(1,11)/10).astype(int)
    labels = np.column_stack([np.searchsorted(limits,ranks(x),side="right") for x in duration.T])
    counts,delay,delay2 = (np.asarray(coefficients[k]) for k in ("counts","delay","delay2"))
    np.testing.assert_array_equal(pair_table(labels[:,0]*10+labels[:,1]),counts)
    measured = arrivals.var(axis=0)
    predicted = np.asarray(baseline)
    rows = []
    cumulative = dict.fromkeys((*TERMS,"memory_sort","memory_within_pair"),0.)
    for section in range(2,8):
        result = seam(arrivals[:,section-1],duration[:,section],labels[:,section-1],labels[:,section],
                      counts,delay,delay2,predicted[section]-predicted[section-1],nearest=nearest)
        for key in cumulative:
            cumulative[key] += result[key]
        error = float(measured[section]-predicted[section])
        np.testing.assert_allclose(sum(cumulative[k] for k in TERMS),error,rtol=0.,atol=3e-15)
        rows.append(dict(section=section+1,increment_s2=result,cumulative_s2=cumulative.copy(),
                         measured_variance_s2=float(measured[section]),frozen_variance_s2=float(predicted[section]),
                         signed_error_s2=error))
    return rows


def evaluate():
    if sha(HOLDOUT)!=HOLDOUT_SHA:
        raise ValueError("The frozen holdout summary changed")
    holdout = validate_record(HOLDOUT)
    path = HOLDOUT.parent/"frozen_prediction.json"
    if sha(path)!=holdout["frozen_prediction_sha256"]:
        raise ValueError("The frozen prediction changed")
    frozen = read_json(path)
    summary,plan = parent_inputs()
    cases = {}
    sensitivity = {}
    for key,case in frozen["cases"].items():
        arrivals = archive_case(summary,plan,key)[:,:8]
        baseline = holdout["analysis"]["cases"][key]["pair"]["profiles"]["variance_s2"]
        cases[key] = audit_case(arrivals,baseline,case)
        alternate = audit_case(arrivals,baseline,case,nearest=True)
        sensitivity[key] = [dict(section=row["section"],cumulative_s2=row["cumulative_s2"])
                            for row in alternate]
    return dict(holdout_sha256=HOLDOUT_SHA,frozen_prediction_sha256=sha(path),
                parent_sha256=frozen["parent_sha256"],cases=cases,nearest_donor_sensitivity=sensitivity)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check",action="store_true")
    args = parser.parse_args()
    record = evaluate()
    target = HERE/"variance_reference.json"
    if args.check:
        if read_json(target)!=record:
            raise ValueError("The inert variance attribution changed on replay")
    else:
        write_json(target,record)
    print(target.relative_to(ROOT))
    for key,rows in record["cases"].items():
        for row in (rows[4],rows[-1]):
            scale = 100/row["measured_variance_s2"]
            print(key,row["section"],{k:round(scale*row["cumulative_s2"][k],2) for k in TERMS},
                  "net",round(scale*row["signed_error_s2"],2))


if __name__=="__main__":
    main()
