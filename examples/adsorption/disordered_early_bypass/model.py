"""Conservative tracer dose distributions and clean-film bounds."""
import numpy as np


def archive_lags(arrivals):
    delays = np.diff(arrivals,axis=1,prepend=0.)
    if not np.isfinite(delays).all() or not np.all(delays>0):
        raise ValueError("Incomplete or unordered first passages")
    return delays


def wall_stages(arrivals,hits,states):
    """Build operators from independent local delay/contact probes."""
    n = int(states.max())+1
    result = {key:[] for key in ("inert","receiver","pair","joint")}
    delays = archive_lags(arrivals)
    if hits.shape!=arrivals.shape or states.shape!=arrivals.shape or np.any(states<0):
        raise ValueError("Unexpected probe data shape")
    if np.isinf(hits).any() or np.any(hits>arrivals+1e-12) or np.any(hits<arrivals-delays-1e-12):
        raise ValueError("The first wall hit falls outside its section")
    clear = np.isnan(hits)
    donor = np.zeros(len(arrivals),dtype=int)
    for j,receiver in enumerate(states.T):
        pair = n*donor+receiver
        flow = 1/np.bincount(donor,minlength=n)[donor]
        counts = np.bincount(pair,minlength=n*n)
        successes = np.bincount(pair,weights=clear[:,j],minlength=n*n)
        receiver_count = np.bincount(receiver,minlength=n)
        receiver_clear = np.bincount(receiver,weights=clear[:,j],minlength=n)
        law = dict(inert=flow,receiver=flow*receiver_clear[receiver]/receiver_count[receiver],
                   pair=flow*successes[pair]/counts[pair],joint=flow*clear[:,j])
        for key in law:
            result[key].append((donor,receiver,delays[:,j],law[key]))
        donor = receiver
    return result


def nonnegative(values):
    if np.min(values)<-1e-12:
        raise ValueError("Nonphysical negative tracer mass")
    return np.maximum(values,0.)


def propagation(stages,step):
    """Advance pulse histograms one section at a time with linear convolutions."""
    n = max(max(i.max(),j.max()) for i,j,_,_ in stages)+1
    lag = [np.rint(t/step).astype(int) for _,_,t,_ in stages]
    length = 1 + sum(int(x.max()) for x in lag)
    nfft = 1 << (length-1).bit_length()
    spectrum = np.zeros((n,nfft//2+1),dtype=complex)
    spectrum[0] = 1.
    result = []
    state = np.eye(1,n).ravel()
    limit = 1
    for (i,j,_,w),bins in zip(stages,lag,strict=True):
        if not np.isfinite(w).all() or np.any(w<0):
            raise ValueError("Invalid conditional wall or carrier weight")
        linear = (i*n+j)*nfft+bins
        hist = np.bincount(linear,weights=w,minlength=n*n*nfft).reshape(n,n,nfft)
        kernel = np.fft.rfft(hist,axis=-1)
        spectrum = np.einsum("ik,ijk->jk",spectrum,kernel)
        limit += int(bins.max())
        probability = nonnegative(np.fft.irfft(spectrum.sum(axis=0),nfft)[:limit])
        amount = np.bincount(j,weights=state[i]*w,minlength=n)
        if amount.sum()>state.sum()+1e-12:
            raise ValueError("A wall operator created tracer")
        np.testing.assert_allclose(probability.sum(),amount.sum(),atol=1e-11,rtol=0.)
        result.append(probability)
        state = amount
    return result


def during_holdout(raw,step):
    arrival,hits = np.split(raw,2,axis=1)
    delays = archive_lags(arrival)
    if np.isinf(hits).any() or np.any(hits>arrival+1e-12) or np.any(hits<arrival-delays-1e-12):
        raise ValueError("Invalid held-out wall crossing times")
    no_hit = ~np.maximum.accumulate(np.isfinite(hits),axis=1)
    return [np.bincount(np.rint(t/step).astype(int),weights=clear/len(raw))
            for t,clear in zip(arrival.T,no_hit.T,strict=True)]


def compose(pattern,blocks,tail):
    whole,part = pattern[-1],pattern[tail-1] if tail else np.ones(1)
    length = blocks*(len(whole)-1)+len(part)
    nfft = 1 << (length-1).bit_length()
    mass = nonnegative(np.fft.irfft(np.fft.rfft(whole,nfft)**blocks*np.fft.rfft(part,nfft),nfft)[:length])
    np.testing.assert_allclose(mass.sum(),whole.sum()**blocks*part.sum(),atol=1e-11,rtol=0.)
    return mass


def cdf_at(probability,step,time):
    cdf = np.cumsum(probability)
    return np.interp(time,step*np.arange(len(cdf)),cdf,left=0.)


def maximum_cdf_error(a,b,step_a,step_b):
    time = np.union1d(step_a*np.arange(len(a)),step_b*np.arange(len(b)))
    return float(np.max(abs(cdf_at(a,step_a,time)-cdf_at(b,step_b,time))))


def plateau_bound(ntu,sections):
    """Minimum survival over contact allocations at fixed total film conductance."""
    if ntu<0 or sections<1:
        raise ValueError("Invalid total conductance or section count")
    return float(np.exp(-sections*np.log1p(ntu/sections)))


def required_ntu(target,sections):
    return float(sections*np.expm1(-np.log(target)/sections))
