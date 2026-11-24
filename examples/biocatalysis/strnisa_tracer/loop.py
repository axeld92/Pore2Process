"""Measured blank-loop convolution, preserving absolute dose and shape."""
import math
from pathlib import Path
import tomllib

import numpy as np

from pore2process_demo.paths import ROOT

MARKERS=ROOT/"literature/sources/strnisa2018/rtd_fig3b.toml"


def digitized(name):
    with MARKERS.open("rb") as f:
        record=tomllib.load(f)
    i={"blank":1,"A":2,"C":3}[name]
    a=np.array(record["markers"],dtype=float)
    t=np.r_[0.,a[:,0]*record["marker_spacing_min"]*60]
    pixels=np.r_[record["axis_conc_0_mM_px"],a[:,i]]
    unit=4/(record["axis_conc_0_mM_px"]-record["axis_conc_4_mM_px"])
    c=(record["axis_conc_0_mM_px"]-pixels)*unit
    return t,c,unit


def measured(name,dt_s=.05,background_pixels=0):
    """Common reading rule; do not estimate any correction from the holdout."""
    time,c,pixel= digitized(name)
    # The source baseline is the zero-mM tick. Positive tail readings are
    # retained; the alternative 2/5 px backgrounds are separate scenarios.
    c=np.maximum(c-background_pixels*pixel,0.)
    t=np.arange(round(time[-1]/dt_s)+1,dtype=float)*dt_s
    return t,np.interp(t,time,c)


def moments(t, signal, flow_ul_min=100.):
    t=np.asarray(t,dtype=float)
    signal=np.asarray(signal,dtype=float)
    if t.shape!=signal.shape or np.any(np.diff(t)<=0) or not np.isfinite(signal).all() or np.min(signal)<0:
        raise ValueError("Invalid nonnegative pulse data")
    area=float(np.trapezoid(signal,t))
    if area<=0:
        raise ValueError("A pulse must have positive area")
    mean=float(np.trapezoid(t*signal,t)/area)
    var=float(np.trapezoid((t-mean)**2*signal,t)/area)
    area_min=area/60
    return dict(area_mM_min=area_min,observed_glucose_umol=area_min*flow_ul_min/1000,
                mean_s=mean,variance_s2=var,
                mean_min=mean/60,variance_min2=var/3600,
                skewness=float(np.trapezoid((t-mean)**3*signal,t)/area/var**1.5) if var>0 else 0.,
                peak_mM=float(signal.max()),peak_s=float(t[np.argmax(signal)]))


def tracer_kernel(times,dt_s):
    """Assign each resolved marker to two adjacent time bins conserving its mass and first moment."""
    times=np.asarray(times,dtype=float)
    valid=times[np.isfinite(times)]
    if len(times)==0 or len(valid)==0 or min(valid)<0 or dt_s<=0:
        raise ValueError("Invalid resolved tracer ages")
    z=valid/dt_s
    lo=np.floor(z).astype(int)
    fraction=z-lo
    weight=1/len(times)
    size=max(lo)+2
    k=np.bincount(lo,weights=(1-fraction)*weight,minlength=size)
    k+=np.bincount(lo+1,weights=fraction*weight,minlength=size)
    np.testing.assert_allclose(k.sum(),len(valid)/len(times),rtol=1e-13)
    np.testing.assert_allclose(k @ (np.arange(size)*dt_s),valid.sum()/len(times),rtol=1e-13,atol=1e-14)
    return k


def convolve(blank, kernel, dt_s):
    """Blank already includes the injected pulse, tubing and detector."""
    b=np.asarray(blank,dtype=float)
    k=np.asarray(kernel,dtype=float)
    if min(k)<0 or not 0<k.sum()<=1+1e-12:
        raise ValueError("Kernel must be a physical exit-mass distribution")
    result=np.convolve(b,k)
    return np.arange(len(result))*dt_s,result


def cumulative(t,c):
    increments=(c[:-1]+c[1:])*np.diff(t)/2
    integral=np.r_[0,np.cumsum(increments)]
    return integral/integral[-1]


def compare_prediction(t_pred,pred,name,dt_s=.05,background_pixels=0):
    t_obs,obs=measured(name,dt_s,background_pixels)
    t_blank,blank=measured("blank",dt_s,background_pixels)
    observed=moments(t_obs,obs)
    predicted=moments(t_pred,pred)
    b=moments(t_blank,blank)
    t=np.arange(max(len(pred),len(obs)))*dt_s
    cdfp=np.interp(t,t_pred,cumulative(t_pred,pred),left=0.,right=1.)
    cdfo=np.interp(t,t_obs,cumulative(t_obs,obs),left=0.,right=1.)
    def relative(x,y):
        return x/y-1 if y>0 else None
    delta_actual={k:observed[k]-b[k] for k in ("mean_s","variance_s2")}
    delta_model={k:predicted[k]-b[k] for k in ("mean_s","variance_s2")}
    return dict(observed=observed,predicted=predicted,blank=b,
                whole_mean_relative=relative(predicted["mean_s"],observed["mean_s"]),
                whole_variance_relative=relative(predicted["variance_s2"],observed["variance_s2"]),
                whole_CDF_max=float(np.max(np.abs(cdfp-cdfo))),
                observed_CDF_ahead_of_prediction_max=float(np.max(cdfo-cdfp)),
                prediction_CDF_ahead_of_observation_max=float(np.max(cdfp-cdfo)),
                predicted_mass_after_measured_window_fraction=float(1-np.interp(t_obs[-1],t_pred,cumulative(t_pred,pred),right=1.)),
                incremental_observed=delta_actual,incremental_predicted=delta_model,
                incremental_mean_relative=relative(delta_model["mean_s"],delta_actual["mean_s"]),
                incremental_variance_relative=relative(delta_model["variance_s2"],delta_actual["variance_s2"]),
                predicted_to_observed_dose_ratio=relative(predicted["observed_glucose_umol"],
                                                      observed["observed_glucose_umol"])+1)


def compare_times(times,name,dt_s=.05,background_pixels=0,extra_kernel=None):
    _,blank=measured("blank",dt_s,background_pixels)
    k=tracer_kernel(times,dt_s)
    if extra_kernel is not None:
        k=np.convolve(k,extra_kernel)
    t,pred=convolve(blank,k,dt_s)
    return t,pred,compare_prediction(t,pred,name,dt_s,background_pixels)


def exponential_connection(mean_s,dt_s):
    """One additional fully mixed volume, on top of the measured blank.

    For a zero volume return an impulse; otherwise bin exact exponential
    probabilities. This is only an A-calibrated sensitivity when needed,
    not an independently measured connector.
    """
    if mean_s<0 or dt_s<=0:
        raise ValueError("Nonphysical extra volume")
    if mean_s==0:
        return np.ones(1)
    n=max(1,math.ceil(25*mean_s/dt_s))
    # Exact integrated mass centered at bin k; truncate and keep mean
    # accurate compared with the much coarser experimental markers.
    edges=np.maximum(0.,(np.arange(n+1)-.5)*dt_s)
    k=np.exp(-edges[:-1]/mean_s)-np.exp(-edges[1:]/mean_s)
    k[-1]+=np.exp(-edges[-1]/mean_s)
    return k


def rare_delay(rate,mean_event_s,dt_s):
    """Conditional diagnostic: Poisson number of independent exponential pauses.

    Mean extra delay = rate * mean_event_s; variance = 2 * rate * mean_event_s**2.
    A single shared process and a process proportional to solid volume can use
    the very same parameters on A. Neither identifies the physical location of
    a pause. Bin the point mass and gamma mixtures, without dropping late mass.
    """
    if rate<0 or mean_event_s<=0 or dt_s<=0 or not math.isfinite(rate+mean_event_s+dt_s):
        raise ValueError("Invalid rare-delay parameters")
    if rate==0:
        return np.ones(1)
    # This diagnostic is limited to rare delays, not an arbitrary fitting family.
    if rate>10:
        raise ValueError("Not a rare-delay model")
    weights=[math.exp(-rate)]
    while 1-math.fsum(weights)>1e-13:
        weights.append(weights[-1]*rate/len(weights))
    size=math.ceil((len(weights)+12*math.sqrt(len(weights))+25)*mean_event_s/dt_s)+1
    edges=np.maximum(0.,(np.arange(size+1)-.5)*dt_s)/mean_event_s
    exponential=np.exp(-edges)
    term=exponential.copy()
    survival=term.copy()
    k=np.zeros(size)
    k[0]=weights[0]
    for count,weight in enumerate(weights[1:],1):
        probabilities=np.maximum(survival[:-1]-survival[1:],0.)
        probabilities[-1]+=survival[-1]
        k+=weight*probabilities
        term*=edges/count
        survival+=term
    # The omitted Poisson-count probability is less than 1e-13.
    k/=k.sum()
    return k
