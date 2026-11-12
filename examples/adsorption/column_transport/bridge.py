"""Main-thread reference calculations and an explicit native qualification gate."""
from functools import lru_cache
import math
import numpy as np

from pore2process_demo.bridge import runtime, _check_thread
from pore2process_demo.paths import ROOT
from examples.adsorption.film_particle.bridge import bed


@lru_cache(maxsize=1)
def _reference():
    jl=runtime()
    jl.include(str(ROOT/'julia/adsorption/ColumnBoundaryReference.jl'))
    return jl.ColumnBoundaryReference


def boundary_moments(kind, **inputs):
    _check_thread()
    if kind not in ('dirichlet_neumann','danckwerts'):
        raise ValueError('Select dirichlet_neumann or danckwerts inlet/outlet reference.')
    jl=runtime()
    r=getattr(_reference(), kind+'_moments')(**inputs)
    return {str(k):float(getattr(r,str(k))) for k in jl.fieldnames(jl.typeof(r))}


def reference_profile(**inputs):
    _check_thread()
    jl=runtime()
    r=_reference().solve_dirichlet_neumann(**inputs)
    return {str(k):np.array(getattr(r,str(k)),dtype=float,copy=True).tolist()
            if bool(jl.isa(getattr(r,str(k)),jl.AbstractArray)) else float(getattr(r,str(k)))
            for k in jl.fieldnames(jl.typeof(r))}


REQUIRED_SCREENS = ('uniform_preservation','mass','bounds','mesh','time',
                    'convective_outlet','longer_segment_prediction')


def qualified_dispersion(closure):
    """A saved native closure must qualify before reaching particle kinetics."""
    if closure.get('boundary') != 'convective-column':
        raise ValueError('Native closure must have qualified convective-column boundaries.')
    checks=closure.get('checks', {})
    failed=[k for k in REQUIRED_SCREENS if checks.get(k) is not True]
    if failed:
        raise ValueError('Native transport is unqualified: '+', '.join(failed))
    if closure.get('calibration_boundary') != 'native-inert-tracer-only':
        raise ValueError('Column curves cannot independently calibrate the transport closure.')
    D=closure.get('axial_dispersion_m2_s')
    if isinstance(D,bool) or not isinstance(D,(int,float)) or not math.isfinite(D) or D<=0:
        raise ValueError('Qualified axial dispersion must be a positive finite SI value.')
    if not closure.get('native_run_id'):
        raise ValueError('Native source run must remain traceable.')
    return float(D)


def adsorption_column(closure, **inputs):
    """Reuse the frozen Julia particle/column component; no alternate solver."""
    _check_thread()
    D=qualified_dispersion(closure)
    if inputs.get('active_inventory_fraction',1.) != 1.:
        raise ValueError('This handoff preserves the full declared dry adsorbent inventory.')
    if 'axial_dispersion_m2_s' in inputs:
        raise ValueError('Supply dispersion only through the qualified native closure.')
    return bed(**inputs, axial_dispersion_m2_s=D)
