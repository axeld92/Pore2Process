"""Use the existing single Julia runtime and project; never mutate old modules."""
from functools import lru_cache

import numpy as np

from pore2process_demo.bridge import ROOT, runtime, _check_thread


@lru_cache(maxsize=1)
def _loaded():
    jl=runtime()
    jl.include(str(ROOT/'julia/adsorption/FilmParticle.jl'))
    return jl.FilmParticle


def module():
    _check_thread()
    return _loaded()


def call(function, **inputs):
    _check_thread()
    jl=runtime()
    result=getattr(module(),function)(**inputs)
    out={}
    for field in jl.fieldnames(jl.typeof(result)):
        key=str(field);value=getattr(result,key)
        if value is None:
            out[key]=None
        elif bool(jl.isa(value,jl.AbstractArray)):
            out[key]=np.array(value,dtype=np.float64,copy=True).tolist()
        elif key in ('axial_cells','radial_cells','steps','max_Newton_iterations'):
            out[key]=int(value)
        else:
            out[key]=float(value)
    return out


def batch(**inputs):
    return call('solve_bath',**inputs)


def bed(**inputs):
    return call('solve_column',**inputs)
