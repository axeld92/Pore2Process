"""One explicit Julia runtime, called on the main thread."""
from functools import lru_cache
import numpy as np
from pore2process_demo.bridge import ROOT,runtime,_check_thread

@lru_cache(maxsize=1)
def _loaded():
    jl=runtime();jl.include(str(ROOT/'julia/adsorption/TracerContact.jl'));return jl.TracerContact

def call(function,**inputs):
    _check_thread();jl=runtime();r=getattr(_loaded(),function)(**inputs);out={}
    for field in jl.fieldnames(jl.typeof(r)):
        key=str(field);value=getattr(r,key)
        out[key]=np.array(value,dtype=float,copy=True).tolist() if bool(jl.isa(value,jl.AbstractArray)) else float(value)
    return out

def moments(time_s,flux_kg_s):
    _check_thread();jl=runtime();r=_loaded().step_moments(jl.Vector[jl.Float64](time_s),jl.Vector[jl.Float64](flux_kg_s))
    return {str(k):float(getattr(r,str(k))) for k in jl.fieldnames(jl.typeof(r))}
