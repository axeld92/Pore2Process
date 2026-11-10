"""Reuse the selected main-thread Julia runtime and frozen particle component."""
from functools import lru_cache
import numpy as np
from pore2process_demo.bridge import ROOT,runtime,_check_thread
from examples.adsorption.film_particle.bridge import module as particle_module


@lru_cache(maxsize=1)
def _loaded():
    particle_module()
    jl=runtime()
    jl.include(str(ROOT/'julia/adsorption/ColumnContact.jl'))
    return jl.ColumnContact


def module():
    _check_thread()
    return _loaded()


def contact(**inputs):
    _check_thread()
    jl=runtime();r=module().solve_contact(**inputs);result={}
    for field in jl.fieldnames(jl.typeof(r)):
        key=str(field);value=getattr(r,key)
        if value is None:result[key]=None
        elif bool(jl.isa(value,jl.AbstractArray)):
            result[key]=np.array(value,dtype=np.float64,copy=True).tolist()
        elif key in ('sections','axial_cells','radial_cells','steps','max_Newton_iterations'):
            result[key]=int(value)
        else:result[key]=float(value)
    return result
