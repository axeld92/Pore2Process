"""Finer, separately guarded reference after the sampled N64 refinement screen."""
import numpy as np
from .reference import concentration,mode,boundary_checks


def converged_reference(x_m,y_m,*,amplitude_kg_m3=.25,**physical):
    coarse=concentration(x_m,y_m,amplitude_kg_m3=amplitude_kg_m3,points=2048,steps=2048,**physical)
    spatial=concentration(x_m,y_m,amplitude_kg_m3=amplitude_kg_m3,points=4096,steps=2048,**physical)
    fine=concentration(x_m,y_m,amplitude_kg_m3=amplitude_kg_m3,points=4096,steps=4096,**physical)
    scale=max(amplitude_kg_m3,1e-30)
    checks=boundary_checks(mode(points=4096,steps=4096,**physical),physical['height_m'])
    result=dict(spatial_max_over_amplitude=float(np.max(abs(coarse-spatial)))/scale,
                temporal_max_over_amplitude=float(np.max(abs(spatial-fine)))/scale,
                points_on_physical_height=4096,points_on_even_period=8192,steps=4096,tolerance=1e-4,**checks)
    result['passed']=max(result['spatial_max_over_amplitude'],result['temporal_max_over_amplitude'])<=1e-4 and checks['even_reflection_error']<=1e-11 and checks['max_wall_normal_mode_derivative_times_height']<=1e-9
    if not result['passed']:raise ValueError('Fine Neumann reference not qualified: '+str(result))
    return fine,result
