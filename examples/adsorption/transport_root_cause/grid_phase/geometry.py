"""Analytical FCC surface distances for read-only diagnostics, not meshing."""
import numpy as np


def sphere_distance(xyz_m, *, side_m, diameter_m=.0008, phase_x_m=0.):
    """Nearest periodic sphere's signed surface gap and outward normal.

    OpenLB stock indicators alone construct the simulation geometry. This
    independent distance calculation audits exported point locations.
    """
    points=np.asarray(xyz_m,dtype=float)
    if points.ndim!=2 or points.shape[1]!=3 or not np.all(np.isfinite(points)):
        raise ValueError('Coordinates must be a finite (N,3) array in metres.')
    if not all(np.isfinite(v) for v in (side_m,diameter_m,phase_x_m)) or side_m<=0 or diameter_m<=0:
        raise ValueError('Positive finite side and diameter and finite phase required.')
    centres=side_m*np.array([[0,0,0],[0,.5,.5],[.5,0,.5],[.5,.5,0]],dtype=float)
    centres[:,0]+=phase_x_m
    displacement=points[:,None,:]-centres[None,:,:]
    displacement-=side_m*np.round(displacement/side_m)
    distance=np.linalg.norm(displacement,axis=2);nearest=np.argmin(distance,axis=1)
    vector=displacement[np.arange(len(points)),nearest];radius=distance[np.arange(len(points)),nearest]
    normal=np.divide(vector,radius[:,None],out=np.zeros_like(vector),where=radius[:,None]>0)
    return radius-diameter_m/2,normal,points-vector


def location(row, *, side_m, dx_m, phase_x_m):
    xyz=np.array([[float(row[k+'_m']) for k in ('x','y','z')]])
    gap,normal,centre=sphere_distance(xyz,side_m=side_m,phase_x_m=phase_x_m)
    return dict(ijk=[int(row[k]) for k in ('i','j','k')],xyz_m=xyz[0].tolist(),
        concentration_kg_m3=float(row['concentration_kg_m3']),
        gap_to_analytical_surface_m=float(gap[0]),gap_over_dx=float(gap[0]/dx_m),
        nearest_sphere_centre_image_m=centre[0].tolist(),outward_normal=normal[0].tolist(),
        downstream_pole_cosine=float(normal[0,0]))
