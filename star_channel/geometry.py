"""Declare geometry and write a watertight STL; OpenLB owns all geometry queries."""
import math
from pathlib import Path

import numpy as np

def polygon(outer=0.0023, inner=0.0014, points=6):
    if not all(math.isfinite(v) and v > 0 for v in (outer, inner)) or inner >= outer or points != 6:
        raise ValueError('Require finite 0 < inner radius < outer radius and six star points.')
    angles = np.arange(2*points)*math.pi/points
    r = np.where(np.arange(2*points) % 2 == 0, outer, inner)
    return np.column_stack((r*np.cos(angles), r*np.sin(angles)))

def dimensions(vertices):
    v = np.asarray(vertices, dtype=float)
    w = np.roll(v, -1, axis=0)
    area = float(abs(np.sum(v[:,0]*w[:,1]-v[:,1]*w[:,0]))/2)
    perimeter = float(np.sum(np.linalg.norm(w-v, axis=1)))
    return {'open_area_m2': area, 'wetted_perimeter_m': perimeter, 'hydraulic_diameter_m': 4*area/perimeter}

def write_stl(path, vertices, xmin, xmax, offset):
    """Outward triangles for a star-shaped extrusion along x; SI coordinates."""
    v = np.asarray(vertices)+offset
    center = np.array([offset, offset])
    triangles = []
    for a,b in zip(v, np.roll(v,-1,axis=0)):
        p,q,r,s = [np.array(t) for t in ((xmin,*a),(xmax,*a),(xmax,*b),(xmin,*b))]
        triangles.extend(((p,s,r),(p,r,q)))
        triangles.append((np.array([xmin,*center]),s,p))
        triangles.append((np.array([xmax,*center]),q,r))
    with Path(path).open('x') as f:
        f.write('solid explicit_star\n')
        for triangle in triangles:
            a,b,c = triangle
            normal = np.cross(b-a,c-a); normal /= np.linalg.norm(normal)
            f.write(' facet normal '+' '.join(f'{x:.17g}' for x in normal)+'\n  outer loop\n')
            for point in triangle:
                f.write('   vertex '+' '.join(f'{x:.17g}' for x in point)+'\n')
            f.write('  endloop\n endfacet\n')
        f.write('endsolid explicit_star\n')
