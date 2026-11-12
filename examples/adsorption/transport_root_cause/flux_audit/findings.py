"""Recompute local diagnostic findings from the preserved actual native files."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from pore2process_demo.io import sha, write_json
from pore2process_demo.paths import ROOT
from .analysis import summarize, summarize_uniform_stencil

HERE=Path(__file__).parent
BASE=ROOT/'results/adsorption/transport_root_cause/runs'


def peak_statistics(path: Path, *, side_m: float, particle_radius_m: float) -> dict:
    """Read-only annotations of the declared smooth spheres and voxel mask.

    The 26-neighbour geometric test labels one-cell wall proximity only; this
    is not a recreated transport descriptor or numerical boundary scheme.
    Distances/normals annotate the specified spheres, not a reconstructed
    experimental specimen. Ties can occur away from the extrema.
    """
    field=np.genfromtxt(path,delimiter=',',names=True)
    concentration=field['concentration_kg_m3'];fluid=field['material']!=2
    n=int(np.max(field['i']))+1;dx=side_m/n
    i,j,k=(field[name].astype(int) for name in ('i','j','k'))
    geometry=np.empty((n,n,n),dtype=int);geometry[i,j,k]=field['material'].astype(int)
    adjacent=np.zeros(len(field),dtype=int)
    for a in (-1,0,1):
        for b in (-1,0,1):
            for c in (-1,0,1):
                if a or b or c:adjacent+=(geometry[(i+a)%n,(j+b)%n,(k+c)%n]==2)
    points=np.column_stack([field[name] for name in ('x_m','y_m','z_m')])
    centres=np.array([(a*side_m/2,b*side_m/2,c*side_m/2) for a in range(-1,4)
                      for b in range(-1,4) for c in range(-1,4) if (a+b+c)%2==0])
    difference=points[:,None,:]-centres[None,:,:]
    distance=np.linalg.norm(difference,axis=2);nearest=np.argmin(distance,axis=1)
    norm=distance[np.arange(len(field)),nearest]
    normals=difference[np.arange(len(field)),nearest,:]/norm[:,None]
    gap=norm-particle_radius_m
    excursion=np.maximum(.25-concentration,concentration-.75)

    def group(mask):
        count=int(np.count_nonzero(mask))
        return dict(nodes=count,wall_adjacent_nodes=int(np.count_nonzero(mask&(adjacent>0))),
                    bulk_nodes=int(np.count_nonzero(mask&(adjacent==0))),
                    gap_over_dx_quantiles=[float(x) for x in np.quantile(gap[mask]/dx,[0,.25,.5,.75,1])] if count else None,
                    downstream_pole_nodes=int(np.count_nonzero(mask&(normals[:,0]>.9))),
                    upstream_pole_nodes=int(np.count_nonzero(mask&(normals[:,0]<-.9))))

    def point(index):
        return dict(node_id_ijk=[int(field[name][index]) for name in ('i','j','k')],
                    coordinates_m=points[index].tolist(),concentration_kg_m3=float(concentration[index]),
                    speed_m_s=float(np.linalg.norm([field[name][index] for name in ('ux_m_s','uy_m_s','uz_m_s')])),
                    voxel_solid_neighbours=int(adjacent[index]),nearest_sphere_centre_m=centres[nearest[index]].tolist(),
                    smooth_surface_gap_m=float(gap[index]),gap_over_dx=float(gap[index]/dx),
                    sphere_outward_normal=normals[index].tolist())
    ids=np.flatnonzero(fluid)
    return dict(file=str(path.relative_to(ROOT)),sha256=sha(path),n=n,
                fluid=group(fluid),engineering_violations=group(fluid&(excursion>.02)),
                severe_violations_over_0p1_kg_m3=group(fluid&(excursion>.1)),
                minimum=point(ids[np.argmin(concentration[fluid])]),maximum=point(ids[np.argmax(concentration[fluid])]),
                caveat='Peak localization is association, not proof of the first origin. Smooth solid spheres are convex; narrow void passages and voxel corners are separate geometry notions.')


def run():
    constant=BASE/'20261110T205538Z-4e56ceb6/native/n32_constant'
    raw=json.loads((constant/'raw.json').read_text())
    evidence=[constant/name for name in ('raw.json','uniform_stencil.csv','audit_1.csv','audit_2.csv')]
    result=dict(schema='pore2process.root-cause-population-findings.v1',
                constant_stencil=summarize_uniform_stencil(constant/'uniform_stencil.csv',
                    first_step_audit=constant/'audit_1.csv',initial_concentration=.5),
                constant_actual_first_step=summarize(constant/'audit_1.csv',dx_m=raw['dx_m'],dt_s=raw['dt_s']),
                peaks={})
    for name in ('n32_sine','n48_sine'):
        path=BASE/'20261110T205424Z-ed34b563/native'/name/'peak.csv'
        if path.exists():
            result['peaks'][name]=peak_statistics(path,side_m=.0013891725543141453,particle_radius_m=.0004)
            evidence += [path,path.parent/'peak.json']
    peak=BASE/'20261110T205757Z-c62abc03/native/n32_sine_audit4145'
    result['actual_n32_peak_population_audit']=summarize(peak/'audit_4145.csv',dx_m=raw['dx_m'],dt_s=raw['dt_s'])
    evidence += [peak/'audit_4145.csv',peak/'raw.json']
    result['evidence_sha256']={str(path.relative_to(ROOT)):sha(path) for path in evidence}
    write_json(HERE/'findings.json',result)
    print('Saved',HERE/'findings.json')
    return result


if __name__=='__main__':run()
