"""Read-only N64 peak localization and exact-population comparison."""
import argparse
import json
from pathlib import Path

import numpy as np

from pore2process_demo.paths import ROOT
from pore2process_demo.io import sha,write_json
from .findings import peak_statistics
from .analysis import summarize

HERE=Path(__file__).parent
FULL=ROOT/'results/adsorption/transport_root_cause/parallel/runs/20261110T212155Z-fb52fcc9'
REPLAY=ROOT/'results/adsorption/transport_root_cause/replay/runs/20261110T213016Z-68e75fda'


def compute():
    evidence={}
    for base,paths in [(FULL,['native/n64_sine/peak.csv','native/n64_sine/peak.json','native/n64_sine/raw.json','native/n64_sine/inventory.csv']),
                       (REPLAY,['native/n64_sine_audit7510/audit_7510.csv','native/n64_sine_audit7510/raw.json'])]:
        summary=json.loads((base/'summary.json').read_text())
        evidence[str((base/'summary.json').relative_to(ROOT))]=sha(base/'summary.json')
        for name in paths:
            path=base/name;digest=sha(path)
            if summary['artifacts_sha256'][name]!=digest:raise ValueError('Changed native evidence: '+str(path))
            evidence[str(path.relative_to(ROOT))]=digest
    p=FULL/'native/n64_sine';r=json.loads((p/'raw.json').read_text());H=r['side_m']
    peak=peak_statistics(p/'peak.csv',side_m=H,particle_radius_m=.0004)
    audit_path=REPLAY/'native/n64_sine_audit7510/audit_7510.csv'
    balance=summarize(audit_path,dx_m=r['dx_m'],dt_s=r['dt_s'])
    f=np.genfromtxt(p/'peak.csv',delimiter=',',names=True)
    a=np.genfromtxt(audit_path,delimiter=',',names=True)
    for name in ('i','j','k'):np.testing.assert_array_equal(f[name],a[name])
    np.testing.assert_array_equal(f['material'],a['material'])
    replay_error=float(np.max(abs(f['concentration_kg_m3']-a['concentration_after_stream'])))
    if replay_error>1e-12:raise ValueError('Peak replay concentration differs')
    fluid=f['material']!=2;excursion=np.maximum(.25-f['concentration_kg_m3'],f['concentration_kg_m3']-.75)
    severe=fluid&(excursion>.1);negative=fluid&(a['negative_post_populations']>0)
    points=np.column_stack([f[name][negative] for name in ('x_m','y_m','z_m')])
    centres=np.array([(i*H/2,j*H/2,k*H/2) for i in range(-1,4) for j in range(-1,4) for k in range(-1,4) if(i+j+k)%2==0])
    diff=points[:,None,:]-centres[None,:,:];distance=np.linalg.norm(diff,axis=2);nearest=np.argmin(distance,axis=1)
    radial=distance[np.arange(len(points)),nearest];normal=diff[np.arange(len(points)),nearest,:]/radial[:,None]
    rows=[]
    for index,row_index in enumerate(np.flatnonzero(negative)):
        rows.append(dict(ijk=[int(f[name][row_index]) for name in ('i','j','k')],
            C_after_stream_kg_m3=float(f['concentration_kg_m3'][row_index]),
            min_unshifted_post_population=float(a['min_unshifted_post_population'][row_index]),
            solid_links=int(a['wall_adjacent_links'][row_index]),smooth_sphere_gap_m=float(radial[index]-.0004),
            gap_over_dx=float((radial[index]-.0004)/r['dx_m']),sphere_outward_normal=normal[index].tolist()))
    previous=json.loads((HERE/'findings.json').read_text())
    comparison={}
    for name,point in {**previous['peaks'],'n64_sine':peak}.items():
        maximum=point['maximum']
        comparison[name]=dict(maximum_C_kg_m3=maximum['concentration_kg_m3'],
            selected_peak_maximum_excursion_kg_m3=maximum['concentration_kg_m3']-.75,
            smooth_sphere_gap_m=maximum['smooth_surface_gap_m'],gap_over_dx=maximum['gap_over_dx'],
            normal_x=maximum['sphere_outward_normal'][0],speed_m_s=maximum['speed_m_s'],
            solid_neighbours=maximum['voxel_solid_neighbours'])
    h=np.genfromtxt(p/'inventory.csv',delimiter=',',names=True)
    evidence[str((HERE/'findings.json').relative_to(ROOT))]=sha(HERE/'findings.json')
    return dict(schema='pore2process.n64-peak-localization.v1',peak=peak,
        peak_record=json.loads((p/'peak.json').read_text()),exact_population_balance=balance,
        whole_history_minimum_C_kg_m3=float(np.min(h['min_concentration_kg_m3'])),
        whole_history_minimum_step=int(h['step'][np.argmin(h['min_concentration_kg_m3'])]),
        peak_replay_maximum_concentration_difference_kg_m3=replay_error,
        negative_population_nodes=rows,negative_node_severe_excursion_overlap=int(np.count_nonzero(negative&severe)),
        mesh_comparison=comparison,evidence_sha256=evidence,
        analysis_source_sha256={str(path.relative_to(ROOT)):sha(path) for path in [Path(__file__),HERE/'findings.py',HERE/'analysis.py']},
        new_native_runs=0,new_DWSIM_runs=0,
        interpretation='Analytic smooth-surface proximity and voxel wall/stagnation structure are correlated with the selected dominant peaks. The actual full-way wall is a voxel midlink, not the analytic sphere cut. These three meshes do not prove causality or a monotonic gap/error law; grid-phase sensitivity is a separate test.')


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--check',action='store_true');args=parser.parse_args()
    result=compute();path=HERE/'findings64.json'
    if args.check:
        if json.loads(path.read_text())!=result:raise ValueError('N64 findings do not replay')
        print('N64 findings replay; zero new native/DWSIM runs')
    else:write_json(path,result);print(path)


if __name__=='__main__':main()
