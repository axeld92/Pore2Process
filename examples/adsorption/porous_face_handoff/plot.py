"""Plot saved throughput, bounds and tracer-spreading diagnostics."""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from examples.adsorption.stationary_curved_wall import reference as ref
from examples.adsorption.transport_accuracy.run import validate_record
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT

from .run import HERE, OUTPUT, SNAPSHOTS


def profiles(center_path, center, bounded_path, bounded):
    fig, axes = plt.subplots(1,3,figsize=(13,3.9),layout='constrained')
    blue, green, orange = '#35698b', '#18766e', '#b35635'
    case = center['analysis']['cases']['96']
    native_q = case['source']['carrier_directional_flow_m3_s'][0]
    x = np.arange(1,97)/96
    for key,color,label in [('before',blue,'Interpolated'),('after',green,'Balanced')]:
        planes = np.array(case['flow'][key]['plane_flows_m3_s'][0])
        axes[0].plot(x,100*(planes/native_q-1),color=color,label=label)
    axes[0].set(title='Axial throughput, N96',xlabel='Axial position / cell length',
                ylabel='Departure from OpenLB mean (%)')
    axes[0].legend(fontsize=8)
    paths = [center_path.parent/'N96/varying_history.csv',
             bounded_path.parent/'N96/varying/varying_history.csv']
    for path,color,label in zip(paths,(orange,green),('Centered','Flux limited'),strict=True):
        h = np.genfromtxt(path,delimiter=',',names=True)
        axes[1].fill_between(h['time_s'],h['min_kg_m3'],h['max_kg_m3'],color=color,alpha=.14)
        axes[1].plot(h['time_s'],h['max_kg_m3'],color=color,label=label)
        axes[1].plot(h['time_s'],h['min_kg_m3'],color=color)
    for bound in (.25,.75):
        axes[1].axhline(bound,color='#333333',ls='--',lw=.8)
    axes[1].set(title='Concentration range, N96',xlabel='Time (s)',ylabel='Concentration (kg/m³)')
    axes[1].legend(fontsize=8)
    for n,color in ((64,blue),(96,green)):
        metrics = bounded['analysis']['cases'][f'N{n}_varying']['snapshots']
        axes[2].plot((0.,*ref.TIMES),(100.,*[100*metrics[k]['variance_retained'] for k in SNAPSHOTS]),
                     marker='o',color=color,label=f'N{n}')
    axes[2].set(title='Flux-limited spreading',xlabel='Time (s)',ylabel='Concentration variance retained (%)')
    axes[2].legend(fontsize=8)
    folder = HERE/'figures'
    folder.mkdir(exist_ok=True)
    fig.savefig(folder/'handoff.png',dpi=175)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--centered',type=Path)
    parser.add_argument('--bounded',type=Path)
    args=parser.parse_args()
    center=(args.centered or ROOT/read_json(OUTPUT/'latest.json')['summary']).resolve()
    bounded=(args.bounded or ROOT/read_json(OUTPUT/'latest_bounded.json')['summary']).resolve()
    a,b=validate_record(center),validate_record(bounded)
    if b['centered_summary_sha256']!=sha(center):
        raise ValueError('Wrong centered parent')
    profiles(center,a,bounded,b)
    write_json(HERE/'reference_summary.json',dict(
        centered_summary=str(center.relative_to(ROOT)),centered_summary_sha256=sha(center),
        bounded_summary=str(bounded.relative_to(ROOT)),bounded_summary_sha256=sha(bounded),
        centered=a['analysis'],bounded=b['analysis'],new_openlb_runs=0,new_bed_runs=0,new_DWSIM_runs=0))
    print('Saved',HERE/'figures/handoff.png')


if __name__=='__main__':
    main()
