"""Plot late pulse recovery and cumulative exit-time variance on both meshes."""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from examples.adsorption.transport_accuracy.run import validate_record
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT

from .n384 import HERE, PARENT, RESULTS
from .run import history


def plot(path,record):
    a = record['analysis']
    fig,(tail,variance) = plt.subplots(1,2,figsize=(9.4,3.6),layout='constrained')
    for n,folder,color in ((288,PARENT.parent/'N288','#39759e'),(384,path.parent/'N384','#167963')):
        h = history(folder/'history.csv')
        stats = a[f'N{n}']
        time = 1e3*h['time_s']
        recovered = h['recovered_kg']
        f = recovered/h['injected_kg'][-1]
        tail.semilogy(time,1-f,color=color,label=f'N{n}')
        if a['tail'][str(n)]['time_to_99pct_s'] is not None:
            t99 = 1e3*a['tail'][str(n)]['time_to_99pct_s']
            tail.plot(t99,.01,'o',color=color,ms=4)
        mu = stats['mean_exit_s']
        central = (h['second_moment_kg_s2']-2*mu*h['first_moment_kg_s']+mu**2*recovered)/recovered[-1]
        variance.plot(time,1e3*central,color=color,label=f'N{n}')
    tail.axhline(.001,ls='--',lw=1,color='#444444',label='99.9% recovered')
    tail.set(title='Late tracer recovery',xlabel='Time from pulse onset (ms)',
             ylabel='Fraction of injected tracer remaining',xlim=(150,None))
    variance.axvline(300,ls='--',lw=1,color='#444444')
    variance.set(title='Accumulated exit-time variance',xlabel='Time from pulse onset (ms)',
                 ylabel='Contribution to variance (10⁻³ s²)',xlim=(0,None))
    for ax in (tail,variance):
        ax.legend(fontsize=8)
    folder = HERE/'figures'
    folder.mkdir(exist_ok=True)
    fig.savefig(folder/'n384_tail.png',dpi=180)
    plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--summary',type=Path)
    args=parser.parse_args()
    path=(args.summary or ROOT/read_json(RESULTS/'latest.json')['summary']).resolve()
    record=validate_record(path)
    if record['parent_summary_sha256']!=sha(PARENT):
        raise ValueError('Wrong N288 comparison')
    plot(path,record)
    write_json(HERE/'n384_reference_summary.json',dict(summary=str(path.relative_to(ROOT)),summary_sha256=sha(path),
        mesh=record['mesh'],analysis=record['analysis'],new_openlb_runs=0,new_bed_runs=0,new_DWSIM_runs=0))
    print('Saved',HERE/'figures/n384_tail.png')


if __name__=='__main__':
    main()
