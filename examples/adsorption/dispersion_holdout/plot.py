"""Plot the one-cell fit and two-cell held-out cumulative outlet curves."""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from examples.adsorption.porous_outlet_tracer.run import history
from examples.adsorption.transport_accuracy.run import validate_record
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT

from .run import HERE, OUTPUT, PARENT, normalized


def render(path,record):
    base = path.parent
    dose = record['request']['dose_kg']
    fig,axes = plt.subplots(1,2,figsize=(9.6,3.7),layout='constrained')
    data = (history(PARENT.parent/'N384/history.csv'),history(base/'holdout/history.csv'))
    for i,(ax,measured) in enumerate(zip(axes,data,strict=True),start=1):
        model = history(base/f'model_{i}.csv')
        ax.plot(1000*model['time_s'],normalized(model,dose),ls='--',color='#b05a31',label='Frozen dispersion model')
        ax.plot(1000*measured['time_s'],normalized(measured,dose),color='#237b78',label='Pore tracer')
        ax.set(title='One-cell calibration response' if i==1 else 'Two-cell held-out response',
               xlabel='Time from pulse onset (ms)',ylabel='Fraction of injected tracer recovered',ylim=(0,1.02))
        ax.set_xlim(0,400 if i==1 else max(500,min(900,1000*measured['time_s'][-1])))
        if i==2 and record['analysis']['status']=='profile_rejected':
            ax.axvline(1000*measured['time_s'][-1],color='#4c4c4c',lw=1,label='Rejection reached')
        ax.legend(fontsize=8)
    folder = HERE/'figures'
    folder.mkdir(exist_ok=True)
    fig.savefig(folder/'outlet_holdout.png',dpi=180)
    plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--summary',type=Path)
    args=parser.parse_args()
    path=(args.summary or ROOT/read_json(OUTPUT/'latest.json')['summary']).resolve()
    record=validate_record(path)
    render(path,record)
    write_json(HERE/'reference_summary.json',dict(summary=str(path.relative_to(ROOT)),summary_sha256=sha(path),
        model=read_json(path.parent/'frozen_prediction.json')['model'],analysis=record['analysis']))
    print('Saved',HERE/'figures/outlet_holdout.png')


if __name__=='__main__':
    main()
