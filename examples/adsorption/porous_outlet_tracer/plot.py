"""Plot the measured outlet pulse and its late mass recovery."""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from examples.adsorption.transport_accuracy.run import validate_record
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT

from .run import HERE, OUTPUT, history


def render(path,record):
    fig,(signal,tail)=plt.subplots(1,2,figsize=(9.7,3.75),layout='constrained')
    for n,color in ((192,'#39688d'),(288,'#147766')):
        h=history(path.parent/f'N{n}/history.csv')
        t=h['time_s']
        entered=h['injected_kg'][-1]
        label=f'N{n}'
        midpoint=(t[1:]+t[:-1])/2
        signal.plot(1e3*midpoint,h['interval_exit_kg'][1:]/(entered*np.diff(t)),color=color,label=label)
        recovered=h['recovered_kg']/entered
        tail.semilogy(1e3*t,1-recovered,color=color,label=label)
        t05=record['analysis']['cases'][str(n)]['early_arrival_s']['0.05']
        if t05 is not None:
            signal.axvline(t05*1e3,color=color,ls=':',lw=.9)
    span=record['meshes']['288']['transit_s']
    signal.set(title='Measured outlet pulse',xlabel='Time from injection onset (ms)',
               ylabel='Exit flux / injected mass (s⁻¹)',xlim=(0,3*span*1e3))
    tail.axhline(.001,color='#393939',ls='--',lw=1,label='99.9% recovered')
    tail.set(title='Tracer still to be recovered',xlabel='Time from injection onset (ms)',
             ylabel='Fraction of injected tracer',ylim=(None,1.2))
    for a in (signal,tail):
        a.legend(fontsize=8)
    folder=HERE/'figures'
    folder.mkdir(exist_ok=True)
    fig.savefig(folder/'outlet.png',dpi=180)
    plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--summary',type=Path)
    args=parser.parse_args()
    path=(args.summary or ROOT/read_json(OUTPUT/'latest.json')['summary']).resolve()
    record=validate_record(path)
    render(path,record)
    write_json(HERE/'reference_summary.json',dict(summary=str(path.relative_to(ROOT)),summary_sha256=sha(path),
        nominal_transits=record['nominal_transits'],meshes=record['meshes'],analysis=record['analysis'],
        new_openlb_runs=0,new_bed_runs=0,new_DWSIM_runs=0))
    print('Saved',HERE/'figures/outlet.png')


if __name__=='__main__':
    main()
