"""Plot tracer variance retention and the final variance budget."""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from examples.adsorption.porous_face_handoff.run import SNAPSHOTS
from examples.adsorption.stationary_curved_wall import reference as ref
from examples.adsorption.transport_accuracy.run import validate_record
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT

from .time_half import OUTPUT

HERE = Path(__file__).parent


def render(analysis):
    matched = analysis['mesh']['matched_case']
    styles = [('N96_original','N96 original step','#7b8290','--'),
              ('N96_half_step','N96 half step','#607da4',':'),
              (matched,'N96 matched step','#2b5c80','-'),
              ('N192','N192 matched step','#08796a','-')]
    styles = [style for i,style in enumerate(styles) if style[0] not in [s[0] for s in styles[:i]]]
    fig, (profile,budget) = plt.subplots(1,2,figsize=(9.6,3.9),layout='constrained')
    for name,label,color,style in styles:
        stats = analysis['cases'][name]['snapshots']
        profile.plot((0.,*ref.TIMES),(100.,*[100*stats[t]['variance_retained'] for t in SNAPSHOTS]),
                     label=label,color=color,ls=style,marker='o')
    profile.set(xlabel='Time (s)',ylabel='Concentration variance retained (%)',title='Varying tracer')
    profile.legend(fontsize=8)
    final = [analysis['variance_budget'][name]['final'] for name,_,_,_ in styles]
    x = np.arange(len(styles))
    molecular = [1000*a['molecular_loss_kg2_m6'] for a in final]
    numerical = [1000*a['numerical_loss_kg2_m6'] for a in final]
    budget.bar(x,molecular,color='#32758e',label='Molecular diffusion')
    budget.bar(x,numerical,bottom=molecular,color='#dc955c',label='Numerical remainder')
    budget.set_xticks(x,[label.replace(' ', '\n', 1) for _,label,_,_ in styles],fontsize=8)
    budget.set(title='Final loss of half-variance',ylabel='Loss (10⁻³ kg²/m⁶)',
               ylim=(0,1.24*max(m+n for m,n in zip(molecular,numerical,strict=True))))
    budget.legend(fontsize=8)
    out = HERE/'figures'
    out.mkdir(exist_ok=True)
    fig.savefig(out/'spreading.png',dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--summary',type=Path)
    args = parser.parse_args()
    path = (args.summary or ROOT/read_json(OUTPUT/'latest.json')['summary']).resolve()
    record = validate_record(path)
    render(record['analysis'])
    write_json(HERE/'reference_summary.json',dict(summary=str(path.relative_to(ROOT)),summary_sha256=sha(path),
        analysis=record['analysis'],new_openlb_runs=0,new_bed_runs=0,new_DWSIM_runs=0))
    print('Saved',HERE/'figures/spreading.png')


if __name__=='__main__':
    main()
