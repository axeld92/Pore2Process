"""Plot the verified alignment sensitivity without interpreting it as a design."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from pore2process_demo.paths import ROOT
from pore2process_demo.io import sha,write_json
from .run import HERE,OUTPUT,replay


def main():
    path=ROOT/json.loads((OUTPUT/'latest.json').read_text())['summary'];s=replay(path)
    baseline=ROOT/s['baseline_summary'];old=s['unshifted'];new=s['shifted']['shifted_sine']
    fig,axes=plt.subplots(1,2,figsize=(11,4.5),layout='constrained')
    for directory,color,label in [(baseline.parent/'native/n32_sine','#0072B2','Original grid phase'),
            (path.parent/'native/shifted_sine','#D55E00','Translated by half a cell')]:
        a=np.genfromtxt(directory/'inventory.csv',delimiter=',',names=True)
        axes[0].plot(a['time_s'],a['min_concentration_kg_m3'],color=color,label=label)
        axes[0].plot(a['time_s'],a['max_concentration_kg_m3'],color=color)
    for value in (.25,.75):axes[0].axhline(value,color='.4',ls='--',lw=1)
    axes[0].set(xlabel='Time (s)',ylabel='Minimum / maximum concentration (kg m⁻³)',title='Same continuum problem, different voxel alignment')
    axes[0].legend(fontsize=9)
    x=np.arange(2);width=.32
    for data,offset,color,label in [(old,-width/2,'#0072B2','Original phase'),(new,width/2,'#D55E00','Half-cell phase')]:
        gaps=[data['extrema'][key]['gap_to_analytical_surface_m']*1e6 for key in ('minimum','maximum')]
        bars=axes[1].bar(x+offset,gaps,width,color=color,label=label)
        for bar,key in zip(bars,('minimum','maximum')):
            axes[1].annotate(f"nₓ={data['extrema'][key]['downstream_pole_cosine']:.3f}",
                (bar.get_x()+bar.get_width()/2,bar.get_height()),xytext=(0,4),textcoords='offset points',ha='center',fontsize=8)
    axes[1].set_xticks(x,labels=['Minimum C location','Maximum C location'])
    axes[1].set(ylabel='Gap to analytical sphere surface (µm)',title='Locations in each saved peak-error field')
    axes[1].set_ylim(bottom=0);axes[1].legend(fontsize=9)
    dest=HERE/'figures';dest.mkdir(exist_ok=True);output=dest/'grid_phase.png';fig.savefig(output,dpi=160);plt.close(fig)
    record=dict(summary=str(path.relative_to(ROOT)),summary_sha256=sha(path),baseline_summary=str(baseline.relative_to(ROOT)),
                baseline_summary_sha256=sha(baseline),figure=str(output.relative_to(ROOT)),figure_sha256=sha(output),plot_sha256=sha(Path(__file__)))
    write_json(OUTPUT/'figure_sources.json',record);print(output)

if __name__=='__main__':main()
