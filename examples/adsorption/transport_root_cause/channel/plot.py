"""Plot independently verified periodic shear controls."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from pore2process_demo.paths import ROOT
from pore2process_demo.io import sha,write_json
from .run import OUTPUT,HERE,check
from .reference import concentration


def main():
    verification=check()
    summary=ROOT/verification['source_summary'];record=json.loads(summary.read_text());base=summary.parent
    fig,axes=plt.subplots(2,2,figsize=(11,8),layout='constrained')
    for n,color in [(32,'#0072B2'),(48,'#D55E00'),(64,'#009E73')]:
        folder=base/'native'/f'shear{n}';a=np.genfromtxt(folder/'inventory.csv',delimiter=',',names=True)
        axes[0,0].plot(a['time_s'],a['min_concentration_kg_m3'],color=color,label=f'N={n}')
        axes[0,0].plot(a['time_s'],a['max_concentration_kg_m3'],color=color)
    for value in (.25,.75):axes[0,0].axhline(value,color='.3',ls='--',lw=1)
    axes[0,0].set(xlabel='Time (s)',ylabel='Minimum / maximum concentration (kg m⁻³)',title='Periodic shear: bounds throughout the run')
    axes[0,0].legend()
    ns=np.array([32,48,64]);cases=record['cases']
    axes[0,1].loglog(ns,[cases[f'shear{n}']['RMS_error_over_amplitude'] for n in ns],'o-',label='RMS / initial amplitude')
    axes[0,1].loglog(ns,[cases[f'shear{n}']['maximum_error_over_amplitude'] for n in ns],'s-',label='Maximum / initial amplitude')
    axes[0,1].axhline(.01,color='.3',ls='--',label='1% RMS comparison')
    axes[0,1].set(xlabel='Cells across the periodic cube',ylabel='Relative error',title='Independent continuum comparison at final time')
    axes[0,1].set_xticks(ns,labels=[str(n) for n in ns]);axes[0,1].legend()
    axes[0,1].xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    for ax,n in zip(axes[1],(32,64)):
        f=np.genfromtxt(base/'native'/f'shear{n}'/'final.csv',delimiter=',',names=True)
        raw=cases[f'shear{n}']['raw']
        physical=dict(side_m=raw['side_m'],speed_m_s=raw['speed_m_s'],shear_fraction=raw['shear_fraction'],
                      diffusivity_m2_s=raw['diffusivity_m2_s'],duration_s=raw['duration_s'])
        exact=concentration(f['x_m'],f['y_m'],**physical)
        field=(f['concentration_kg_m3']-exact).reshape(n,n).T
        im=ax.imshow(field,origin='lower',extent=(0,raw['side_m']*1000,0,raw['side_m']*1000),cmap='RdBu_r',aspect='equal')
        ax.set(xlabel='x (mm)',ylabel='y (mm)',title=f'N={n}: OpenLB − continuum (kg m⁻³)')
        fig.colorbar(im,ax=ax,label='Concentration difference (kg m⁻³)')
    dest=HERE/'figures';dest.mkdir(exist_ok=True);output=dest/'periodic_shear.png';fig.savefig(output,dpi=160);plt.close(fig)
    write_json(OUTPUT/'figure_sources.json',dict(source_summary=str(summary.relative_to(ROOT)),source_summary_sha256=sha(summary),
        figure=str(output.relative_to(ROOT)),figure_sha256=sha(output),plot_sha256=sha(Path(__file__)),
        verification=verification['run_id']))
    print(output)

if __name__=='__main__':main()
