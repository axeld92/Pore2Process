"""Combine guarded full/thin shear refinement and the independent reference."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from pore2process_demo.paths import ROOT
from pore2process_demo.io import sha,write_json
from ..reference import concentration
from .run import HERE,OUTPUT,check


def main():
    verification=check();path=ROOT/verification['summary'];s=json.loads(path.read_text());base=path.parent
    parent_path=ROOT/s['parent_summary'];parent=json.loads(parent_path.read_text())
    ns=[32,48,64,96,128]
    cases=[parent['cases'][f'shear{n}'] for n in ns[:3]]+[s['cases'][f'thin{n}'] for n in ns[3:]]
    fig,axes=plt.subplots(2,2,figsize=(11,8),layout='constrained')
    axes[0,0].loglog(ns,[c['RMS_error_over_amplitude'] for c in cases],'o-',label='Matched dt; thin z above N64')
    axes[0,0].loglog([128],[s['cases']['thin128_half_dt']['RMS_error_over_amplitude']],'s',label='N128 half dt')
    axes[0,0].axhline(.01,color='.3',ls='--',label='1% profile tolerance')
    axes[0,0].set(xlabel='Cells across x/y',ylabel='Final RMS / initial amplitude',title='Bulk-model error versus continuum reference')
    axes[0,0].set_xticks(ns,labels=[str(n) for n in ns]);axes[0,0].legend()
    axes[0,0].xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    axes[0,1].plot(ns,[c['maximum_bound_excursion_kg_m3'] for c in cases],'o-',label='Matched dt; thin z above N64')
    axes[0,1].plot([128],[s['cases']['thin128_half_dt']['maximum_bound_excursion_kg_m3']],'s',label='N128 half dt')
    axes[0,1].axhline(.02,color='.3',ls='--',label='Engineering allowance')
    axes[0,1].axhline(0,color='.5',lw=.7)
    axes[0,1].set(xlabel='Cells across x/y',ylabel='Largest bound excursion (kg m⁻³)',title='Concentration bounds throughout the run')
    axes[0,1].legend()
    f=np.genfromtxt(base/'native/thin128/final.csv',delimiter=',',names=True);raw=s['cases']['thin128']['raw'];n=128
    exact=concentration(f['x_m'],f['y_m'],side_m=raw['side_m'],speed_m_s=raw['speed_m_s'],shear_fraction=raw['shear_fraction'],
                        diffusivity_m2_s=raw['diffusivity_m2_s'],duration_s=raw['duration_s'])
    for ax,values,title in [(axes[1,0],f['concentration_kg_m3'],'N128 OpenLB concentration'),(axes[1,1],f['concentration_kg_m3']-exact,'N128 OpenLB − continuum')]:
        im=ax.imshow(values.reshape(n,n).T,origin='lower',extent=(0,raw['side_m']*1000,0,raw['side_m']*1000),aspect='equal',cmap='viridis' if ax==axes[1,0] else 'RdBu_r')
        ax.set(xlabel='x (mm)',ylabel='y (mm)',title=title);fig.colorbar(im,ax=ax,label='Concentration (kg m⁻³)')
    dest=HERE/'figures';dest.mkdir(exist_ok=True);output=dest/'shear_refinement.png';fig.savefig(output,dpi=160);plt.close(fig)
    write_json(OUTPUT/'figure_sources.json',dict(summary=str(path.relative_to(ROOT)),summary_sha256=sha(path),
        parent_summary=str(parent_path.relative_to(ROOT)),parent_summary_sha256=sha(parent_path),
        figure=str(output.relative_to(ROOT)),figure_sha256=sha(output),plot_sha256=sha(Path(__file__))))
    print(output)

if __name__=='__main__':main()
