"""Static scientific plots from checked native flat-wall data and its oracle."""
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from pore2process_demo.paths import ROOT
from pore2process_demo.io import sha,write_json
from .run import HERE,OUTPUT
from .continuation import check
from .reference_qualified import converged_reference


def main():
    verified=check();summary=ROOT/verified['source_summary'];s=json.loads(summary.read_text());base=summary.parent
    names=['thin32','thin64','thin96'];n=[s['cases'][name]['raw']['n'] for name in names]
    rms=[100*s['cases'][name]['RMS_error_over_amplitude'] for name in names]
    raw=s['cases']['thin96']['raw'];field_path=ROOT/s['case_directories']['thin96']/'final.csv'
    f=np.genfromtxt(field_path,delimiter=',',names=True)
    exact,reference=converged_reference(f['x_m'],f['y_m'],height_m=raw['height_m'],
        mean_speed_m_s=raw['mean_speed_m_s'],diffusivity_m2_s=raw['diffusivity_m2_s'],duration_s=raw['duration_s'])
    difference=f['concentration_kg_m3']-exact;h=raw['height_m']*1e3
    fig,axes=plt.subplots(1,3,figsize=(13,3.6),layout='constrained')
    axes[0].semilogy(n,rms,'o-',label='Native versus continuum')
    axes[0].axhline(1.,color='0.4',ls='--',label='1% RMS screen')
    axes[0].set(xlabel='Fluid cells across channel, N',ylabel='RMS error / initial amplitude (%)',title='Resolution at fixed physical conditions')
    axes[0].set_xticks(n);axes[0].legend(fontsize=8)
    bound=float(np.max(abs(difference)));im=axes[1].imshow(difference.reshape((96,96)).T,origin='lower',extent=(0,h,0,h),cmap='coolwarm',vmin=-bound,vmax=bound,aspect='equal')
    axes[1].axhline(0,color='k',lw=2);axes[1].axhline(h,color='k',lw=2)
    axes[1].set(xlabel='x (mm)',ylabel='y (mm)',title='N96: numerical minus continuum')
    fig.colorbar(im,ax=axes[1],label='Concentration error (kg/m³)',shrink=.8)
    mask=f['j']==0
    axes[2].plot(f['x_m'][mask]*1e3,exact[mask],'k-',label='Qualified continuum oracle')
    axes[2].plot(f['x_m'][mask]*1e3,f['concentration_kg_m3'][mask],'o',ms=3,label='Stock OpenLB')
    axes[2].set(xlabel='x (mm)',ylabel='Concentration (kg/m³)',title='N96: first fluid row next to wall')
    axes[2].legend(fontsize=8)
    image=HERE/'flat_wall_diagnosis.png';fig.savefig(image,dpi=180);plt.close(fig)
    write_json(OUTPUT/'figure_sources.json',dict(source_summary=str(summary.relative_to(ROOT)),source_summary_sha256=sha(summary),
        final_field=str(field_path.relative_to(ROOT)),final_field_sha256=sha(field_path),oracle=reference,
        image=str(image.relative_to(ROOT)),image_sha256=sha(image),new_native_runs=0,new_DWSIM_runs=0))
    print(image)


if __name__=='__main__':main()
