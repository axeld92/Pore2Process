"""Static scientific figures from checked native peak fields and mesh results."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from pore2process_demo.paths import ROOT
from pore2process_demo.io import sha, write_json
from .study import HERE, OUTPUT


def main():
    source=OUTPUT/'assessment.json';a=json.loads(source.read_text())
    names=[k for k in ('n32_sine','n48_sine','n64_sine','n96_sine') if k in a['cases']]
    figures=HERE/'figures';figures.mkdir(exist_ok=True)
    rows,cols=(2,2) if len(names)==4 else (1,len(names))
    fig,axes=plt.subplots(rows,cols,figsize=(4.3*cols,4*rows),squeeze=False,layout='constrained')
    span=max([.4]+[abs(a['cases'][name]['result'][key]-.5) for name in names
        for key in ('minimum_concentration_kg_m3','maximum_concentration_kg_m3')])
    span=max(.4,np.ceil(20*span)/20)
    fields=[]
    for axis,name in zip(axes.flat,names):
        item=a['cases'][name];r=item['result'];p=ROOT/item['native_directory']/'peak.csv';fields.append(p)
        f=np.genfromtxt(p,delimiter=',',names=True);n=r['raw']['n'];C=f['concentration_kg_m3'].reshape((n,n,n));fluid=(f['material']!=2).reshape(C.shape)
        peak=r['peak'];point=peak['minimum'] if .25-peak['minimum']['concentration_kg_m3']>peak['maximum']['concentration_kg_m3']-.75 else peak['maximum']
        j=point['ijk'][1];s=np.where(fluid[:,j,:],C[:,j,:],np.nan)
        side=r['raw']['side_m']*1e3
        cmap=plt.get_cmap('coolwarm').copy();cmap.set_bad('#50555b')
        image=axis.imshow(s.T,origin='lower',extent=[0,side,0,side],cmap=cmap,vmin=.5-span,vmax=.5+span,interpolation='nearest')
        axis.set(title=f"N={n}; peak at {r['peak']['time_s']:.4f} s\nwhole-run range {r['minimum_concentration_kg_m3']:.3f}–{r['maximum_concentration_kg_m3']:.3f} kg/m³",xlabel='x (mm)',ylabel='z (mm)')
        x,_,z=point['xyz_m'];axis.plot(x*1e3,z*1e3,'ko',ms=5,mfc='none')
        axis.text(.02,.02,f'y = {(j+.5)*side/n:.3f} mm',transform=axis.transAxes,fontsize=9,bbox=dict(facecolor='white',alpha=.8,edgecolor='none'))
    fig.colorbar(image,ax=list(axes.flat),label='Concentration (kg/m³)',shrink=.8)
    fig.savefig(figures/'peak_fields.png',dpi=180,bbox_inches='tight',pad_inches=.15);plt.close(fig)
    fig,axes=plt.subplots(2,2,figsize=(10.2,7.2),layout='constrained');axes=axes.ravel()
    ns=[a['cases'][k]['result']['raw']['n'] for k in names]
    excursions=[a['cases'][k]['result']['peak']['maximum_excursion_kg_m3'] for k in names]
    axes[0].plot(ns,excursions,'o-');axes[0].axhline(.02,color='firebrick',ls='--',label='Chosen allowance')
    axes[0].set(xlabel='Cells per unit-cell side',ylabel='Largest excursion (kg/m³)',title='Extrema across meshes');axes[0].legend(fontsize=8)
    for name in names:
        item=a['cases'][name];p=ROOT/item['native_directory']/'moments.csv';fields.append(p)
        m=np.genfromtxt(p,delimiter=',',names=True)
        q=p.parent/'inventory.csv';fields.append(q);h=np.genfromtxt(q,delimiter=',',names=True)
        label=name.split('_')[0].upper()
        e=np.maximum(0.,np.maximum(.25-h['min_concentration_kg_m3'],h['max_concentration_kg_m3']-.75))
        axes[1].plot(h['time_s'],e,label=label)
        axes[2].plot(m['time_s'],m['sd_C_kg_m3'],label=label)
        axes[3].plot(m['time_s'],np.hypot(m['sin_coefficient_kg_m3'],m['cos_coefficient_kg_m3']),label=label)
    axes[1].axhline(.02,color='firebrick',ls='--',label='Chosen allowance')
    axes[1].set(xlabel='Time (s)',ylabel='Excursion (kg/m³)',title='Concentration bounds over time')
    axes[2].set(xlabel='Time (s)',ylabel='Liquid concentration SD (kg/m³)',title='Mixing in physical liquid')
    axes[3].set(xlabel='Time (s)',ylabel='First-harmonic projection (kg/m³)',title='Longitudinal wave evolution')
    for axis in axes[1:]:axis.legend(fontsize=8)
    fig.savefig(figures/'mesh_diagnostics.png',dpi=180,bbox_inches='tight',pad_inches=.15);plt.close(fig)
    write_json(OUTPUT/'figure_sources.json',dict(assessment=str(source.relative_to(ROOT)),assessment_sha256=sha(source),
        script_sha256=sha(Path(__file__)),input_sha256={str(p.relative_to(ROOT)):sha(p) for p in fields},
        output_sha256={str(p.relative_to(ROOT)):sha(p) for p in figures.glob('*.png')},
        visual_inspection='pending'))
    print('Saved',figures,flush=True)


if __name__=='__main__':main()
