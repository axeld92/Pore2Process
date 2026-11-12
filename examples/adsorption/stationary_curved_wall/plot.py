"""Render static comparisons and a compact record from verified results."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import numpy as np

from pore2process_demo.io import sha,write_json
from pore2process_demo.paths import ROOT
from examples.adsorption.transport_accuracy.run import validate_record
from examples.adsorption.transport_verification.analysis import load
from . import reference as ref
from .run import HERE,OUTPUT


def fields(base,record,output):
    raw,mask,_,f,_,q=load(base/'native/shear128')
    s,modes=ref.read_mode(base/'reference'/record['reference']['reference_file']);mode=(s,modes[-1])
    fine=(np.arange(420)+.5)*ref.L/420
    x,y=np.meshgrid(fine,fine,indexing='ij');r=np.hypot(x-ref.L/2,y-ref.L/2)
    exact=np.where(r<=ref.R,ref.profile(x,y,ref.TIMES[-1],mode=mode),np.nan)
    p=ref.profile(q[0][:,None],q[1][None,:],ref.TIMES[-1],mode=mode)
    actual=np.where(mask[:,:,0]==1,f['final'][:,:,0],np.nan);diff=actual-p
    fig,axes=plt.subplots(1,3,figsize=(11.8,4.05),layout='constrained')
    cmap=plt.get_cmap('coolwarm').copy();cmap.set_bad('#535a60');side=ref.L*1e3
    cspan=max(.25,float(np.nanmax(abs(actual-.5))));espan=max(float(np.nanmax(abs(diff))),1e-8)
    for i,(ax,a,title) in enumerate(zip(axes,(exact,actual,diff),('Julia radial reference','OpenLB TRT','OpenLB − reference'))):
        limit=espan if i==2 else cspan;center=0 if i==2 else .5
        image=ax.imshow(a.T,extent=(0,side,0,side),origin='lower',cmap=cmap,
                        vmin=center-limit,vmax=center+limit,interpolation='nearest')
        ax.add_patch(Circle((side/2,side/2),ref.R*1e3,fill=False,ec='black',lw=.6))
        lo,hi=side/2-.44,side/2+.44
        ax.set(title=title,xlabel='x (mm)',ylabel='y (mm)',xlim=(lo,hi),ylim=(lo,hi))
        if i==1:fig.colorbar(image,ax=list(axes[:2]),shrink=.8,label='Concentration (kg/m³)')
        if i==2:fig.colorbar(image,ax=ax,shrink=.8,label='Difference (kg/m³)')
    fig.savefig(output,dpi=175);plt.close(fig)


def history(directory,amplitude):
    a=np.genfromtxt(directory/'inventory.csv',delimiter=',',names=True)
    error=np.maximum(0,np.maximum(.5-amplitude-a['min_concentration_kg_m3'],a['max_concentration_kg_m3']-(.5+amplitude)))
    return a['time_s'],error


def comparison(base,record,output):
    fig,ax=plt.subplots(2,2,figsize=(10.4,7.5),layout='constrained');ns=[32,64,96,128]
    cases=record['cases']
    for key,ls in [('middle','--'),('final','-')]:
        for metric,color,label in [('RMS_over_amplitude','#2166ac','All liquid'),('shell_RMS_over_amplitude','#b2182b','Outer 50 µm')]:
            ax[0,0].plot(ns,[100*cases[f'shear{n}']['fields'][key][metric] for n in ns],marker='o',ls=ls,color=color,label=f'{label}, {key}')
    ax[0,0].axhline(1,ls=':',color='black',label='1% screen')
    ax[0,0].set(yscale='log',xlabel='Cells across square',ylabel='RMS (% of initial amplitude)',title='Concentration field error')
    for n in ns:
        t,e=history(base/'native'/f'shear{n}',ref.A);ax[0,1].plot(t,e,label=f'N{n}')
    ax[0,1].axhline(.02,color='black',ls=':',label='Allowance')
    ax[0,1].set(xlabel='Time (s)',ylabel='Excursion (kg/m³)',title='Transient bounds for varying tracer')
    for n in (64,128):
        t,e=history(base/'native'/f'constant{n}',0.);ax[1,0].plot(t,e,label=f'Stationary N{n}')
    if 'rigid_comparison' in record:
        cmp=record['rigid_comparison'];parent=ROOT/cmp['parent_summary']
        t,e=history(parent.parent/'native/constant128',0.);ax[1,0].plot(t,e,label='Rotating N128',ls='--')
    ax[1,0].set(xlabel='Time (s)',ylabel='Maximum departure (kg/m³)',title='Uniform tracer drift')
    for n in ns:
        x=cases[f'shear{n}']['fields']['final']
        ax[1,1].plot(n,100*x['link_RMS_over_amplitude'],'o',color='#b2182b',label='Wall-link cells' if n==32 else None)
        ax[1,1].plot(n,100*x['interior_RMS_over_amplitude'],'s',color='#2166ac',label='Interior, outside shell' if n==32 else None)
    ax[1,1].axhline(1,color='black',ls=':',label='1% screen')
    ax[1,1].set(yscale='log',xlabel='Cells across square',ylabel='RMS (% of initial amplitude)',title='Final error by region')
    for a in ax.flat:a.legend(fontsize=8)
    fig.savefig(output,dpi=175);plt.close(fig)


def main():
    p=argparse.ArgumentParser();p.add_argument('--summary',type=Path);a=p.parse_args()
    path=(a.summary or ROOT/json.loads((OUTPUT/'latest.json').read_text())['summary']).resolve();s=validate_record(path)
    if 'rigid_comparison' in s:
        cmp=s['rigid_comparison'];parent=ROOT/cmp['parent_summary']
        if sha(parent)!=cmp['parent_summary_sha256']:raise ValueError('Changed rotating-wall record')
        validate_record(parent)
    out=HERE/'figures';out.mkdir(exist_ok=True)
    fields(path.parent,s,out/'fields.png');comparison(path.parent,s,out/'comparison.png')
    compact={k:s[k] for k in ('reference','cases','new_native_runs','depth_comparison','background','new_bed_runs','new_DWSIM_runs','application_transport_qualified')}
    compact['reference']={k:v for k,v in compact['reference'].items() if k!='commands'}
    if 'rigid_comparison' in s:compact['rigid_comparison']=s['rigid_comparison']
    compact.update(summary=str(path.relative_to(ROOT)),summary_sha256=sha(path))
    write_json(HERE/'reference_summary.json',compact)
    write_json(OUTPUT/'figures.json',dict(summary_sha256=sha(path),plot_source_sha256=sha(Path(__file__)),
        outputs={str(p.relative_to(ROOT)):sha(p) for p in out.glob('*.png')},visual_inspection='pending'))
    print('Saved',out)


if __name__=='__main__':main()
