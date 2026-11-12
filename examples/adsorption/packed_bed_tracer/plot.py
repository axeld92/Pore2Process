"""Scientific plots from saved native signals and Julia predictions."""
import json
from pathlib import Path
import numpy as np
from pore2process_demo.paths import ROOT


def plot(base):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    base=Path(base);s=json.loads((base/'summary.json').read_text())
    fig,axes=plt.subplots(1,2,figsize=(11,4.5),layout='constrained',sharey=True)
    for ax,keys,model,title in [(axes[0],['bed16','bed24','bed32'],'bed32','Short segment: derive contact moments'),
                              (axes[1],['length32'],'length32','Twice the length: no refit')]:
        for key in keys:
            data=np.genfromtxt(base/'native'/key/'flux.csv',delimiter=',',names=True)
            p=s['cases'][key];ax.plot(data['time_s'],data['outlet_flux_kg_s']/p['moments']['steady_exit_flux_kg_s'],
                label=f'OpenLB, N={p["raw"]["n"]}',lw=1.4)
        if (base/f'julia_{model}.json').exists():
            r=json.loads((base/f'julia_{model}.json').read_text())
            ax.plot(r['result']['time_s'],np.array(r['result']['outlet_flux_per_area_m_s'])/
                r['analytical_moments']['steady_flux_per_area_m_s'],'k--',lw=1.6,label='Julia ADE, common coefficients')
        ax.set(xlabel='Time (s)',title=title,ylim=(-.02,1.03));ax.grid(alpha=.2);ax.legend(fontsize=8,loc='lower right')
        ax.set_xlim(0,s['cases'][keys[-1]]['moments']['mean_s']*2.5)
    axes[0].set_ylabel('Exit tracer flux / steady exit flux (1)')
    fig.suptitle('Synthetic FCC bed: inert-tracer contact, absorbing outlet')
    fig.savefig(base/'tracer_comparison.png',dpi=180);plt.close(fig)
    # A pore velocity section illustrates actual native data, without rendering
    # a fictional specimen or adding a second geometry generator.
    p=s['cases']['bed32']['raw'];a=np.genfromtxt(base/'native/bed32/final.csv',delimiter=',',names=True)
    selected=a[a['i']==p['n']//4];values=np.full((p['n'],p['n']),np.nan)
    for r in selected:
        if r['material']!=2:values[int(r['j']),int(r['k'])]=r['ux_m_s']*1000
    fig,ax=plt.subplots(figsize=(5.6,4.6),layout='constrained')
    im=ax.imshow(values.T,origin='lower',extent=(0,p['side_m']*1000,0,p['side_m']*1000),cmap='viridis')
    fig.colorbar(im,ax=ax,label='Native pore velocity along x (mm/s)')
    ax.set(xlabel='y (mm)',ylabel='z (mm)',title='FCC packing: flow through one pore section')
    fig.savefig(base/'pore_velocity.png',dpi=180);plt.close(fig)

def main():
    p=json.loads((ROOT/'results/adsorption/packed_bed_tracer/latest.json').read_text());plot((ROOT/p['summary']).parent)

if __name__=='__main__':main()
