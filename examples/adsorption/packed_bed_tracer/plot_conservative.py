"""Absolute native exit flux and the selected Q/V-preserving Julia prediction."""
import json
import numpy as np
from pore2process_demo.paths import ROOT

def main():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    pointer=json.loads((ROOT/'results/adsorption/packed_bed_tracer/latest.json').read_text())
    native=(ROOT/pointer['summary']).parent
    p=ROOT/json.loads((native/'conservative_latest.json').read_text())['summary'];s=json.loads(p.read_text())
    fig,axes=plt.subplots(1,2,figsize=(11,4.3),layout='constrained',sharey=True)
    for ax,key,title in zip(axes,('bed32','length32'),('Short segment: derive D from mean only','Longer segment: same D, known Q and V')):
        r=json.loads((p.parent/f'{key}.json').read_text())
        f=np.genfromtxt(native/'native'/key/'flux.csv',delimiter=',',names=True)
        ax.plot(f['time_s'],f['outlet_flux_kg_s']*1e9,label='OpenLB pore-scale exit flux',lw=2)
        ax.plot(r['result']['time_s'],np.array(r['result']['outlet_flux_per_area_m_s'])*r['physical_fluid_area_m2']*1e9,
            'k--',label='Julia, fixed physical flow and volume',lw=1.5)
        ax.set(xlabel='Time (s)',title=title,xlim=(0,3*r['analytical_moments']['mean_s']))
        ax.grid(alpha=.2);ax.legend(loc='lower right',fontsize=8)
    axes[0].set_ylabel('Exit tracer flux (micrograms/s)')
    fig.suptitle('Pore-informed contact model: absolute flux, absorbing boundaries')
    fig.savefig(native/'conservative_comparison.png',dpi=180);plt.close(fig)

if __name__=='__main__':main()
