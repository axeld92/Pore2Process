"""Scientific figures from saved native fields and actual process results."""
import json
from pathlib import Path
import shutil

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from pore2process_demo.paths import ROOT
from pore2process_demo.io import sha, write_json
from .check import load_latest
from .bridge import reference_profile


HERE=Path(__file__).parent


def transport_figures():
    summary=load_latest('transport');s=json.loads(summary.read_text());base=summary.parent
    fig,ax=plt.subplots(figsize=(7.3,4.4),layout='constrained')
    c=s['cases']['outlet_step_zero32'];p=c['reference_inputs']
    f=np.genfromtxt(base/'native/outlet_step_zero32/flux.csv',delimiter=',',names=True)
    r=reference_profile(**p,t_end_s=float(f['time_s'][-1]),dt_s=.002,cells=800)
    ax.plot(r['time_s'],np.array(r['outlet_flux_per_area_m_s'])/p['velocity_m_s'],lw=2,label='Independent Julia boundary reference')
    ax.plot(f['time_s'][::80],f['outlet_flux_kg_s'][::80]/c['raw']['flow_m3_s'],'o',ms=3,mfc='none',label='OpenLB, stock zero-gradient outlet')
    ax.set(xlim=(0,25),ylim=(-.02,1.04),xlabel='Time (s)',ylabel='Exit solute flux / (Q × inlet concentration)',
        title='Convective outlet passes the uniform-flow benchmark')
    ax.grid(alpha=.2);ax.legend(fontsize=9)
    fig.savefig(base/'outlet_verification.png',dpi=180);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(11,4.5),layout='constrained')
    palette=['#1261a0','#c55315','#7063a1','#18826b']
    for model,label,color in zip(('trt27','regularized27','mrt7','trt7'),
            ('TRT D3Q27','Regularized TRT D3Q27','MRT D3Q7','TRT D3Q7'),palette):
        p=base/'native'/('application_constant_'+model)
        a=np.genfromtxt(p/'inventory.csv',delimiter=',',names=True)
        axes[0].plot(a['time_s'],a['max_concentration_kg_m3'],color=color,label=label)
        axes[0].plot(a['time_s'],a['min_concentration_kg_m3'],color=color,ls='--')
        if model=='regularized27':
            axes[1].plot(a['time_s'],a['max_concentration_kg_m3'],color=color,label='Maximum')
            axes[1].plot(a['time_s'],a['min_concentration_kg_m3'],color=color,ls='--',label='Minimum')
    for ax in axes:
        ax.axhline(1.,color='black',lw=1,label='Required uniform concentration')
        ax.set(xlabel='Time (s)',ylabel='Local concentration (kg/m³)');ax.grid(alpha=.2);ax.legend(fontsize=8)
    axes[0].set_title('Initially uniform tracer: no adsorption or imposed gradient')
    axes[1].set_title('Best tested candidate still fails preservation')
    fig.suptitle('Packed-bed transport at the declared application regime remains unqualified')
    fig.savefig(base/'transport_failure.png',dpi=180);plt.close(fig)
    for name in ('outlet_verification.png','transport_failure.png'):
        shutil.copy2(base/name,HERE/'figures'/name)
    return dict(source_summary_sha256=sha(summary),outputs={name:sha(base/name) for name in ('outlet_verification.png','transport_failure.png')})


def hydraulic_figure():
    summary=load_latest('process');s=json.loads(summary.read_text());base=summary.parent
    values=[dict(pressure_drop_Pa=s['cases'][k]['julia_output']['pressure_drop_Pa'],
        specific_pump_energy_kWh_m3=s['cases'][k]['specific_shaft_energy_Wh_m3']/1000)
        for k in ('d800','d450')]
    fig,axes=plt.subplots(1,2,figsize=(9,4),layout='constrained')
    labels=['0.80 mm','0.45 mm']
    for ax,key,unit in [(axes[0],'pressure_drop_Pa','Pressure drop (Pa)'),(axes[1],'specific_pump_energy_kWh_m3','Pump energy (Wh/m³)')]:
        y=[v[key]*(1000 if key.startswith('specific') else 1) for v in values]
        bars=ax.bar(labels,y,color=['#1261a0','#c55315'])
        ax.bar_label(bars,fmt='%.3g',padding=3)
        ax.set(ylabel=unit,ylim=(0,max(y)*1.2));ax.grid(axis='y',alpha=.2)
    fig.suptitle('Actual Julia → DWSIM duties: same flow and dry inventory\nSynthetic FCC bed; friction only')
    fig.savefig(base/'pumping_comparison.png',dpi=180);plt.close(fig)
    shutil.copy2(base/'pumping_comparison.png',HERE/'figures/pumping_comparison.png')
    return dict(source_summary_sha256=sha(summary),output_sha256=sha(base/'pumping_comparison.png'))


def main():
    write_json(ROOT/'results/adsorption/column_transport/figure_sources.json',
        dict(transport=transport_figures(),hydraulics=hydraulic_figure(),visual_inspection='pending'))


if __name__=='__main__':main()
