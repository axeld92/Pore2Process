"""Export scientific figures with measured, simulated and assumed roles explicit."""
import matplotlib
matplotlib.use('Agg')
from matplotlib.figure import Figure
import numpy as np


def plot(native,cases,rows,chosen,e,p,path):
    rows=[r for r in rows if r['viscosity_multiplier']==1.]
    fig=Figure(figsize=(13,12),layout='constrained');axes=fig.subplots(3,2)
    ax=axes[0,0];ns=[32,48,64]
    ax.plot(ns,[native['cases'][f'grid{n}']['geometry_coefficient_m4']*1e12 for n in ns],'o-',label='OpenLB D3Q27 / TRT')
    ax.axhline(native['reference_C_m4']*1e12,color='gray',ls='--',label='Independent Julia reference')
    ax.set(xlabel='Cells across outer diameter',ylabel='C (10⁻¹² m⁴)',title='Native spatial refinement');ax.legend(fontsize=8)
    ax=axes[0,1]; names=[f'high{re}_n64' for re in (1000,1388,1600,1850,2000)]
    ax.plot([native['cases'][n]['Re_mean'] for n in names],[native['cases'][n]['geometry_coefficient_m4']/native['cases']['grid64']['geometry_coefficient_m4'] for n in names],'o-')
    ax.axhspan(.98,1.02,color='#167c80',alpha=.08);ax.axhline(1,color='gray',ls='--');ax.set(ylim=(.975,1.025),xlabel='Measured native Reynolds number',ylabel='C / low-force C',title='Verified fully developed laminar branch');ax.ticklabel_format(useOffset=False,axis='y')
    ax=axes[1,0];points=e['points'];re=np.array([v['source_Re'] for v in points]);f=np.array([v['critical_flux_LMH'] for v in points]);lo=np.array([v['reading_interval_LMH'][0] for v in points]);hi=np.array([v['reading_interval_LMH'][1] for v in points])
    ax.errorbar(re,f,yerr=np.vstack((f-lo,hi-f)),fmt='o',capsize=3,label='Published biological feed + reading bounds')
    ax.plot(re,f,'--',color='gray',label='Linear scenario, not a transport law')
    ax.axvspan(*native['verified_Re_range'],alpha=.12,color='#167c80',label='Native laminar-branch envelope')
    ax.set(xlabel='Source Reynolds number',ylabel='Critical flux (LMH)',title='Chiu 2006: separate experimental evidence');ax.legend(fontsize=7)
    ax=axes[1,1];rs=p['process']['Re_targets'];js=p['process']['fluxes_LMH']
    subset=[r for r in rows if r['support']=='straight' and r['fixed_multiplier']==1.]
    energy=np.array([[next(r['specific_shaft_energy_kWh_m3'] for r in subset if r['inlet_Re']==rre and r['flux_LMH']==j) for rre in rs] for j in js])
    im=ax.imshow(energy,aspect='auto',origin='lower',cmap='viridis')
    ax.set_xticks(range(len(rs)),rs);ax.set_yticks(range(len(js)),js)
    for i in range(len(js)):
        for j in range(len(rs)):ax.text(j,i,f'{energy[i,j]:.4f}',ha='center',va='center',color='white',fontsize=9)
    ax.set(xlabel='Inlet Reynolds number',ylabel='Mean flux (LMH)',title='Actual DWSIM: straight support, fixed factor 1')
    fig.colorbar(im,ax=ax,label='Pump shaft energy (kWh/m³ filtrate)')
    ax=axes[2,0]
    for support,label in [('widened_5um','5 µm interior'),('widened_6um','6 µm interior')]:
        rr=sorted([r for r in rows if r['support']==support and r['fixed_multiplier']==1 and r['flux_LMH']==12],key=lambda r:r['inlet_Re'])
        ax.plot([r['inlet_Re'] for r in rr],[100*r['support_energy_saving_fraction'] for r in rr],'o-',label=label)
    ax.set(xlabel='Inlet Reynolds number',ylabel='Total shaft-energy saving (%)',title='Fixed 12 LMH, unchanged non-support impedance');ax.legend(fontsize=8)
    ax=axes[2,1]
    for support,label in [('straight','4 µm'),('widened_5um','5 µm'),('widened_6um','6 µm')]:
        rr=sorted([r for r in rows if r['support']==support and r['fixed_multiplier']==1 and r['flux_LMH']==15],key=lambda r:r['inlet_Re'])
        ax.plot([r['inlet_Re'] for r in rr],[r['maximum_local_flux_LMH'] for r in rr],'o-',label=label+' support: peak flux')
    queries=np.linspace(min(rs),max(rs),100)
    ax.plot(queries,.8*np.interp(queries,re,f),color='black',ls='--',label='0.8 × transferred linear Jcrit')
    ax.fill_between(queries,.8*np.interp(queries,re,lo),.8*np.interp(queries,re,hi),color='gray',alpha=.2,label='Figure-reading sensitivity')
    ax.set(xlabel='Inlet / equally mapped source Re',ylabel='Flux (LMH)',title='Local-flux screen at 15 LMH; transfer factor 1');ax.legend(fontsize=7)
    for ax in axes.flat:ax.grid(alpha=.15)
    fig.suptitle('Pore2Process: hydraulic operating map with a conditional literature capacity screen\nTRT native hydraulics → Julia scale translation → actual DWSIM recycle',fontsize=13)
    fig.savefig(path,dpi=160)
