"""Scientific plots of verified fields and actual process results."""
import matplotlib
matplotlib.use('Agg')
from matplotlib.figure import Figure
import numpy as np
from .geometry import polygon

def profile_figure(t):
    a=t['axial'];f=t['feed'];fig=Figure(figsize=(9,3.5),layout='constrained')
    ax,bx=fig.subplots(1,2)
    x=np.asarray(a['x_m'])*1e3
    ax.plot(x,np.asarray(a['pressure_Pa'])/1e3,label='Local TMP')
    ax.set(xlabel='Axial position (mm)',ylabel='TMP (kPa)');ax.grid(alpha=.2)
    by=bx.twinx();bx.plot(x,np.asarray(a['local_flux_m_s'])*3.6e6,color='#167c80',label='Flux')
    by.plot(x,np.asarray(f['axial_COD_kg_m3'])*1000,color='#ae5743',linestyle='--',label='COD')
    if t['critical_flux_limit_LMH'] is not None:
        bx.axhline(t['critical_flux_limit_LMH'],color='gray',linestyle=':',label='Transferred limit')
    bx.set(xlabel='Axial position (mm)',ylabel='Local flux (LMH)');by.set_ylabel('COD (mg/L)')
    lines,labels=bx.get_legend_handles_labels();other_lines,other_labels=by.get_legend_handles_labels()
    bx.legend(lines+other_lines,labels+other_labels,loc='lower left',fontsize=8);bx.grid(alpha=.2)
    fig.suptitle(f"Clean initial hydraulic case · inlet Re={a['inlet_Re']:g} · illustrative COD retention")
    return fig

def plot_study(native,cases,rows,path):
    fig=Figure(figsize=(11,8),layout='constrained');axes=fig.subplots(2,2)
    v=polygon()*1e3;v=np.vstack((v,v[0]));ax=axes[0,0]
    ax.fill(v[:,0],v[:,1],color='#dbedf0');ax.plot(v[:,0],v[:,1],color='#176d7c')
    ax.set(xlabel='y relative to centre (mm)',ylabel='z relative to centre (mm)',title='Declared six-point synthetic star',aspect='equal')
    ax=axes[0,1];ns=native['protocol']['native']['resolutions'];cs=[native['cases'][f'grid{n}']['geometry_coefficient_m4']*1e12 for n in ns]
    ax.plot(ns,cs,'o-',label='OpenLB')
    ax.axhline(native['references'][-1]['conductance_geometry_m4']*1e12,linestyle='--',color='gray',label='Independent Julia reference')
    ax.set(xlabel='Cells across outer diameter',ylabel='Geometry coefficient (10⁻¹² m⁴)',title='Hydraulic verification');ax.legend(fontsize=8)
    ax=axes[1,0]
    nominal=max(r['flux_LMH'] for r in rows)
    subset=[r for r in rows if r['flux_LMH']==nominal and r['rejection']==.65 and r['viscosity_multiplier']==r['fixed_multiplier']==1 and r['inlet_Re']==20]
    for i,r in enumerate(subset):
        scale=1/(nominal/3.6e6)/3.6e6
        ax.bar(i,r['feed_power_W']*scale,color='#176d7c',label='Fresh-feed pump' if i==0 else None)
        ax.bar(i,r['circulation_power_W']*scale,bottom=r['feed_power_W']*scale,color='#dd9c58',label='Circulation pump' if i==0 else None)
    ax.set_xticks(range(len(subset)),['4 µm','5 µm','6 µm']);ax.set(xlabel='Support interior radius; fixed 4 µm openings',ylabel='Pump shaft energy (kWh/m³)',title=f'{nominal:g} LMH biological-feed scenario');ax.legend(fontsize=8)
    ax=axes[1,1];c=next(c for c in cases if c['request']['support']=='straight' and c['request']['flux_LMH']==nominal)
    a=c['translation']['axial'];ax.plot(np.asarray(a['x_m'])*1e3,np.asarray(a['local_flux_m_s'])*3.6e6,label='Local flux')
    if c['translation']['critical_flux_limit_LMH'] is not None:
        ax.axhline(c['translation']['critical_flux_limit_LMH'],color='#ae5743',linestyle='--',label='Transferred conservative limit')
    ax.set(xlabel='Axial position (mm)',ylabel='Flux (LMH)',title='Distributed membrane calculation');ax.legend(fontsize=8)
    for ax in axes.flat:ax.grid(alpha=.2)
    fig.suptitle('OpenLB → Julia → DWSIM · verified low-Re hydraulics, illustrative biological-feed retention')
    fig.savefig(path,dpi=160)
