"""Scientific figures from completed Julia trajectories and digitized markers."""
import numpy as np


def figures(base,batch_records,column_records,observed):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from examples.adsorption.alhamed.study import observations as column_observations
    fig,axes=plt.subplots(2,2,figsize=(10.5,7),layout='constrained')
    by_id={r['id']:r for r in batch_records}
    for ax,(key,data) in zip(axes.flat,observed.items()):
        r=by_id[f'batch_{key}']['result'];t=np.array(r['time_s'])/60;q=np.array(r['mean_loading_kg_kg'])*1000
        ax.plot(t,q,color='tab:blue',lw=1.5,label='Shared film + particle model')
        ax.plot(data['time_s']/60,data['loading_kg_kg']*1000,'o',mfc='none',ms=4,color='black',label='Published markers')
        ax.set(xlabel='Time (min)',ylabel='Loading (mg/g dry carbon)',xlim=(0,3000),ylim=(0,52),
               title=f"{data['diameter_m']*1000:g} mm: {data['role']}")
        ax.grid(alpha=.2)
        inset=ax.inset_axes([.46,.17,.50,.42])
        inset.plot(t,q,color='tab:blue',lw=1)
        inset.plot(data['time_s']/60,data['loading_kg_kg']*1000,'o',mfc='none',ms=3,color='black')
        inset.set(xlim=(0,300),ylim=(30,50.5));inset.tick_params(labelsize=7);inset.grid(alpha=.15)
    axes[0,0].legend(fontsize=8,loc='lower left')
    fig.suptitle('Finite bath: source isotherm retained; two shared transport parameters',fontsize=12)
    fig.savefig(base/'batch_comparison.png',dpi=180);plt.close(fig)
    records={r['id']:r for r in column_records}
    fig,axes=plt.subplots(1,2,figsize=(10.8,4.3),layout='constrained',sharey=True)
    for ax,series in zip(axes,['calibration','holdout']):
        t,y=column_observations(series)
        ax.plot(t/60,y,'o',mfc='none',color='black',ms=4,label='Published markers')
        for prefix,label,style in [('primary','No column fitting','-'),('anchored','One shared film multiplier','--'),
                                   ('dispersion','Chosen positive dispersion',':')]:
            r=records[f'{prefix}_{series}']['result']
            ax.plot(np.array(r['time_s'])/60,np.array(r['outlet_concentration_kg_m3'])/.15,style,lw=1.5,label=label)
        ax.axhline(.05,color='0.6',ls=':',lw=1)
        ax.set(xlabel='Time (min)',xlim=(0,300),ylim=(0,1.02),
            title='0.8 mm: anchor only for dashed model' if series=='calibration' else '0.45 mm: withheld from all column fitting')
        ax.grid(alpha=.2);ax.legend(fontsize=7.5,loc='lower right')
    axes[0].set_ylabel('Outlet / inlet phenol concentration (1)')
    fig.suptitle('Column transfer: batch transport and source inventory frozen',fontsize=12)
    fig.savefig(base/'column_comparison.png',dpi=180);plt.close(fig)
