"""Scientific comparison from actual DWSIM stream results (reference volumes)."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

def plot(cases,path,recovery):
    fig,axes=plt.subplots(1,3,figsize=(13,4.2),layout='constrained')
    names=['straight','widened_5um','widened_6um'];labels=['4 µm','5 µm','6 µm']
    colors=['#184e77','#168aad']
    for channel,color in zip((0.,20000.),colors,strict=True):
        group=[c for c in cases if c['recovery_fraction']==recovery and c['channel_pressure_drop_Pa']==channel]
        x=np.arange(3)+(-.18 if channel==0 else .18)
        energies=[c['validation']['actual_specific_shaft_energy_kWh_m3_reference_volume'] for c in group]
        axes[0].bar(x,energies,width=.34,label=f'{channel/1000:g} kPa channel loss',color=color)
        for c,e in zip(group,energies,strict=True):
            axes[1].plot(c['julia_membrane']['measured_body_resistance_ratio'],100*c['energy_saving_vs_straight_fraction'],'o',color=color)
        for name,ls in zip(names,('-','--',':'),strict=True):
            g=[c for c in cases if c['support']==name and c['channel_pressure_drop_Pa']==channel]
            axes[2].plot([100*c['recovery_fraction'] for c in g],
                [c['validation']['actual_specific_shaft_energy_kWh_m3_reference_volume'] for c in g],
                ls+'o',color=color,label=f'{labels[names.index(name)]}, {channel/1000:g} kPa')
    axes[0].set_xticks(np.arange(3),labels);axes[0].set_xlabel('Synthetic pore interior radius (fixed 4 µm mouths)')
    axes[0].set_ylabel('Pump shaft energy (kWh/m³ permeate at feed density)')
    axes[0].set_title(f'{100*recovery:g}% recovery');axes[0].legend(fontsize=8)
    axes[1].set_xlabel('Native OpenLB body resistance / straight body')
    axes[1].set_ylabel('Pump shaft energy saving (%)');axes[1].set_title('Separate losses dilute material benefit')
    axes[2].set_xlabel('Permeate/feed recovery (%)');axes[2].set_ylabel('Pump shaft energy (kWh/m³ permeate)')
    axes[2].set_title('All fresh feed is pressurized');axes[2].legend(fontsize=7)
    for ax in axes: ax.grid(alpha=.2,axis='y');ax.set_axisbelow(True)
    fig.suptitle('Once-through pure-water scenario • 100 LMH, 1 m² • no recycle or pressure recovery',fontsize=11)
    fig.savefig(path,dpi=180);plt.close(fig)
