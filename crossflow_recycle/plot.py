"""Static scientific figures from actual DWSIM recycle results."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

def plot(cases,path):
    fig,axes=plt.subplots(1,3,figsize=(13,4.4),layout='constrained')
    selected=[v for v in cases if not v['request']['validation_only']]
    colors=['#184e77','#168aad','#58b6a0'];names=['straight','widened_5um','widened_6um']
    for name,color in zip(names,colors,strict=True):
        g=[v for v in selected if v['request']['support']==name and (v['request']['condition']=='nominal' or v['request']['condition'].startswith('circulation_'))]
        x=[v['julia_process']['module_to_permeate_ratio'] for v in g]
        axes[0].plot(x,[v['validation']['specific_shaft_energy_kWh_m3_reference_volume'] for v in g],'-o',color=color,label=f'{4+names.index(name)} µm interior')
        axes[1].plot(x,[100*v['energy_saving_vs_straight_fraction'] for v in g],'-o',color=color)
    g=[v for v in selected if v['request']['support']=='straight' and (v['request']['condition']=='nominal' or v['request']['condition'].startswith('circulation_'))]
    x=np.arange(len(g));feed=[v['dwsim']['powers']['feed_W'] for v in g];circ=[v['dwsim']['powers']['circulation_W'] for v in g]
    axes[2].bar(x,feed,color=colors[0],label='Fresh-feed pump')
    axes[2].bar(x,circ,bottom=feed,color=colors[1],label='Circulation pump')
    axes[2].plot(x,[v['dwsim']['powers']['cooling_removed_W'] for v in g],'D--',color='#b35c1e',label='Cooling heat removed')
    axes[2].set_xticks(x,[f'{v["julia_process"]["module_to_permeate_ratio"]:g}' for v in g])
    axes[0].set_ylabel('Shaft energy (kWh/m³ permeate at feed density)');axes[0].legend(fontsize=8)
    axes[1].set_ylabel('Total shaft energy saving vs straight (%)')
    axes[2].set_ylabel('Power / heat duty (W)');axes[2].legend(fontsize=8)
    axes[0].set_title('Actual DWSIM pump duties');axes[1].set_title('Same pore benefit, larger circulation load')
    axes[2].set_title('Straight support: separate work and heat')
    for ax in axes:ax.set_xlabel('Module flow / permeate flow');ax.grid(axis='y',alpha=.2);ax.set_axisbelow(True)
    fig.suptitle('Clean-water recycle scenario • 100 LMH, 1 m², 50% overall recovery • fixed channel resistance',fontsize=11)
    fig.savefig(path,dpi=180);plt.close(fig)
