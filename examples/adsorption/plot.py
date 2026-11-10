import matplotlib
matplotlib.use('Agg')
from matplotlib.figure import Figure

def plot(cases,native,transport,path):
    fig=Figure(figsize=(11,8),layout='constrained');axes=fig.subplots(2,2)
    for case in cases:
        r=case['column'];label=f"{case['particle_diameter_m']*1000:g} mm particles"
        axes[0,0].plot([x/86400 for x in r['time_s']],[x/case['inputs']['inlet_concentration_kg_m3'] for x in r['outlet_concentration_kg_m3']],label=label)
    axes[0,0].axhline(.05,color='gray',ls='--',label='Chosen 5% breakthrough')
    axes[0,0].set(xlabel='Time (days)',ylabel='Outlet / inlet phenol concentration',title='Julia: conditional breakthrough');axes[0,0].legend(fontsize=8)
    sizes=[c['particle_diameter_m']*1000 for c in cases]
    axes[0,1].plot(sizes,[c['column']['treated_volume_m3'] for c in cases],'o-')
    axes[0,1].set(xlabel='Particle diameter (mm)',ylabel='Treated volume at 5% breakthrough (m³)',title='Fixed dry adsorbent inventory')
    axes[1,0].plot(sizes,[c['dwsim']['shaft_power_W']*1000 for c in cases],'o-')
    axes[1,0].set(xlabel='Particle diameter (mm)',ylabel='Actual DWSIM pump shaft power (mW)',title='Native FCC hydraulic load → actual pump')
    ns=[24,32,48]
    axes[1,1].plot(ns,[native['cases'][f'grid{n}']['permeability_m2']*1e9 for n in ns],'o-')
    axes[1,1].set(xlabel='Cells per FCC unit-cell side',ylabel='Native permeability (10⁻⁹ m²)',title='OpenLB: stock D3Q27 / TRT refinement')
    for ax in axes.flat:ax.grid(alpha=.2)
    fig.suptitle('Pore2Process adsorption: first isothermal synthetic example\nLiterature-informed uptake; experimental validation and native bed mass-transfer closures remain pending')
    fig.savefig(path,dpi=160)
