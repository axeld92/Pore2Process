"""Scientific source comparison and conditional process-energy figure."""
def plot(extraction, audit, cases, protocol, path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), constrained_layout=True)
    rows = extraction['garcia2008_table1']['rows']
    axes[0].scatter([r['velocity_m_s'] for r in rows], [r['critical_flux_LMH'] for r in rows],
                    marker='s', label='García table (reused critical-flux data)')
    readings = audit['figure_readings']['chiu2005_fig3']
    x = [p['x']['nominal'] for p in readings]; y = [p['critical_flux_LMH']['nominal'] for p in readings]
    dx = [[p['x']['nominal'] - p['x']['reading_interval'][0] for p in readings],
          [p['x']['reading_interval'][1] - p['x']['nominal'] for p in readings]]
    dy = [[p['critical_flux_LMH']['nominal'] - p['critical_flux_LMH']['reading_interval'][0] for p in readings],
          [p['critical_flux_LMH']['reading_interval'][1] - p['critical_flux_LMH']['nominal'] for p in readings]]
    axes[0].errorbar(x, y, xerr=dx, yerr=dy, fmt='d', capsize=3, label='Chiu figure; pixel-reading bounds')
    axes[0].axhline(protocol['target_flux_LMH'], color='black', linestyle=':', label='Target flux')
    axes[0].set(xlabel='Reported crossflow velocity (m/s)', ylabel='Critical flux (L m$^{-2}$ h$^{-1}$)', title='Source comparison; no transition interpolation')
    for convention in protocol['flow_conventions']:
        for support in sorted({c['request']['support'] for c in cases}):
            group = [c for c in cases if not c['request']['validation_only'] and c['request']['flow_convention'] == convention['id'] and c['request']['support'] == support]
            axes[1].plot([c['translation']['source_point']['velocity_m_s'] for c in group],
                [c['validation']['specific_shaft_energy_kWh_m3_reference_volume'] for c in group],
                marker='o', linestyle='--' if convention['factor'] == 7 else '-',
                label=f"{support}, flow ×{convention['factor']:g}")
    axes[1].set(xlabel='Reported crossflow velocity (m/s)', ylabel='Two-pump shaft energy (kWh/m³ permeate)',
                title='Conditional water-carrier scenarios')
    for ax in axes:
        ax.grid(alpha=.2); ax.legend(fontsize=7)
    fig.savefig(path, dpi=170)
    plt.close(fig)
