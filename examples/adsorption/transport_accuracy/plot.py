"""Scientific plots of saved native histories; no new simulation or UI."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from pore2process_demo.paths import ROOT
from pore2process_demo.io import sha, write_json
from .run import HERE, OUTPUT, validate_record


def stage(name):
    pointer = json.loads((OUTPUT / (name + '_latest.json')).read_text())
    summary = ROOT / pointer['summary']
    return summary.parent, validate_record(summary)


def main():
    figures = HERE / 'figures'
    figures.mkdir(exist_ok=True)
    sources = {}
    base, native = stage('screen')
    sources['screen'] = dict(path=str((base / 'summary.json').relative_to(ROOT)), sha256=sha(base / 'summary.json'))
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.0), layout='constrained')
    for diameter, label, colour in [('d800', '0.80 mm', '#2166ac'), ('d450', '0.45 mm', '#b35806')]:
        p = base / 'native' / (diameter + '_constant_base')
        r = json.loads((p / 'raw.json').read_text())
        a = np.genfromtxt(p / 'inventory.csv', delimiter=',', names=True)
        t = a['time_s'] / (r['fluid_volume_m3'] / r['flow_m3_s'])
        axes[0].plot(t, 100*(a['min_concentration_kg_m3']-1), color=colour, label=label+' minimum')
        axes[0].plot(t, 100*(a['max_concentration_kg_m3']-1), color=colour, ls='--', label=label+' maximum')
        axes[1].plot(t, 100*(a['fluid_tracer_kg']/a['fluid_tracer_kg'][0]-1), color=colour, label=label)
    for ax in axes:
        ax.axhline(0, color='.5', lw=.7)
        ax.set_xlabel('Time / physical void time (–)')
        ax.grid(alpha=.2)
        ax.legend(fontsize=8)
    axes[0].set_ylabel('Local concentration departure (%)')
    axes[1].set_ylabel('Average liquid concentration change (%)')
    axes[0].set_title('Initially uniform inert tracer')
    axes[1].set_title('Physical liquid only; wall storage excluded')
    fig.savefig(figures / 'uniform_accuracy.png', dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(10.4, 7.0), layout='constrained')
    variants = [('screen', 'Sharp inlet', '-', '#b2182b'),
                ('ramp', 'Gradual inlet: 1 transit', '--', '#2166ac'),
                ('regularized_inlet', 'Gradual + regularized inlet', '-.', '#1b7837')]
    for folder, label, style, colour in variants:
        base, s = stage(folder)
        sources[folder] = dict(path=str((base / 'summary.json').relative_to(ROOT)), sha256=sha(base / 'summary.json'))
        for column, diameter in enumerate(('d800', 'd450')):
            candidates = [p for p in s['specifications'] if p['particle'] == diameter and not p['periodic'] and
                          (folder == 'screen' or p['ramp_transits'] == 1.)]
            spec = candidates[0]
            p = base / 'native' / spec['id']
            r = json.loads((p / 'raw.json').read_text())
            a = np.genfromtxt(p / 'inventory.csv', delimiter=',', names=True)
            f = np.genfromtxt(p / 'flux.csv', delimiter=',', names=True)
            tau = r['fluid_volume_m3'] / r['flow_m3_s']
            axes[0, column].plot(a['time_s']/tau, a['min_concentration_kg_m3'], color=colour, ls=style, label=label+' min')
            axes[0, column].plot(a['time_s']/tau, a['max_concentration_kg_m3'], color=colour, ls=style, alpha=.55)
            axes[1, column].plot(f['time_s']/tau, f['outlet_flux_kg_s']/r['flow_m3_s'], color=colour, ls=style, label=label)
    for column, diameter in enumerate(('0.80 mm', '0.45 mm')):
        axes[0, column].set_title(diameter+' particles')
        for row in (0, 1):
            ax = axes[row, column]
            ax.axhline(0, color='.5', lw=.8)
            ax.axhline(1, color='.5', lw=.8)
            ax.set_xlim(0, 4)
            ax.set_xlabel('Time / physical void time (–)')
            ax.grid(alpha=.2)
            ax.legend(fontsize=7, loc='lower right')
    axes[0, 0].set_ylabel('Local C / final feed C (–)')
    axes[1, 0].set_ylabel('Outlet solute flux / (water Q × feed C) (–)')
    fig.suptitle('Native front diagnostics; later tails retained in saved data')
    fig.savefig(figures / 'front_diagnostics.png', dpi=180)
    plt.close(fig)

    base, s = stage('profile')
    sources['profile'] = dict(path=str((base / 'summary.json').relative_to(ROOT)), sha256=sha(base / 'summary.json'))
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.0), layout='constrained')
    p = base / 'native/d800_free'
    r = json.loads((p / 'raw.json').read_text())
    f = np.genfromtxt(p / 'final.csv', delimiter=',', names=True)
    shape = (r['nx'], r['n'], r['n'])
    x = f['x_m'].reshape(shape)[:, 0, 0]
    values = f['concentration_kg_m3'].reshape(shape).mean(axis=(1, 2))
    t = r['steps'] * r['dt_s']
    u = f['ux_m_s'][0]
    exact = .5+.25*np.exp(-r['diffusivity_m2_s']*(2*np.pi/r['length_m'])**2*t)*np.sin(
        2*np.pi*(x-u*t)/r['length_m'])
    axes[0].plot(x*1000, exact, color='.25', label='Analytical solution')
    axes[0].plot(x*1000, values, 'o', ms=3.5, color='#2166ac', label='Actual OpenLB, 0.80 mm scale')
    axes[0].set_xlabel('Position x (mm)')
    axes[0].set_ylabel('Concentration (kg/m³)')
    axes[0].set_title('Smooth periodic profile in free flow')
    axes[0].legend(fontsize=8)
    for diameter, label, colour in [('d800', '0.80 mm', '#2166ac'), ('d450', '0.45 mm', '#b35806')]:
        p = base / 'native' / (diameter + '_porous')
        r = json.loads((p / 'raw.json').read_text())
        a = np.genfromtxt(p / 'inventory.csv', delimiter=',', names=True)
        t = a['time_s']/(r['fluid_volume_m3']/r['flow_m3_s'])
        axes[1].plot(t, a['min_concentration_kg_m3'], color=colour, label=label+' minimum')
        axes[1].plot(t, a['max_concentration_kg_m3'], color=colour, ls='--', label=label+' maximum')
    for level in (.25, .75):
        axes[1].axhline(level, color='.4', lw=1., ls=':')
    axes[1].set_xlabel('Time / physical void time (–)')
    axes[1].set_ylabel('Concentration (kg/m³)')
    axes[1].set_title('Same profile in the periodic porous bed')
    axes[1].legend(fontsize=8)
    for ax in axes:
        ax.grid(alpha=.2)
    fig.savefig(figures / 'periodic_isolation.png', dpi=180)
    plt.close(fig)
    write_json(OUTPUT / 'figure_sources.json', dict(sources=sources,
        plotting_source_sha256=sha(Path(__file__)),
        files={p.name: sha(p) for p in figures.glob('*.png')}, visual_inspection='pending'))
    print('Saved figures', figures)


if __name__ == '__main__':
    main()
