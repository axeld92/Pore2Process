"""Plot the verified stationary-wall tracer comparisons from saved fields."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from examples.adsorption.stationary_curved_wall import reference as ref
from examples.adsorption.stationary_curved_wall.plot import history
from examples.adsorption.transport_accuracy.run import validate_record
from examples.adsorption.transport_verification.analysis import load, RMS_LIMIT
from pore2process_demo.io import sha, write_json
from pore2process_demo.paths import ROOT

from .run import HERE, OUTPUT, PARENT, parent_record

LABELS = {
    'baseline': 'Regularized TRT\nFull-way',
    'direct_fullway': 'Direct TRT\nFull-way',
    'regularized_halfway': 'Regularized TRT\nHalf-way',
    'direct_halfway': 'Direct TRT\nHalf-way',
}
COLORS = ('#42566b', '#be6337', '#007a78', '#7961a6')


def plot_fields(cases, radial, output):
    _, mask, _, _, _, coords = load(cases['baseline'])
    x, y = coords[0][:, None], coords[1][None, :]
    s, modes = radial
    fields = [ref.profile(x, y, ref.TIMES[-1], mode=(s, modes[-1]))]
    fields.extend(load(path)[3]['final'][:, :, 0] for path in cases.values())
    fields = [np.where(mask[:, :, 0] == 1, field, np.nan) for field in fields]
    span = max(np.nanmax(abs(field-.5)) for field in fields)
    fig, axes = plt.subplots(1, 5, figsize=(13, 3.2), sharex=True, sharey=True, layout='constrained')
    cmap = plt.get_cmap('coolwarm').copy()
    cmap.set_bad('#f0f2f4')
    width = ref.L*1e3/2
    for ax, label, field in zip(axes, ('Julia reference', *LABELS.values()), fields, strict=True):
        image = ax.imshow(field.T, extent=(-width, width, -width, width), origin='lower',
                          vmin=.5-span, vmax=.5+span, cmap=cmap, interpolation='nearest')
        ax.set(title=label, xlabel='x (mm)', xlim=(-.41, .41), ylim=(-.41, .41))
        ax.set_xticks([-.4, 0., .4])
    axes[0].set_ylabel('y (mm)')
    fig.suptitle(f'Varying tracer at {ref.TIMES[-1]:.3f} s')
    fig.colorbar(image, ax=list(axes), shrink=.8, label='Concentration (kg/m³)')
    fig.savefig(output, dpi=180)
    plt.close(fig)


def plot_comparison(cases, constants, output):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), layout='constrained')
    for setup, color in zip(LABELS, COLORS, strict=True):
        time, error = history(constants[setup], 0.)
        axes[0].plot(time, error, color=color, label=LABELS[setup].replace('\n', ', '))
    axes[0].set(xlabel='Time (s)', ylabel='Maximum departure from 0.5 (kg/m³)', title='Initially uniform tracer')
    axes[0].legend(fontsize=8)
    positions = np.arange(len(LABELS))
    metrics = [('middle', 'RMS_over_amplitude', 'Midpoint, all liquid'),
               ('middle', 'shell_RMS_over_amplitude', 'Midpoint, near wall'),
               ('final', 'RMS_over_amplitude', 'Final, all liquid'),
               ('final', 'shell_RMS_over_amplitude', 'Final, near wall')]
    for index, (time, metric, label) in enumerate(metrics):
        axes[1].bar(positions+(.19*(index-1.5)), [100*cases[s]['fields'][time][metric] for s in LABELS],
                    width=.18, color=('#9bcad3', '#438493', '#f4b880', '#c66529')[index], label=label)
    axes[1].axhline(100*RMS_LIMIT, color='#222222', ls='--', lw=1, label='1% screen')
    axes[1].set_xticks(positions, list(LABELS.values()), fontsize=8)
    axes[1].set(yscale='log', ylabel='RMS error (% of initial amplitude)', title='Varying tracer versus Julia')
    axes[1].legend(fontsize=8)
    fig.savefig(output, dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--summary', type=Path)
    args = parser.parse_args()
    path = (args.summary or ROOT / json.loads((OUTPUT / 'latest.json').read_text())['summary']).resolve()
    record = validate_record(path)
    parent, mode = parent_record()
    if ROOT / record['reference_path'] != mode:
        raise ValueError('Comparison used a different reference')
    cases = {'baseline': PARENT.parent / 'native/shear128'}
    constants = {'baseline': PARENT.parent / 'native/constant128'}
    varying = {'baseline': parent['cases']['shear128']}
    uniform = {'baseline': parent['cases']['constant128']}
    for setup in list(LABELS)[1:]:
        cases[setup] = path.parent / 'native' / f'{setup}_shear128'
        constants[setup] = path.parent / 'native' / f'{setup}_constant128'
        varying[setup] = record['cases'][f'{setup}_shear128']
        uniform[setup] = record['cases'][f'{setup}_constant128']
    folder = HERE / 'figures'
    folder.mkdir(exist_ok=True)
    plot_fields(cases, ref.read_mode(mode), folder / 'final_fields.png')
    plot_comparison(varying, constants, folder / 'checks.png')
    results = {
        setup: dict(uniform_maximum_departure_kg_m3=uniform[setup]['maximum_excursion_kg_m3'],
            uniform_range_kg_m3=uniform[setup]['range_kg_m3'], uniform_checks=uniform[setup]['checks'],
            varying_range_kg_m3=varying[setup]['range_kg_m3'], varying_checks=varying[setup]['checks'],
            varying_maximum_excursion_kg_m3=varying[setup]['maximum_excursion_kg_m3'],
            closed_inventory_relative_drift=max(uniform[setup]['closed_inventory_relative_drift'],
                                                varying[setup]['closed_inventory_relative_drift']),
            profile=varying[setup]['fields'])
        for setup in LABELS
    }
    write_json(HERE / 'reference_summary.json', dict(summary=str(path.relative_to(ROOT)), summary_sha256=sha(path),
        reference_path=record['reference_path'], reference_sha256=record['reference']['reference_sha256'],
        new_native_runs=record['new_native_runs'], new_bed_runs=record['new_bed_runs'],
        new_DWSIM_runs=record['new_DWSIM_runs'], passes_at_128=record['passes_at_128'], cases=results))
    print('Saved', folder)


if __name__ == '__main__':
    main()
