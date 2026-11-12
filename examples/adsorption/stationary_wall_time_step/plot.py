"""Summarize the saved time-step comparison and plot tracer errors."""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from examples.adsorption.stationary_curved_wall.plot import history
from examples.adsorption.stationary_curved_wall import reference as ref
from examples.adsorption.transport_accuracy.run import validate_record
from examples.adsorption.transport_verification.analysis import EXCURSION_LIMIT, RMS_LIMIT, load
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT

from .run import CONTROL, HERE, NAMES, OUTPUT, control_record

COLORS = ('#466886', '#008078')
LABELS = ('Original step', 'Half step')


def varying_profile(axes, original, half):
    positions = np.arange(4)
    for offset, cases, label, color in zip((-.17, .17), (original, half), LABELS, COLORS, strict=True):
        values = [100*cases['fields'][time][region] for time in ('middle', 'final')
                  for region in ('RMS_over_amplitude', 'shell_RMS_over_amplitude')]
        axes.bar(positions+offset, values, width=.31, label=label, color=color)
    axes.axhline(100*RMS_LIMIT, color='#323232', ls='--', lw=1, label='1% target')
    axes.set_xticks(positions, ['Midpoint\nAll liquid', 'Midpoint\nNear wall',
                                'Final\nAll liquid', 'Final\nNear wall'], fontsize=8)
    axes.set(title='Varying tracer vs Julia', ylabel='RMS error (% of initial amplitude)')
    axes.legend(fontsize=8)


def comparison(path, summary, control):
    base = path.parent / 'native'
    old = CONTROL.parent / 'native'
    fig, axes = plt.subplots(1, 3, figsize=(13.2, 3.7), layout='constrained')
    uniform = axes[0]
    uniform.axhline(.5, color='#323232', ls='--', lw=1, label='Exact')
    for directory, label, color in zip((old / NAMES['constant'], base / 'constant'), LABELS, COLORS, strict=True):
        _, _, _, _, h, _ = load(directory)
        uniform.fill_between(h['time_s'], h['min_concentration_kg_m3'], h['max_concentration_kg_m3'],
                             color=color, alpha=.17)
        uniform.plot(h['time_s'], h['max_concentration_kg_m3'], color=color, label=label)
        uniform.plot(h['time_s'], h['min_concentration_kg_m3'], color=color)
    uniform.set(title='Initially uniform tracer', xlabel='Time (s)', ylabel='Liquid concentration range (kg/m³)')
    uniform.legend(fontsize=8)
    bounds = axes[1]
    for directory, label, color in zip((old / NAMES['varying'], base / 'varying'), LABELS, COLORS, strict=True):
        time, excursion = history(directory, ref.A)
        bounds.plot(time, excursion, color=color, label=label)
    bounds.axhline(EXCURSION_LIMIT, color='#323232', ls='--', lw=1, label='Allowance')
    bounds.set(title='Varying tracer bounds', xlabel='Time (s)', ylabel='Excursion (kg/m³)')
    bounds.legend(fontsize=8)
    varying_profile(axes[2], control['cases'][NAMES['varying']], summary['analysis']['varying'])
    out = HERE / 'figures'
    out.mkdir(exist_ok=True)
    fig.savefig(out / 'comparison.png', dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--summary', type=Path)
    args = parser.parse_args()
    path = (args.summary or ROOT / read_json(OUTPUT / 'latest.json')['summary']).resolve()
    summary = validate_record(path)
    control, mode, _ = control_record()
    if ROOT / summary['reference_path'] != mode or summary['control_summary_sha256'] != sha(CONTROL):
        raise ValueError('Comparison used a different reference or control')
    comparison(path, summary, control)
    write_json(HERE / 'reference_summary.json', dict(summary=str(path.relative_to(ROOT)), summary_sha256=sha(path),
        control_summary=summary['control_summary'], control_summary_sha256=summary['control_summary_sha256'],
        reference_sha256=summary['reference']['reference_sha256'],
        original={name: control['cases'][original] for name, original in NAMES.items()},
        half_step=summary['analysis'], new_native_runs=summary['new_native_runs'],
        new_bed_runs=summary['new_bed_runs'], new_DWSIM_runs=summary['new_DWSIM_runs']))
    print('Saved', HERE / 'figures/comparison.png')


if __name__ == '__main__':
    main()
