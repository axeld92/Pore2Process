"""Plot the working-mesh disk solution and its error at the archived OpenLB points."""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from examples.adsorption.stationary_curved_wall import reference as ref
from examples.adsorption.transport_accuracy.run import validate_record
from examples.adsorption.transport_verification.analysis import load
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT

from .run import CHECK_TIMES, HERE, OUTPUT, PARENT, sample, spatial_grid


def fields(path, record):
    chosen = record['analysis']['selected']
    case = record['analysis']['cases'][chosen]['settings']
    edges, theta = spatial_grid(case)
    mode_path = path.parent / 'reference.csv'
    if sha(mode_path) != record['reference']['reference_sha256']:
        raise ValueError('Wrong archived radial reference')
    r, modes = ref.read_mode(mode_path)
    _, mask, _, _, _, coords = load(PARENT.parent / 'native/shear128')
    fluid = mask[:, :, 0] == 1
    x, y = np.meshgrid(coords[0]-ref.L/2, coords[1]-ref.L/2, indexing='ij')
    result = []
    for label, time, radial in zip(CHECK_TIMES, ref.TIMES, modes, strict=True):
        target = ref.profile(x[fluid]+ref.L/2, y[fluid]+ref.L/2, time, mode=(r,radial))
        native = np.fromfile(path.parent / 'cases' / chosen / (label+'.f64'), dtype='<f8')
        native = native.reshape((case['nr'],case['na']), order='F')
        sampled = sample(native, edges, theta, x[fluid], y[fluid])
        row = []
        for a in (target, sampled, sampled-target):
            field = np.full(fluid.shape, np.nan)
            field[fluid] = a
            row.append(field)
        result.append(row)
    return result


def plot(data):
    fig, axes = plt.subplots(2, 3, figsize=(9.9, 6.5), layout='constrained', sharex=True, sharey=True)
    cmap = plt.get_cmap('coolwarm').copy()
    cmap.set_bad('#eaeef1')
    width = ref.L*1e3/2
    error_limit = max(np.nanmax(abs(row[2])) for row in data)
    for i, row in enumerate(data):
        for j, field in enumerate(row):
            center, span = (0., error_limit) if j == 2 else (.5, ref.A)
            im = axes[i,j].imshow(field.T, extent=(-width,width,-width,width), origin='lower',
                                 cmap=cmap, vmin=center-span, vmax=center+span, interpolation='nearest')
            axes[i,j].set(xlim=(-.415,.415), ylim=(-.415,.415), xticks=[-.4,0.,.4], yticks=[-.4,0.,.4])
            if i == 0:
                axes[i,j].set_title(('Radial reference', 'Balanced face flows', 'Difference')[j])
                if j in (1,2):
                    fig.colorbar(im, ax=list(axes[:, :2].flat) if j == 1 else list(axes[:,2]),
                                 label=('Concentration (kg/m³)' if j == 1 else 'Difference (kg/m³)'), shrink=.72)
            if i == 1:
                axes[i,j].set_xlabel('x (mm)')
        axes[i,0].set_ylabel(f'{ref.TIMES[i]:.3f} s\ny (mm)')
    out = HERE / 'figures'
    out.mkdir(exist_ok=True)
    fig.savefig(out / 'profiles.png', dpi=175)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--summary', type=Path)
    parser.add_argument('--confirmation', type=Path)
    args = parser.parse_args()
    path = (args.summary or ROOT / read_json(OUTPUT / 'latest.json')['summary']).resolve()
    confirm = (args.confirmation or ROOT / read_json(OUTPUT / 'latest_confirmation.json')['summary']).resolve()
    initial, finer = validate_record(path), validate_record(confirm)
    if finer['initial_summary_sha256'] != sha(path):
        raise ValueError('Confirmation used a different working study')
    plot(fields(path, initial))
    write_json(HERE / 'reference_summary.json', dict(
        initial_summary=str(path.relative_to(ROOT)), initial_summary_sha256=sha(path),
        confirmation_summary=str(confirm.relative_to(ROOT)), confirmation_summary_sha256=sha(confirm),
        radial_reference_sha256=initial['reference']['reference_sha256'],
        initial=initial['analysis'], confirmation=finer['analysis'],
        new_openlb_runs=initial['new_openlb_runs']+finer['new_openlb_runs'],
        new_bed_runs=initial['new_bed_runs']+finer['new_bed_runs'],
        new_DWSIM_runs=initial['new_DWSIM_runs']+finer['new_DWSIM_runs']))
    print('Saved', HERE / 'figures/profiles.png')


if __name__ == '__main__':
    main()
