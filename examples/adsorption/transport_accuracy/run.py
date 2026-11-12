"""Execute native accuracy trials; no change to earlier physics or tolerances."""
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import argparse
import json
from pathlib import Path
import shutil
from uuid import uuid4

import numpy as np

from pore2process_demo.paths import ROOT
from pore2process_demo.io import sha, write_json
from pore2process_demo.provenance import source_manifest, verified_source_path
from pore2process_demo.runtime.openlb import build_case, verify_vendor
from examples.adsorption.bed import analyze as flow_analysis
from examples.adsorption.column_transport.run import execute

HERE = Path(__file__).parent
OUTPUT = ROOT / 'results/adsorption/transport_accuracy'
OLD = ROOT / 'examples/adsorption/column_transport'


def validate_record(path):
    path = Path(path).resolve()
    saved = json.loads(path.read_text())
    for name, digest in saved['source_sha256'].items():
        if sha(verified_source_path(name, digest)) != digest:
            raise ValueError('Changed parent source: ' + name)
    for name, digest in saved['artifacts_sha256'].items():
        if sha(path.parent / name) != digest:
            raise ValueError('Changed parent artifact: ' + name)
    return saved


def source_point(hydraulic, diameter, n):
    paths = sorted((hydraulic / f'native_flow/{diameter}_grid{n}').glob('iteration*/raw.json'),
                   key=lambda p: int(p.parent.name.removeprefix('iteration')))
    if not paths:
        raise ValueError('Completed native field is required.')
    native = paths[-1].parent
    flow = flow_analysis(native)
    carrier = json.loads((hydraulic / 'inputs.json').read_text())
    if abs(flow['superficial_velocity_m_s'] / carrier['native_target_velocity_m_s'] - 1) > 1e-4:
        raise ValueError('Native field does not impose the declared flow.')
    scale = carrier['carrier_viscosity_Pa_s'] / carrier['carrier_density_kg_m3'] / 1e-6
    field = np.genfromtxt(native / 'field.csv', delimiter=',', names=True)
    maxu = float(np.max(np.sqrt(sum(field[k] ** 2 for k in ('ux_m_s', 'uy_m_s', 'uz_m_s'))))) * scale
    return dict(field=str(native / 'field.csv'), velocity_scale=scale,
                side_m=flow['raw']['unit_cell_side_m'],
                speed_m_s=flow['superficial_velocity_m_s'] * scale / flow['voxel_porosity'],
                max_speed_m_s=maxu, voxel_porosity=flow['voxel_porosity'],
                field_sha256=sha(native / 'field.csv'), raw_sha256=sha(native / 'raw.json'))


def specifications(hydraulic, stage, parent=None, transits=16):
    jobs = []
    for diameter in ('d800', 'd450'):
        coarse = source_point(hydraulic, diameter, 32)
        base_dt = .02 * (coarse['side_m'] / 32) / coarse['max_speed_m_s']
        versions = [(32, 1., 'base')] if stage == 'screen' else (
            [(32, .5, 'half_dt'), (48, .5, 'grid48')] if stage == 'refine' else [(32, 1., 'long')])
        for n, dt_factor, label in versions:
            point = source_point(hydraulic, diameter, n)
            dt = base_dt * dt_factor
            for kind in (('constant', 'front') if stage != 'holdout' else ('front',)):
                periodic = kind == 'constant'
                cells = 2 if stage == 'holdout' else 1
                # Finite control-volume length excludes the two reservoir planes.
                nominal_length = (cells * 32 - (0 if periodic else 1)) * coarse['side_m'] / 32
                end = transits * nominal_length / coarse['speed_m_s']
                jobs.append(dict(id=f'{diameter}_{kind}_{label}', model=1, n=n, cells=cells,
                    D_m2_s=1e-9, tau_diffusive=.5 + 3e-9 * dt / (point['side_m'] / n) ** 2,
                    end_s=end, initial_concentration_kg_m3=1. if periodic else 0.,
                    periodic=periodic, field=point['field'], speed_m_s=point['speed_m_s'],
                    velocity_scale=point['velocity_scale'], outlet='zero', side_m=point['side_m'],
                    particle=diameter, mesh_role=label, prescribed_transits=transits,
                    source_point=point, requested_dt_s=dt))
    return jobs


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--stage', choices=['screen', 'refine', 'holdout'], default='screen')
    parser.add_argument('--hydraulic-run', type=Path)
    parser.add_argument('--parent', type=Path)
    parser.add_argument('--transits', type=int, choices=[16, 32, 64], default=16)
    parser.add_argument('--workers', type=int, choices=[1, 2], default=2)
    args = parser.parse_args()
    if args.hydraulic_run:
        hydraulic = args.hydraulic_run.resolve()
    else:
        p = ROOT / 'results/adsorption/column_transport/hydraulics/latest.json'
        hydraulic = (ROOT / json.loads(p.read_text())['summary']).parent
    parent_hydraulic = validate_record(hydraulic / 'summary.json')
    if not parent_hydraulic['hydraulic_checks_passed']:
        raise ValueError('Qualified hydraulics are required.')
    if args.stage != 'screen' and args.parent is None:
        raise ValueError('Refinement/holdout requires an explicit screening parent directory.')
    parent = None
    if args.parent:
        parent = validate_record(args.parent / 'summary.json')
        if parent['hydraulic_summary_sha256'] != sha(hydraulic / 'summary.json'):
            raise ValueError('The parent uses different hydraulic evidence.')
    jobs = specifications(hydraulic, args.stage, parent, args.transits)
    rid = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ') + '-' + uuid4().hex[:8]
    base = OUTPUT / 'runs' / rid
    base.mkdir(parents=True)
    files = [Path(__file__), HERE / 'protocol.md', OLD / 'run.py',
             ROOT / 'examples/adsorption/bed.py', *list((OLD / 'openlb').iterdir()),
             ROOT / 'uv.lock', ROOT / 'julia/Project.toml', ROOT / 'julia/Manifest.toml']
    sources = source_manifest(files)
    for name in sources:
        target = base / 'source_snapshot' / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT / name, target)
    generated = base / 'generated_source'
    shutil.copytree(OLD / 'openlb', generated)
    config = generated / 'config.mk'
    config.write_text(config.read_text() + '\nCXXFLAGS += -DP2P_SCALAR_MODEL=1\n')
    build = ROOT / '.local/adsorption/transport_accuracy'
    print('Native accuracy run', base, flush=True)
    metadata = build_case(generated, build, base / 'build')
    binary = build / 'capillary3d'
    record = dict(schema='pore2process.transport-accuracy-native.v1', run_id=rid,
        stage=args.stage, transits=args.transits, source_sha256=sources,
        hydraulic_summary=str((hydraulic / 'summary.json').relative_to(ROOT)),
        hydraulic_summary_sha256=sha(hydraulic / 'summary.json'),
        parent_summary=None if args.parent is None else str((args.parent.resolve() / 'summary.json').relative_to(ROOT)),
        parent_summary_sha256=None if args.parent is None else sha(args.parent / 'summary.json'),
        specifications=jobs, build=metadata, new_native_runs=len(jobs), new_DWSIM_runs=0,
        new_adsorption_predictions=0, physics_changed=False)
    write_json(base / 'plan.json', record)
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        list(pool.map(lambda spec: execute(base, spec, binary), jobs))
    record['vendor'] = verify_vendor()
    record['artifacts_sha256'] = {str(p.relative_to(base)): sha(p) for p in base.rglob('*') if p.is_file()}
    write_json(base / 'summary.json', record)
    write_json(OUTPUT / (args.stage + '_latest.json'), dict(summary=str((base / 'summary.json').relative_to(ROOT))))
    print('Native stage completed', base, 'actual cases', len(jobs), flush=True)


if __name__ == '__main__':
    main()
