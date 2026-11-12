"""Diagnostic of the installed stock regularized-temperature/TRT inlet."""
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import argparse
import json
from pathlib import Path
import shutil
import subprocess
from time import perf_counter
from uuid import uuid4

from pore2process_demo.paths import ROOT
from pore2process_demo.io import sha, write_json
from pore2process_demo.provenance import source_manifest
from pore2process_demo.runtime.openlb import build_case, verify_vendor
from .run import HERE, OUTPUT, source_point, validate_record


def execute(base, spec, binary):
    out = base / 'native' / spec['id']
    out.mkdir(parents=True)
    command = [str(binary), str(spec['n']), str(spec['cells']), str(spec['side_m']),
        str(spec['D_m2_s']), str(spec['tau_diffusive']), str(spec['end_s']), spec['field'],
        str(spec['speed_m_s']), '0.', str(spec['velocity_scale']), '0', 'zero', str(spec['ramp_s'])]
    write_json(out / 'request.json', dict(spec=spec, command=command,
        binary_sha256=sha(binary), native_field_sha256=sha(Path(spec['field']))))
    print('Native smooth-input start', spec['id'], flush=True)
    start = perf_counter()
    with (out / 'stdout.log').open('w') as stream:
        result = subprocess.run(command, cwd=out, stdout=stream, stderr=subprocess.STDOUT, timeout=1800)
    write_json(out / 'execution.json', dict(returncode=result.returncode, elapsed_s=perf_counter()-start))
    if result.returncode:
        raise RuntimeError('Native run failed: ' + str(out / 'stdout.log'))
    print('Native smooth-input finish', spec['id'], f'{perf_counter()-start:.1f}s', flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--hydraulic-run', type=Path)
    parser.add_argument('--n', type=int, choices=[32, 48], default=32)
    parser.add_argument('--half-dt', action='store_true')
    parser.add_argument('--cells', type=int, choices=[1, 2], default=1)
    parser.add_argument('--transits', type=int, choices=[2, 16, 32, 64], default=2)
    parser.add_argument('--ramp-transits', type=float, nargs='+', choices=[.5, 1.], default=[1.])
    args = parser.parse_args()
    if args.hydraulic_run:
        hydraulic = args.hydraulic_run.resolve()
    else:
        pointer = ROOT / 'results/adsorption/column_transport/hydraulics/latest.json'
        hydraulic = (ROOT / json.loads(pointer.read_text())['summary']).parent
    parent = validate_record(hydraulic / 'summary.json')
    if not parent['hydraulic_checks_passed']:
        raise ValueError('Qualified hydraulics required.')
    rid = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ') + '-' + uuid4().hex[:8]
    base = OUTPUT / 'regularized_inlet/runs' / rid
    base.mkdir(parents=True)
    sources = source_manifest([Path(__file__), HERE / 'run.py', HERE / 'protocol.md',
        HERE / 'ramp_protocol.md', HERE / 'inlet_protocol.md', *list((HERE / 'openlb_regularized_inlet').iterdir()),
        ROOT / 'examples/adsorption/bed.py'])
    for name in sources:
        dest = base / 'source_snapshot' / name
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT / name, dest)
    generated = base / 'generated_source'
    shutil.copytree(HERE / 'openlb_regularized_inlet', generated)
    config = generated / 'config.mk'
    config.write_text(config.read_text() + '\nCXXFLAGS += -DP2P_SCALAR_MODEL=1\n')
    build = ROOT / '.local/adsorption/transport_accuracy_regularized_inlet'
    print('Smooth-input run', base, flush=True)
    metadata = build_case(generated, build, base / 'build')
    jobs = []
    for diameter in ('d800', 'd450'):
        coarse = source_point(hydraulic, diameter, 32)
        point = source_point(hydraulic, diameter, args.n)
        dt = .02 * (coarse['side_m'] / 32) / coarse['max_speed_m_s']
        if args.half_dt:
            dt /= 2
        transit = coarse['side_m'] / coarse['speed_m_s']
        control_transit = ((args.cells * 32 - 1) / 32) * transit
        for ramp in args.ramp_transits:
            jobs.append(dict(id=f'{diameter}_ramp{ramp:g}_n{args.n}' + ('_halfdt' if args.half_dt else '') + f'_cells{args.cells}',
                model=1, n=args.n, cells=args.cells, D_m2_s=1e-9,
                tau_diffusive=.5 + 3e-9 * dt / (point['side_m'] / args.n) ** 2,
                end_s=ramp * transit + args.transits * control_transit,
                initial_concentration_kg_m3=0., periodic=False,
                field=point['field'], speed_m_s=point['speed_m_s'], velocity_scale=point['velocity_scale'],
                outlet='zero', side_m=point['side_m'], particle=diameter, source_point=point,
                ramp_s=ramp*transit, ramp_transits=ramp, prescribed_post_ramp_transits=args.transits,
                requested_dt_s=dt))
    record = dict(schema='pore2process.transport-accuracy-native.v1', run_id=rid,
        stage='regularized_inlet', source_sha256=sources,
        hydraulic_summary=str((hydraulic / 'summary.json').relative_to(ROOT)),
        hydraulic_summary_sha256=sha(hydraulic / 'summary.json'), specifications=jobs, build=metadata,
        new_native_runs=len(jobs), new_DWSIM_runs=0, new_adsorption_predictions=0,
        operators_changed=True, change='Stock regularized-temperature/TRT inlet; other physical inputs and bulk/wall/outlet operators unchanged.')
    write_json(base / 'plan.json', record)
    with ThreadPoolExecutor(max_workers=2) as pool:
        list(pool.map(lambda spec: execute(base, spec, build / 'capillary3d'), jobs))
    record['vendor'] = verify_vendor()
    record['artifacts_sha256'] = {str(p.relative_to(base)): sha(p) for p in base.rglob('*') if p.is_file()}
    write_json(base / 'summary.json', record)
    write_json(OUTPUT / 'regularized_inlet_latest.json', dict(summary=str((base / 'summary.json').relative_to(ROOT))))
    print('Smooth-input trials completed', base, flush=True)


if __name__ == '__main__':
    main()
