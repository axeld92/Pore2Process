"""Run or replay three stock alternatives on the archived stationary-disk benchmark."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import shutil
import subprocess
from time import perf_counter
from uuid import uuid4

from examples.adsorption.stationary_curved_wall import reference as ref
from examples.adsorption.transport_accuracy.run import validate_record
from examples.adsorption.transport_root_cause.hydraulics.stock_omp import build_case_omp
from pore2process_demo.io import sha, write_json
from pore2process_demo.paths import ROOT
from pore2process_demo.provenance import source_manifest
from pore2process_demo.runtime.openlb import verify_vendor

from .analysis import SETUPS, analyze, passes

HERE = Path(__file__).parent
OUTPUT = ROOT / 'results/adsorption/stock_transport_candidate'
PARENT = ROOT / 'results/adsorption/stationary_curved_wall/runs/20261111T122659Z-56628d90/summary.json'
PARENT_SHA = '4c63619355b7072e6263c3c9f72d7b063754181e421392fb1f3a40b06322ecb6'
THREADS = 2
WORKERS = 2


def parent_record():
    if sha(PARENT) != PARENT_SHA:
        raise ValueError('Archived stationary-wall summary differs from the selected benchmark')
    saved = validate_record(PARENT)
    qualification = saved['reference']
    mode = PARENT.parent / 'reference' / qualification['reference_file']
    if not qualification['passed'] or sha(mode) != qualification['reference_sha256']:
        raise ValueError('Archived Julia reference differs or failed its convergence checks')
    return saved, mode


def specifications(n, setups=SETUPS):
    return [dict(setup=setup, n=n, amplitude=amplitude, name=f'{setup}_{kind}{n}')
            for setup in setups for kind, amplitude in [('constant', 0.), ('shear', ref.A)]]


def parent_case(spec):
    return ('shear' if spec['amplitude'] else 'constant') + str(spec['n'])


def analyze_case(directory, spec, radial):
    return analyze(directory, PARENT.parent / 'native' / parent_case(spec), radial)


def execute(base, binary, spec, radial):
    directory = base / 'native' / spec['name']
    directory.mkdir(parents=True)
    args = ['stationary_disk', spec['n'], 4, ref.STEPS, ref.L, ref.DT, ref.D, ref.U,
            spec['amplitude'], spec['setup']]
    cmd = [str(binary), *map(str, args)]
    write_json(directory / 'request.json', dict(command=cmd, spec=spec, binary_sha256=sha(binary),
                                                 OMP_NUM_THREADS=THREADS))
    print('Native start', spec['name'], flush=True)
    start = perf_counter()
    with (directory / 'stdout.log').open('w') as stream:
        process = subprocess.run(cmd, cwd=directory, stdout=stream, stderr=subprocess.STDOUT,
                                 env=dict(os.environ, OMP_NUM_THREADS=str(THREADS), OMP_DYNAMIC='FALSE'),
                                 timeout=3600)
    write_json(directory / 'execution.json', dict(returncode=process.returncode, elapsed_s=perf_counter()-start))
    if process.returncode:
        raise RuntimeError(f"OpenLB failed for {spec['name']}. See {directory / 'stdout.log'}")
    result = analyze_case(directory, spec, radial)
    write_json(directory / 'analysis.json', result)
    print('Native complete', spec['name'], result['range_kg_m3'], result['checks'], flush=True)
    return spec['name'], result


def qualification(cases, n, setups):
    return {setup: passes(cases[f'{setup}_constant{n}'], cases[f'{setup}_shear{n}']) for setup in setups}


def sources(base):
    files = [HERE / name for name in ('run.py', 'analysis.py', 'protocol.md')]
    files.extend((HERE / 'openlb').iterdir())
    files.extend(HERE.parent / name for name in (
        'stationary_curved_wall/analysis.py', 'stationary_curved_wall/reference.py',
        'transport_verification/analysis.py', 'transport_verification/reference.py',
        'transport_root_cause/flux_audit/half_way/halfway_setup.h',
        'transport_root_cause/hydraulics/stock_omp.py', 'transport_accuracy/run.py'))
    hashes = source_manifest(files)
    for path in hashes:
        target = base / 'source_snapshot' / path
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT / path, target)
    return hashes


def run():
    old, mode = parent_record()
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ') + '-' + uuid4().hex[:8]
    base = OUTPUT / 'runs' / stamp
    base.mkdir(parents=True)
    specs = specifications(128)
    record = dict(schema='pore2process.stock-transport-candidate.v1', run_id=stamp,
        parent_summary=str(PARENT.relative_to(ROOT)), parent_summary_sha256=PARENT_SHA,
        reference=old['reference'], reference_path=str(mode.relative_to(ROOT)),
        source_sha256=sources(base), specifications=specs, threads_per_process=THREADS,
        maximum_concurrent_processes=WORKERS, new_bed_runs=0, new_DWSIM_runs=0)
    write_json(base / 'plan.json', record)
    print('Stock stationary-wall comparison', base, flush=True)
    radial = ref.read_mode(mode)
    metadata = build_case_omp(HERE / 'openlb', ROOT / '.local/adsorption/stock_transport_candidate', base / 'build')
    binary = base / 'binary/capillary3d'
    binary.parent.mkdir()
    shutil.copy2(metadata['binary'], binary)
    with ThreadPoolExecutor(max_workers=WORKERS) as pool:
        cases = dict(pool.map(lambda spec: execute(base, binary, spec, radial), specs))
        fine = qualification(cases, 128, SETUPS)
        followup = specifications(64, [setup for setup, passed in fine.items() if passed])
        cases.update(pool.map(lambda spec: execute(base, binary, spec, radial), followup))
    record.update(specifications=specs+followup, cases=cases, build=metadata, binary_sha256=sha(binary),
        passes_at_128=fine, passes_at_64=qualification(cases, 64, [setup for setup in fine if fine[setup]]),
        vendor=verify_vendor(), new_native_runs=len(cases))
    record['artifacts_sha256'] = {str(path.relative_to(base)): sha(path) for path in base.rglob('*') if path.is_file()}
    path = base / 'summary.json'
    write_json(path, record)
    write_json(OUTPUT / 'latest.json', dict(summary=str(path.relative_to(ROOT))))
    print('Completed stock stationary-wall comparison', path, flush=True)
    return path


def check(path):
    path = Path(path).resolve()
    saved = validate_record(path)
    old, mode = parent_record()
    if (saved['parent_summary_sha256'] != PARENT_SHA or saved['reference'] != old['reference']
            or ROOT / saved['reference_path'] != mode):
        raise ValueError('Comparison used a different parent or reference')
    radial = ref.read_mode(mode)
    binary = path.parent / 'binary/capillary3d'
    if sha(binary) != saved['binary_sha256']:
        raise ValueError('Archived binary hash differs')
    replayed = {}
    for spec in saved['specifications']:
        directory = path.parent / 'native' / spec['name']
        execution = json.loads((directory / 'execution.json').read_text())
        request = json.loads((directory / 'request.json').read_text())
        if execution['returncode'] or request['spec'] != spec or request['binary_sha256'] != sha(binary):
            raise ValueError('Native execution evidence differs')
        replayed[spec['name']] = analyze_case(directory, spec, radial)
    fine = qualification(replayed, 128, SETUPS)
    selected = [setup for setup in fine if fine[setup]]
    if saved['specifications'] != specifications(128)+specifications(64, selected):
        raise ValueError('Executed cases differ from the protocol')
    if (replayed != saved['cases'] or fine != saved['passes_at_128']
            or qualification(replayed, 64, selected) != saved['passes_at_64']):
        raise ValueError('Saved analysis differs from replay')
    verification = dict(summary=str(path.relative_to(ROOT)), summary_sha256=sha(path),
        sources=len(saved['source_sha256']), artifacts=len(saved['artifacts_sha256']),
        native_runs=len(replayed), new_native_runs_in_check=0, vendor=verify_vendor())
    write_json(OUTPUT / 'verification.json', verification)
    print(json.dumps(verification, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', type=Path, help='Recheck archived runs without running OpenLB')
    args = parser.parse_args()
    check(args.check) if args.check else run()


if __name__ == '__main__':
    main()
