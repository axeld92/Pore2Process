"""Confirm spatial resolution of the qualified tracer accuracy at 160×1024."""
import argparse
from datetime import datetime, timezone
from pathlib import Path
import shutil
from uuid import uuid4

from examples.adsorption.stationary_curved_wall import reference as ref
from examples.adsorption.transport_accuracy.run import validate_record
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT
from pore2process_demo.provenance import source_manifest

from .run import HERE, OUTPUT, SPACE_LIMIT, archived, evaluate, execute, julia_command, passes, refinement, spec

ORIGINAL = OUTPUT / 'runs/20261111T153035Z-2fda74bc/summary.json'
ORIGINAL_SHA = '400d11c320252b0cb5c54c08fa4dd843b05568e0118f2be6256d2b5fbdb206f2'
SETTINGS = spec(320, 2048, 1e-4)


def control():
    if sha(ORIGINAL) != ORIGINAL_SHA:
        raise ValueError('Initial disk record differs')
    record = validate_record(ORIGINAL)
    qualification, mode_path, samples = archived()
    if record['reference'] != qualification:
        raise ValueError('Initial disk record used a different radial reference')
    return record, mode_path, samples


def compare(directory, original, mode_path, samples):
    old = original['analysis']
    name = old['selected']
    if name != 'fine_half_dt' or old['space']['medium_to_fine']['passed']:
        raise ValueError('The archived study does not match the confirmation protocol')
    case_dir = directory / 'varying'
    runtime = original['julia_command']
    if (read_json(case_dir / 'request.json') != dict(command=julia_command(case_dir, SETTINGS, runtime), case=SETTINGS)
            or read_json(case_dir / 'execution.json')['returncode']):
        raise ValueError('Saved execution differs from the confirmation protocol')
    mode = ref.read_mode(mode_path)
    working_case = old['cases'][name]
    working, fixed = evaluate(ORIGINAL.parent / 'cases' / name, working_case['settings'], mode, samples)
    if working != working_case['metrics']:
        raise ValueError('Working-mesh results differ from the archived analysis')
    metrics, fine = evaluate(case_dir, SETTINGS, mode, samples)
    space = refinement(fixed, fine, samples, SPACE_LIMIT)
    qualified = (passes(working) and passes(metrics) and space['passed'] and old['time']['passed']
                 and passes(old['cases']['uniform']['metrics']))
    return dict(working_case=name, working_mesh=working_case['settings'], confirmation_mesh=SETTINGS,
                confirmation_metrics=metrics, spatial_change=space, disk_qualified=qualified)


def run():
    original, mode, samples = control()
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid4().hex[:8]
    directory = OUTPUT / 'confirmations' / stamp
    directory.mkdir(parents=True)
    hashes = source_manifest([HERE / 'confirm.py', HERE / 'refinement_note.md'])
    for path in hashes:
        target = directory / 'source_snapshot' / path
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT / path, target)
    record = dict(schema='pore2process.disk-face-flow-confirmation.v1', run_id=stamp,
                  initial_summary=str(ORIGINAL.relative_to(ROOT)), initial_summary_sha256=ORIGINAL_SHA,
                  source_sha256=hashes, reference_sha256=sha(mode), settings=SETTINGS,
                  new_bed_runs=0, new_openlb_runs=0, new_DWSIM_runs=0)
    write_json(directory / 'plan.json', record)
    execute(directory / 'varying', SETTINGS, original['julia_command'])
    record['analysis'] = compare(directory, original, mode, samples)
    record['artifacts_sha256'] = {str(p.relative_to(directory)): sha(p)
                                 for p in directory.rglob('*') if p.is_file()}
    path = directory / 'summary.json'
    write_json(path, record)
    write_json(OUTPUT / 'latest_confirmation.json', dict(summary=str(path.relative_to(ROOT))))
    print('Completed spatial confirmation', path, record['analysis']['spatial_change'],
          'qualified', record['analysis']['disk_qualified'], flush=True)


def check(path):
    path = path.resolve()
    saved = validate_record(path)
    original, mode, samples = control()
    if saved['initial_summary_sha256'] != ORIGINAL_SHA or sha(mode) != saved['reference_sha256']:
        raise ValueError('Confirmation reference differs')
    result = compare(path.parent, original, mode, samples)
    if result != saved['analysis']:
        raise ValueError('Replayed spatial confirmation differs')
    report = dict(summary=str(path.relative_to(ROOT)), summary_sha256=sha(path), sources=len(saved['source_sha256']),
                  artifacts=len(saved['artifacts_sha256']), disk_qualified=result['disk_qualified'])
    write_json(OUTPUT / 'confirmation_verification.json', report)
    print(report)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', type=Path)
    args = parser.parse_args()
    check(args.check) if args.check else run()
