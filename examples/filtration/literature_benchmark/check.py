"""Replay Julia, inspect saved DWSIM states and verify the benchmark decisions."""
import argparse
import csv
import json
import math
from pathlib import Path

from examples.filtration.calibration.check import close_tree
from pore2process_demo.io import sha
from examples.filtration.crossflow_recycle.check import verify as verify_parent
from examples.filtration.crossflow_recycle.model import adapter_inputs, topology, validate
from pore2process_demo.bridge import runtime
from pore2process_demo.paths import ROOT; from pore2process_demo.io import read_json as read
from .model import prepare, requests, selections, translate
from .run import comparison_rows


def verify(path, verify_upstream=True):
    path = Path(path).resolve(); base = path.parent; s = read(path)
    if s['status'] != 'passed':
        raise ValueError('Expected a completed benchmark.')
    files = {str(f.relative_to(base)) for f in base.rglob('*') if f.is_file() and
             str(f.relative_to(base)) not in ('summary.json', 'status.json', 'verification.json')}
    if files != set(s['artifact_sha256']):
        raise ValueError('Artifact coverage differs.')
    for name, digest in s['artifact_sha256'].items():
        f = (base / name).resolve()
        if base not in f.parents or sha(f) != digest:
            raise ValueError('Artifact hash mismatch: ' + name)
    snapshot = read(base / 'inputs.json')
    for key, value in snapshot.items():
        close_tree(s[key], value, key)
    if sha(base / 'parent_summary.json') != s['parent_summary_sha256'] or sha(base / 'source_summary.json') != s['source_summary_sha256']:
        raise ValueError('Ancestor snapshot changed.')
    if verify_upstream:
        parent = ROOT / s['parent_summary_path']
        if sha(parent) != s['parent_summary_sha256']:
            raise ValueError('Original ancestor changed.')
        close_tree(verify_parent(parent), s['parent_verification'], 'upstream verification')
    close_tree(read(base / 'parent_verification.json'), s['parent_verification'])
    p = s['protocol']; h = s['hydraulics']; extraction = s['extraction']; source = read(base / 'source_summary.json')
    jl = runtime(); decisions = prepare(jl, extraction, p)
    close_tree(decisions, s['decisions'], 'admissibility decisions')
    expected = list(requests(decisions, p, source))
    if len(expected) != len(s['cases']):
        raise ValueError('Missing/extra process cases.')
    for case, request in zip(s['cases'], expected, strict=True):
        close_tree(case['request'], request, 'case request')
        t = translate(jl, extraction, p, h, source, request)
        close_tree(case['translation'], t, 'Julia translation')
        inputs = adapter_inputs(h, source, t['process'], **request['seed_options'])
        inputs['Scope'] = 'Conditional water-carrier process scenario; literature critical-flux constraint and passive-solute balances are external Julia models. Flow convention is unresolved; no suspension properties or particle species in DWSIM.'
        directory = base / 'process' / request['case_id']
        actual = read(directory / 'summary.json')
        close_tree(case['dwsim'], actual, 'DWSIM evidence')
        close_tree(read(directory / 'inputs.json'), inputs, 'adapter input')
        close_tree(case['validation'], validate(actual, inputs, t['process'], h, source), 'balances')
        close_tree(case['saved_connections'], topology(directory / 'crossflow.dwxmz'), 'topology')
        solute = t['passive_solute_bookkeeping']
        scale = p['module_inlet_concentration_kg_m3'] * t['process']['module_reference_flow_m3_s']
        if max(abs(solute[k]) for k in solute if k.endswith('_residual_kg_s')) > 1e-12 * scale:
            raise ValueError('Passive-solute conservation failed.')
        baseline = next(c for c in s['cases'] if not c['request']['validation_only'] and
            c['request']['flow_convention'] == request['flow_convention'] and
            c['request']['row_id'] == request['row_id'] and c['request']['support'] == 'straight')
        saving = 1 - case['validation']['total_shaft_power_W'] / baseline['validation']['total_shaft_power_W']
        close_tree(case['support_energy_saving_vs_straight_fraction'], saving, 'support saving')
        if request['validation_only']:
            for key in ('feed_W', 'circulation_W', 'cooling_removed_W'):
                if not math.isclose(actual['powers'][key], baseline['dwsim']['powers'][key], rel_tol=h['acceptance']['seed_invariance_relative_power']):
                    raise ValueError('High-circulation seed/tolerance invariance failed.')
        if request['flow_factor'] == 7 and not request['validation_only']:
            literal = next(c for c in s['cases'] if c['request']['flow_factor'] == 1 and
                c['request']['row_id'] == request['row_id'] and c['request']['support'] == request['support'])
            if not math.isclose(actual['powers']['circulation_W'] / literal['dwsim']['powers']['circulation_W'], 7., rel_tol=1e-8):
                raise ValueError('Flow-convention circulation-power oracle failed.')
    close_tree(selections(s['cases'], decisions, p), s['selected_discrete_points'], 'discrete selections')
    with (base / 'comparison.csv').open(newline='') as stream:
        rows = list(csv.DictReader(stream))
    expected_rows = comparison_rows(s['cases'])
    if len(rows) != len(expected_rows):
        raise ValueError('CSV row count differs.')
    for row, expected_row in zip(rows, expected_rows, strict=True):
        for key, value in expected_row.items():
            close_tree(row[key] if isinstance(value, str) else float(row[key]), value, 'CSV.' + key)
    rejected = read(base / 'validation/nonconvergence/rejection.json')
    close_tree(rejected, s['nonconvergence_rejection'])
    if rejected['nonconvergence_rejected'] is not True or rejected['return_code'] == 0:
        raise ValueError('Missing nonconvergence rejection.')
    if 'RECYCLE_NOT_CONVERGED' not in (base / 'validation/nonconvergence/stdout.log').read_text():
        raise ValueError('Missing failed-solve diagnostic.')
    return {'status': 'passed', 'actual_saved_DWSIM_cases_checked': len(s['cases']),
            'design_cases': sum(not c['request']['validation_only'] for c in s['cases']),
            'high_circulation_seed_tolerance_cases': sum(c['request']['validation_only'] for c in s['cases']),
            'artifacts_checked': len(files), 'upstream_reverified': verify_upstream,
            'absolute_flow_basis_validated': False, 'new_external_solves_by_checker': False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('summary', nargs='?', type=Path)
    args = parser.parse_args()
    path = args.summary or Path(read(ROOT / 'results/crossflow_benchmark/latest.json')['summary'])
    print(json.dumps(verify(path), indent=2))


if __name__ == '__main__':
    main()
