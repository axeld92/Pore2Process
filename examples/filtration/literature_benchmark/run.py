"""Run the real Julia benchmark and conditional stock-DWSIM process cases."""
from pore2process_demo.provenance import implementation_files, source_manifest
import argparse
import copy
import csv
from datetime import datetime, timezone
from pathlib import Path
import shutil
from time import perf_counter
from uuid import uuid4

from examples.filtration.calibration.audit_crossflow_sources import audit
from pore2process_demo.io import sha
from examples.filtration.crossflow_recycle import dwsim
from examples.filtration.crossflow_recycle.check import verify as verify_parent
from examples.filtration.crossflow_recycle.model import adapter_inputs, topology, validate
from examples.filtration.crossflow_recycle.run import check_protocol
from pore2process_demo.io import write_json
from pore2process_demo.bridge import runtime, runtime_versions
from pore2process_demo.paths import ROOT; from pore2process_demo.io import read_json as read
from .model import prepare, requests, selections, translate


def comparison_rows(cases):
    rows = []
    for c in cases:
        q = c['request']; t = c['translation']; v = c['validation']
        rows.append({'case_id': q['case_id'], 'support': q['support'],
            'flow_convention': q['flow_convention'], 'row_id': q['row_id'],
            'validation_only': int(q['validation_only']),
            'velocity_m_s': t['source_point']['velocity_m_s'],
            'critical_flux_m_s': t['source_point']['critical_flux_m_s'],
            'channel_drop_Pa': t['process']['channel_pressure_drop_Pa'],
            'module_reference_flow_m3_s': t['process']['module_reference_flow_m3_s'],
            'source_axial_energy_J_m3': t['source_point']['source_specific_axial_energy_J_m3'],
            'scenario_axial_energy_J_m3': t['source_point']['scenario_specific_axial_energy_J_m3'],
            'TMP_Pa': t['membrane']['pressure_drop_Pa'],
            'fresh_concentration_kg_m3': t['passive_solute_bookkeeping']['fresh_concentration_kg_m3'],
            'actual_feed_shaft_W': c['dwsim']['powers']['feed_W'],
            'actual_circulation_shaft_W': c['dwsim']['powers']['circulation_W'],
            'actual_cooling_removed_W': c['dwsim']['powers']['cooling_removed_W'],
            'actual_total_shaft_W': v['total_shaft_power_W'],
            'actual_specific_shaft_energy_kWh_m3': v['specific_shaft_energy_kWh_m3_reference_volume'],
            'support_energy_saving_vs_straight_fraction': c['support_energy_saving_vs_straight_fraction'],
            'recycle_iterations': c['dwsim']['recycle']['iterations_taken']})
    return rows


def execute(output, parent_path, require_sources=False):
    p = read(ROOT / 'examples/filtration/literature_benchmark/protocol.json')
    extraction = read(ROOT / p['extraction'])
    source_audit = audit(extraction, require_sources)
    parent_path = Path(parent_path).resolve()
    parent = read(parent_path)
    source_path = parent_path.parent / 'source_summary.json'
    source = read(source_path)
    hydraulics = copy.deepcopy(parent['protocol'])
    check_protocol(hydraulics, source)
    if p['target_flux_LMH'] != hydraulics['permeate_flux_LMH'] or p['process_area_m2'] != hydraulics['membrane_area_m2']:
        raise ValueError('Revalidate membrane translation before changing the inherited flux or area.')
    hydraulics['recycle']['maximum_iterations'] = p['maximum_recycle_iterations']
    start = perf_counter(); jl = runtime(); cold = perf_counter() - start
    print('Rechecking saved recycle, filtration and native OpenLB evidence...', flush=True)
    parent_verification = verify_parent(parent_path)
    decisions = prepare(jl, extraction, p)
    run_id = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ') + '-' + uuid4().hex[:8]
    base = Path(output) / 'runs' / run_id
    base.mkdir(parents=True)
    write_json(base / 'status.json', {'status': 'running', 'run_id': run_id})
    try:
        shutil.copy2(parent_path, base / 'parent_summary.json')
        shutil.copy2(source_path, base / 'source_summary.json')
        write_json(base / 'source_audit.json', source_audit)
        write_json(base / 'parent_verification.json', parent_verification)
        adapter = dwsim.build(base / 'adapter_build')
        paths = list((ROOT / 'examples/filtration/literature_benchmark').glob('*.py')) + [ROOT / 'examples/filtration/literature_benchmark/protocol.json']
        paths += [ROOT / name for name in ('julia/src/CrossflowBenchmark.jl', 'julia/src/Pore2ProcessDemo.jl',
            'julia/src/CrossflowRecycle.jl', 'julia/src/GradedSupport.jl',
            'examples/filtration/calibration/audit_crossflow_sources.py', 'examples/filtration/crossflow_recycle/model.py',
            'examples/filtration/crossflow_recycle/dwsim.py', 'examples/filtration/crossflow_recycle/dwsim_runner/Program.cs',
            'examples/filtration/crossflow_recycle/dwsim_runner/CrossflowRunner.csproj', 'src/pore2process_demo/bridge.py',
            'uv.lock', 'julia/Project.toml', 'julia/Manifest.toml')]
        paths = implementation_files(paths)
        snapshot = {'protocol': p, 'hydraulics': hydraulics, 'extraction': extraction,
            'parent_summary_path': str(parent_path.relative_to(ROOT)),
            'parent_summary_sha256': sha(parent_path), 'source_summary_sha256': sha(source_path),
            'extraction_sha256': sha(ROOT / p['extraction']),
            'implementation_sha256': source_manifest(sorted(paths)),
            'runtime_sha256': {name: sha(dwsim.DISTRIBUTION / name) for name in
                ('DWSIM.UnitOperations.dll', 'DWSIM.Thermodynamics.dll', 'DWSIM.FlowsheetSolver.dll',
                 'DWSIM.Automation.FluentAPI.dll', 'DWSIM.Interfaces.dll')}}
        write_json(base / 'inputs.json', snapshot)
        cases = []
        for request in requests(decisions, p, source):
            t = perf_counter()
            translated = translate(jl, extraction, p, hydraulics, source, request)
            elapsed = perf_counter() - t
            inputs = adapter_inputs(hydraulics, source, translated['process'], **request['seed_options'])
            inputs['Scope'] = 'Conditional water-carrier process scenario; literature critical-flux constraint and passive-solute balances are external Julia models. Flow convention is unresolved; no suspension properties or particle species in DWSIM.'
            print('Solving/reloading actual DWSIM: ' + request['case_id'], flush=True)
            directory = base / 'process' / request['case_id']
            actual = dwsim.solve(directory, inputs)
            checks = validate(actual, inputs, translated['process'], hydraulics, source)
            baseline = next((c for c in cases if not c['request']['validation_only'] and
                c['request']['flow_convention'] == request['flow_convention'] and
                c['request']['row_id'] == request['row_id'] and c['request']['support'] == 'straight'), None)
            saving = 0. if baseline is None else 1 - checks['total_shaft_power_W'] / baseline['validation']['total_shaft_power_W']
            cases.append({'request': request, 'translation': translated, 'dwsim': actual,
                'validation': checks, 'saved_connections': topology(directory / 'crossflow.dwxmz'),
                'support_energy_saving_vs_straight_fraction': saving, 'Julia_translation_s': elapsed})
        first = cases[0]
        bad = dict(first['dwsim']['inputs']); bad['Maximum_iterations'] = 1
        bad['Seed_module_mass_flow_kg_s'] *= .2
        rejection = dwsim.reject_nonconvergence(base / 'validation/nonconvergence', bad)
        rows = comparison_rows(cases)
        with (base / 'comparison.csv').open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
        selected = selections(cases, decisions, p)
        from .plot import plot
        plot(extraction, source_audit, cases, p, base / 'benchmark.png')
        shutil.copy2(base / 'process' / first['request']['case_id'] / 'flowsheet.png', base / 'flowsheet.png')
        result = {'schema_version': 1, 'status': 'passed', 'run_id': run_id,
            'timestamp_utc': datetime.now(timezone.utc).isoformat(), **snapshot,
            'parent_verification': parent_verification, 'adapter': adapter, 'versions': runtime_versions(),
            'Julia_initialization_s': cold, 'decisions': decisions, 'cases': cases,
            'selected_discrete_points': selected, 'nonconvergence_rejection': rejection,
            'scope': 'Conditional critical-flux-constrained water-carrier scenarios with externally balanced passive solute. No fitted fouling kinetics or experimentally validated absolute flow convention.',
            'artifact_sha256': {str(f.relative_to(base)): sha(f) for f in base.rglob('*') if f.is_file() and f.name != 'status.json'}}
        write_json(base / 'summary.json', result)
        from .check import verify
        verification = verify(base / 'summary.json', verify_upstream=False)
        write_json(base / 'verification.json', verification)
        write_json(base / 'status.json', {'status': 'passed', 'run_id': run_id})
        write_json(Path(output) / 'latest.json', {'summary': str((base / 'summary.json').resolve())})
        print('Passed: ' + str(base / 'summary.json'), flush=True)
        return result
    except BaseException as error:
        write_json(base / 'status.json', {'status': 'failed', 'run_id': run_id, 'error': str(error)})
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT / 'results/crossflow_benchmark')
    parser.add_argument('--parent', type=Path)
    parser.add_argument('--require-sources', action='store_true')
    args = parser.parse_args()
    parent = args.parent or Path(read(ROOT / 'results/crossflow_recycle/latest.json')['summary'])
    execute(args.output, parent, args.require_sources)


if __name__ == '__main__':
    main()
