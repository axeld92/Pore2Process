"""Verify saved native evidence and replay Julia analysis, without simulators."""
import argparse
import json
from pathlib import Path

from pore2process_demo.paths import ROOT
from pore2process_demo.io import sha, write_json
from pore2process_demo.runtime.openlb import verify_vendor
from .run import HERE, OUTPUT, validate_record
from .assess import evaluate


def check(path):
    path = Path(path).resolve()
    s = validate_record(path)
    native_path = ROOT / s['native_summary']
    if sha(native_path) != s['native_summary_sha256']:
        raise ValueError('Changed native parent summary.')
    native = validate_record(native_path)
    hydraulic_path = ROOT / native['hydraulic_summary']
    if sha(hydraulic_path) != native['hydraulic_summary_sha256']:
        raise ValueError('Changed hydraulic parent summary.')
    validate_record(hydraulic_path)
    actual, decision = evaluate(native_path.parent)
    if actual != s['cases'] or decision != s['decision']:
        raise ValueError('Saved analysis does not replay.')
    result = dict(schema=s['schema'], assessment_run_id=s['run_id'],
        source_summary=str(path.relative_to(ROOT)), source_summary_sha256=sha(path),
        native_run_id=native['run_id'], actual_native_runs=native['new_native_runs'],
        actual_new_DWSIM_runs=0, decision=decision,
        cases={key: {k: c[k] for k in ('engineering_checks', 'engineering_screen_passed',
            'actual_duration_s', 'physical_void_time_s', 'observed_transits',
            'maximum_average_liquid_concentration_relative_change', 'final_quarter_extrema_drift',
            'final_spatial_diagnostic', 'final_wall_population_change_over_physical_liquid_plateau_inventory',
            'RTD', 'observed_outlet_moments', 'imposed_inlet_moments',
            'RTD_unavailable_reason', 'absolute_outlet_plateau_relative_error', 'extrema_locations')}
            | dict(uniform_maximum_departure=c['strict_diagnostic']['uniform_maximum_departure'],
                   minimum_concentration=c['strict_diagnostic']['minimum_concentration_kg_m3'],
                   maximum_concentration=c['strict_diagnostic']['maximum_concentration_kg_m3'],
                   normalized_kinetic_mass_residual=c['strict_diagnostic']['normalized_mass_residual'])
            for key, c in actual.items()},
        native_source_files_checked=len(native['source_sha256']),
        native_artifacts_checked=len(native['artifacts_sha256']),
        assessment_source_files_checked=len(s['source_sha256']),
        assessment_artifacts_checked=len(s['artifacts_sha256']),
        new_native_runs_in_check=0, new_DWSIM_runs_in_check=0,
        Julia_RTD_analyses=sum(c['RTD'] is not None for c in actual.values()), vendor=verify_vendor())
    write_json(HERE / ('reference/' + native['stage'] + '_reference.json'), result)
    write_json(path.parent / 'verification.json', {k: v for k, v in result.items() if k != 'cases'})
    print(json.dumps({k: v for k, v in result.items() if k not in ('cases', 'vendor')}, indent=2))
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--assessment', type=Path)
    parser.add_argument('--all', action='store_true')
    args = parser.parse_args()
    if args.all:
        reports = []
        for name in ('screen', 'smooth_input', 'regularized_inlet'):
            pointer = OUTPUT / ('assessment_' + name + '_latest.json')
            if not pointer.is_file():
                raise ValueError('Missing assessment stage: ' + name)
            reports.append(check(ROOT / json.loads(pointer.read_text())['summary']))
        from .profile import check as profile_check
        profile_pointer = json.loads((OUTPUT / 'profile_latest.json').read_text())
        profile_path = ROOT / profile_pointer['summary']
        profile_native = validate_record(profile_path)
        hydraulic = ROOT / profile_native['hydraulic_summary']
        if sha(hydraulic) != profile_native['hydraulic_summary_sha256']:
            raise ValueError('Changed periodic-profile hydraulic parent.')
        validate_record(hydraulic)
        for spec in profile_native['specifications']:
            request = json.loads((profile_path.parent / 'native' / spec['id'] / 'request.json').read_text())
            if spec['field'] != 'uniform' and sha(Path(spec['field'])) != request['native_field_sha256']:
                raise ValueError('Changed periodic-profile velocity field.')
        profile_check()
        write_json(OUTPUT / 'verification.json', dict(stages=[r['native_run_id'] for r in reports],
            periodic_profile_run_id=profile_native['run_id'],
            actual_native_runs=sum(r['actual_native_runs'] for r in reports)+profile_native['new_native_runs'],
            new_native_runs_in_check=0, new_DWSIM_runs_in_check=0,
            service_prediction_qualified=False, checker_source_sha256=sha(Path(__file__))))
    else:
        path = args.assessment or ROOT / json.loads((OUTPUT / 'assessment_latest.json').read_text())['summary']
        check(path)


if __name__ == '__main__':
    main()
