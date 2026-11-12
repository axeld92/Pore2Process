"""Archive a fresh assessment of immutable native outputs; no new simulators."""
from datetime import datetime, timezone
import argparse
import json
from pathlib import Path
import shutil
from time import perf_counter
from uuid import uuid4

from pore2process_demo.paths import ROOT
from pore2process_demo.io import sha, write_json
from pore2process_demo.provenance import source_manifest
from pore2process_demo.bridge import runtime_versions
from .run import HERE, OUTPUT, validate_record
from .analysis import analyze, service_gate, refinement_metrics


def evaluate(native):
    s = validate_record(native / 'summary.json')
    cases = {p['id']: analyze(native / 'native' / p['id']) for p in s['specifications']}
    gate = service_gate(cases, [p['id'] for p in s['specifications']])
    if s['stage'] == 'smooth_input':
        gate['ramp_duration_comparisons'] = {}
        for diameter in ('d800', 'd450'):
            selected = [p['id'] for p in s['specifications'] if p['particle'] == diameter]
            if len(selected) == 2:
                gate['ramp_duration_comparisons'][diameter] = refinement_metrics(*(cases[k] for k in selected))
        gate['eligible_for_closure_assessment'] &= bool(gate['ramp_duration_comparisons']) and all(
            v['passed'] for v in gate['ramp_duration_comparisons'].values())
    return cases, gate


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--run', type=Path, required=True)
    args = parser.parse_args()
    native = args.run.resolve()
    native_record = validate_record(native / 'summary.json')
    rid = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ') + '-' + uuid4().hex[:8]
    out = OUTPUT / 'assessments' / rid
    out.mkdir(parents=True)
    files = [Path(__file__), HERE / 'analysis.py', HERE / 'run.py', HERE / 'protocol.md',
             ROOT / 'examples/adsorption/column_transport/run.py',
             ROOT / 'examples/adsorption/packed_bed_tracer/bridge.py',
             ROOT / 'julia/adsorption/TracerContact.jl']
    sources = source_manifest(files)
    for name in sources:
        dest = out / 'source_snapshot' / name
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT / name, dest)
    start = perf_counter()
    cases, gate = evaluate(native)
    s = dict(schema='pore2process.transport-engineering-accuracy.v1', run_id=rid,
             native_stage=native_record['stage'],
             source_sha256=sources, native_summary=str((native / 'summary.json').relative_to(ROOT)),
             native_summary_sha256=sha(native / 'summary.json'), cases=cases, decision=gate,
             new_native_runs=0, new_DWSIM_runs=0, new_adsorption_predictions=0,
             Julia_RTD_analyses=sum(c['RTD'] is not None for c in cases.values()),
             versions=runtime_versions(), elapsed_s=perf_counter() - start)
    s['artifacts_sha256'] = {str(p.relative_to(out)): sha(p) for p in out.rglob('*') if p.is_file()}
    write_json(out / 'summary.json', s)
    write_json(OUTPUT / 'assessment_latest.json', dict(summary=str((out / 'summary.json').relative_to(ROOT))))
    write_json(OUTPUT / ('assessment_' + native_record['stage'] + '_latest.json'),
               dict(summary=str((out / 'summary.json').relative_to(ROOT))))
    print('Assessment saved', out, flush=True)
    print(json.dumps({k: dict(checks=c['engineering_checks'], RTD=c['RTD'],
        local_departure=c['strict_diagnostic']['uniform_maximum_departure'],
        average_liquid_change=c['maximum_average_liquid_concentration_relative_change'],
        late_drift=c['final_quarter_extrema_drift']) for k, c in cases.items()}, indent=2))


if __name__ == '__main__':
    main()
