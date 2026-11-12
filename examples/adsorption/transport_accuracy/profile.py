"""Native periodic controls with an independent uniform-flow exact solution."""
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
from pore2process_demo.provenance import source_manifest
from pore2process_demo.runtime.openlb import build_case, verify_vendor
from examples.adsorption.column_transport.run import execute, analyze as strict_analysis
from .run import HERE, OUTPUT, source_point, validate_record


def analyze(path):
    path = Path(path)
    s = strict_analysis(path, references=False)
    r = s['raw']
    a = np.genfromtxt(path / 'inventory.csv', delimiter=',', names=True)
    f = np.genfromtxt(path / 'final.csv', delimiter=',', names=True)
    t = r['dt_s'] * r['steps']
    exact_error = None
    if r['uniform']:
        exact = .5 + .25*np.exp(-r['diffusivity_m2_s']*(2*np.pi/r['length_m'])**2*t)*np.sin(
            2*np.pi*(f['x_m']-f['ux_m_s']*t)/r['length_m'])
        exact_error = float(np.sqrt(np.mean((f['concentration_kg_m3']-exact)**2))/.25)
    checks = dict(finite=s['checks']['finite'], kinetic_mass=s['checks']['mass'],
        lower_bound=s['minimum_concentration_kg_m3'] >= .23,
        upper_bound=s['maximum_concentration_kg_m3'] <= .77)
    if exact_error is not None:
        checks['analytical_profile'] = exact_error <= .01
    return dict(checks=checks, passed=all(checks.values()),
        minimum_concentration_kg_m3=s['minimum_concentration_kg_m3'],
        maximum_concentration_kg_m3=s['maximum_concentration_kg_m3'],
        normalized_kinetic_mass_residual=s['normalized_mass_residual'],
        exact_final_profile_RMS_over_initial_amplitude=exact_error,
        actual_duration_s=t, void_transits=t*r['flow_m3_s']/r['fluid_volume_m3'],
        maximum_relative_physical_liquid_inventory_change=float(np.max(abs(a['fluid_tracer_kg']/a['fluid_tracer_kg'][0]-1))),
        raw=r, transport_closure_qualified=False, new_service_prediction=False)


def check():
    p = ROOT / json.loads((OUTPUT / 'profile_latest.json').read_text())['summary']
    s = validate_record(p)
    for key, result in s['cases'].items():
        if analyze(p.parent / 'native' / key) != result:
            raise ValueError('Periodic profile analysis differs.')
    compact = dict(source_summary=str(p.relative_to(ROOT)), source_summary_sha256=sha(p),
        run_id=s['run_id'], cases=s['cases'], actual_native_runs=4,
        native_source_files_checked=len(s['source_sha256']),
        native_artifacts_checked=len(s['artifacts_sha256']),
        new_native_runs_in_check=0, new_DWSIM_runs_in_check=0, service_prediction_qualified=False)
    write_json(HERE / 'reference/profile_reference.json', compact)
    print(json.dumps(compact, indent=2))


def run():
    pointer = ROOT / 'results/adsorption/column_transport/hydraulics/latest.json'
    hydraulic = (ROOT / json.loads(pointer.read_text())['summary']).parent
    parent = validate_record(hydraulic / 'summary.json')
    if not parent['hydraulic_checks_passed']:
        raise ValueError('Qualified hydraulic parent required.')
    rid = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid4().hex[:8]
    base = OUTPUT / 'profile/runs' / rid
    base.mkdir(parents=True)
    sources = source_manifest([Path(__file__), HERE / 'run.py', HERE / 'profile_protocol.md',
        *list((HERE / 'openlb_profile').iterdir()),
        ROOT / 'examples/adsorption/column_transport/run.py', ROOT / 'examples/adsorption/bed.py'])
    for name in sources:
        dest = base / 'source_snapshot' / name
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT / name, dest)
    generated = base / 'generated_source'
    shutil.copytree(HERE / 'openlb_profile', generated)
    config = generated / 'config.mk'
    config.write_text(config.read_text()+'\nCXXFLAGS += -DP2P_SCALAR_MODEL=1\n')
    build = ROOT / '.local/adsorption/transport_accuracy_profile'
    print('Periodic profile run', base, flush=True)
    metadata = build_case(generated, build, base / 'build')
    jobs = []
    for diameter in ('d800', 'd450'):
        point = source_point(hydraulic, diameter, 32)
        h = point['side_m']/32
        dt = .02*h/point['max_speed_m_s']
        for free in (True, False):
            jobs.append(dict(id=diameter+('_free' if free else '_porous'), model=1,
                n=32, cells=1, D_m2_s=1e-9, tau_diffusive=.5+3e-9*dt/h**2,
                end_s=2*point['side_m']/point['speed_m_s'], initial_concentration_kg_m3=.5,
                periodic=True, field='uniform' if free else point['field'],
                speed_m_s=point['speed_m_s'], velocity_scale=point['velocity_scale'],
                outlet='zero', side_m=point['side_m'], source_point=point))
    with ThreadPoolExecutor(max_workers=2) as pool:
        list(pool.map(lambda spec: execute(base, spec, build / 'capillary3d'), jobs))
    record = dict(schema='pore2process.transport-profile-isolation.v1', run_id=rid,
        source_sha256=sources, specifications=jobs, build=metadata, vendor=verify_vendor(),
        hydraulic_summary=str((hydraulic / 'summary.json').relative_to(ROOT)),
        hydraulic_summary_sha256=sha(hydraulic / 'summary.json'),
        cases={p['id']: analyze(base / 'native' / p['id']) for p in jobs},
        new_native_runs=4, new_DWSIM_runs=0, service_prediction_qualified=False)
    record['artifacts_sha256'] = {str(p.relative_to(base)): sha(p) for p in base.rglob('*') if p.is_file()}
    write_json(base / 'summary.json', record)
    write_json(OUTPUT / 'profile_latest.json', dict(summary=str((base / 'summary.json').relative_to(ROOT))))
    print('Periodic profile controls completed', base, flush=True)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--check', action='store_true')
    args = p.parse_args()
    check() if args.check else run()


if __name__ == '__main__':
    main()
