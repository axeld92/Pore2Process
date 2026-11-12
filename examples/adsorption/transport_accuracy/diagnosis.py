"""Bounded native diagnoses using the preserved profile binary, no solver edits."""
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import argparse
import csv
import json
from pathlib import Path
import shutil
from uuid import uuid4

import numpy as np

from pore2process_demo.paths import ROOT
from pore2process_demo.io import sha, write_json
from pore2process_demo.provenance import source_manifest
from pore2process_demo.runtime.openlb import verify_vendor
from examples.adsorption.column_transport.run import execute
from .profile import analyze
from .run import HERE, OUTPUT, validate_record, source_point


def zero_velocity_copy(source, target):
    with Path(source).open() as inp, Path(target).open('w', newline='') as out:
        reader = csv.DictReader(inp)
        writer = csv.DictWriter(out, fieldnames=reader.fieldnames)
        writer.writeheader()
        for row in reader:
            for name in ('ux_m_s', 'uy_m_s', 'uz_m_s'):
                row[name] = '0'
            writer.writerow(row)


def field_audit(point, final_path):
    native = np.genfromtxt(point['field'], delimiter=',', names=True)
    scalar = np.genfromtxt(final_path, delimiter=',', names=True)
    n = round(len(scalar) ** (1/3))
    selected = scalar['material'] != 2
    np.testing.assert_array_equal(native[['i', 'j', 'k']], scalar[selected][['i', 'j', 'k']])
    errors = {axis: float(np.max(abs(native[axis+'_m']-scalar[selected][axis+'_m']))) for axis in ('x', 'y', 'z')}
    for axis in ('x', 'y', 'z'):
        np.testing.assert_allclose(native['u'+axis+'_m_s']*point['velocity_scale'],
                                   scalar[selected]['u'+axis+'_m_s'], rtol=1e-11, atol=1e-16)
    fluid = selected.reshape((n, n, n))
    u = np.stack([scalar['u'+axis+'_m_s'].reshape(fluid.shape) for axis in ('x', 'y', 'z')], axis=-1)
    h = point['side_m']/n
    div = np.zeros(fluid.shape)
    near = np.zeros(fluid.shape, dtype=bool)
    for axis in range(3):
        plus = np.roll(fluid, -1, axis=axis)
        minus = np.roll(fluid, 1, axis=axis)
        near |= (~plus | ~minus) & fluid
        positive_face = .5*(u[..., axis]+np.roll(u[..., axis], -1, axis=axis))*fluid*plus
        div += (positive_face-np.roll(positive_face, 1, axis=axis))/h
    c = scalar['concentration_kg_m3'].reshape(fluid.shape)
    violations = ((c < .23) | (c > .77)) & fluid
    density = native['density_kg_m3']
    return dict(coordinate_maximum_error_m=errors, fluid_node_IDs_match=True,
        physical_velocities_match=True, fluid_nodes=int(np.sum(fluid)),
        native_maximum_relative_density_departure=float(np.max(abs(density/1000-1))),
        velocity_divergence_proxy_RMS_s_inv=float(np.sqrt(np.mean(div[fluid]**2))),
        velocity_divergence_proxy_maximum_s_inv=float(np.max(abs(div[fluid]))),
        proxy_definition='Centered Cartesian face interpolation; zero faces adjacent to solid. Not OpenLB population flux or a local-continuity acceptance test.',
        final_fraction_fluid_nodes_outside_engineering_bounds=float(np.mean(violations[fluid])),
        final_fraction_violating_nodes_next_to_solid=None if not np.any(violations) else float(np.mean(near[violations])),
        instantaneous_diffusion_length_over_one_reference_transit_m=float(np.sqrt(2e-9*point['side_m']/point['speed_m_s'])),
        grid_spacing_m=h)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--check', action='store_true')
    args = p.parse_args()
    if args.check:
        check()
        return
    parent_path = ROOT / json.loads((OUTPUT / 'profile_latest.json').read_text())['summary']
    parent = validate_record(parent_path)
    hydraulic_path = ROOT / parent['hydraulic_summary']
    assert sha(hydraulic_path) == parent['hydraulic_summary_sha256']
    validate_record(hydraulic_path)
    coarse = next(dict(x) for x in parent['specifications'] if x['id']=='d800_porous')
    binary = ROOT / '.local/adsorption/transport_accuracy_profile/capillary3d'
    original_request = json.loads((parent_path.parent/'native/d800_porous/request.json').read_text())
    if sha(binary) != original_request['binary_sha256']:
        raise ValueError('Preserved profile executable changed; do not silently rebuild it.')
    rid = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid4().hex[:8]
    base = OUTPUT/'diagnosis/runs'/rid
    base.mkdir(parents=True)
    sources = source_manifest([Path(__file__), HERE/'diagnosis_protocol.md', HERE/'profile.py',
        HERE/'run.py', *list((HERE/'openlb_profile').iterdir()),
        ROOT/'examples/adsorption/column_transport/run.py'])
    for name in sources:
        target = base/'source_snapshot'/name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT/name, target)
    zero = base/'zero_velocity.csv'
    zero_velocity_copy(coarse['field'], zero)
    jobs = []
    for name, factor in [('no_advection', 1.), ('half_dt', .5), ('tenth_dt', .1)]:
        spec = dict(coarse)
        spec.update(id=name, tau_diffusive=.5+(coarse['tau_diffusive']-.5)*factor)
        if name=='no_advection':
            spec['field'] = str(zero)
        jobs.append(spec)
    fine = source_point(hydraulic_path.parent, 'd800', 48)
    base_dt = parent['cases']['d800_porous']['raw']['dt_s']
    spec = dict(coarse)
    spec.update(id='grid48_matched_dt', n=48, field=fine['field'], source_point=fine,
        speed_m_s=fine['speed_m_s'], tau_diffusive=.5+3e-9*(base_dt/2)/(fine['side_m']/48)**2)
    jobs.append(spec)
    plan = dict(run_id=rid, source_sha256=sources, specifications=jobs,
        parent_summary=str(parent_path.relative_to(ROOT)), parent_summary_sha256=sha(parent_path),
        hydraulic_summary=str(hydraulic_path.relative_to(ROOT)), hydraulic_summary_sha256=sha(hydraulic_path),
        binary_sha256=sha(binary), binary=str(binary), new_native_runs=4,
        new_DWSIM_runs=0, new_adsorption_predictions=0, numerical_source_changed=False)
    write_json(base/'plan.json', plan)
    print('Bounded diagnosis', base, flush=True)
    with ThreadPoolExecutor(max_workers=2) as pool:
        list(pool.map(lambda job: execute(base, job, binary), jobs))
    plan['cases'] = {s['id']:analyze(base/'native'/s['id']) for s in jobs}
    plan['original_case'] = parent['cases']['d800_porous']
    plan['field_audits'] = {s['id']:field_audit(s['source_point'], base/'native'/s['id']/'final.csv')
        for s in jobs if s['id']!='no_advection'}
    plan['vendor'] = verify_vendor()
    plan['service_prediction_qualified'] = False
    plan['artifacts_sha256'] = {str(p.relative_to(base)):sha(p) for p in base.rglob('*') if p.is_file()}
    write_json(base/'summary.json', plan)
    write_json(OUTPUT/'diagnosis_latest.json', dict(summary=str((base/'summary.json').relative_to(ROOT))))
    print('Bounded diagnosis completed', base, {k:v['passed'] for k,v in plan['cases'].items()}, flush=True)


def check():
    path = ROOT / json.loads((OUTPUT/'diagnosis_latest.json').read_text())['summary']
    s = validate_record(path)
    for role in ('parent', 'hydraulic'):
        p = ROOT/s[role+'_summary']
        if sha(p) != s[role+'_summary_sha256']:
            raise ValueError('Changed parent summary.')
        validate_record(p)
    for spec in s['specifications']:
        dest = path.parent/'native'/spec['id']
        request = json.loads((dest/'request.json').read_text())
        if sha(Path(spec['field'])) != request['native_field_sha256']:
            raise ValueError('Changed native velocity input.')
        if analyze(dest) != s['cases'][spec['id']]:
            raise ValueError('Diagnosis does not replay.')
        if spec['id']!='no_advection' and field_audit(spec['source_point'], dest/'final.csv') != s['field_audits'][spec['id']]:
            raise ValueError('Field audit does not replay.')
    compact = {k:s[k] for k in ('run_id','parent_summary','parent_summary_sha256','original_case',
        'cases','field_audits','new_native_runs','new_DWSIM_runs','new_adsorption_predictions',
        'service_prediction_qualified','numerical_source_changed','vendor')}
    compact.update(source_summary=str(path.relative_to(ROOT)), source_summary_sha256=sha(path),
        source_files_checked=len(s['source_sha256']), artifacts_checked=len(s['artifacts_sha256']),
        new_native_runs_in_check=0, new_DWSIM_runs_in_check=0)
    write_json(HERE/'reference/diagnosis_reference.json', compact)
    print(json.dumps(compact, indent=2))


if __name__=='__main__':
    main()
