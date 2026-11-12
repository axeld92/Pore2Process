"""Execute and replay a fixed-physics, fixed-mesh OpenLB disk at half its time step."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import os
from pathlib import Path
import shutil
import subprocess
from time import perf_counter
from uuid import uuid4

import numpy as np

from examples.adsorption.stationary_curved_wall import reference as ref
from examples.adsorption.stationary_curved_wall.analysis import rms
from examples.adsorption.stock_transport_candidate.run import parent_record
from examples.adsorption.transport_accuracy.run import validate_record
from examples.adsorption.transport_verification.analysis import EXCURSION_LIMIT, RMS_LIMIT, SHELL_M, load
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT
from pore2process_demo.provenance import source_manifest
from pore2process_demo.runtime.openlb import verify_vendor

HERE = Path(__file__).parent
OUTPUT = ROOT / 'results/adsorption/stationary_wall_time_step'
CONTROL = ROOT / 'results/adsorption/stock_transport_candidate/runs/20261111T133159Z-7f69a322/summary.json'
CONTROL_SHA = '99f9aaa64ecc977b614d4b1e25051e93e3539f57e72446741729d82961e2c242'
NAMES = {'constant': 'regularized_halfway_constant128', 'varying': 'regularized_halfway_shear128'}


def control_record():
    if sha(CONTROL) != CONTROL_SHA:
        raise ValueError('Archived regularized TRT comparison differs')
    control = validate_record(CONTROL)
    stationary, reference_path = parent_record()
    if control['reference'] != stationary['reference'] or ROOT / control['reference_path'] != reference_path:
        raise ValueError('Archived Julia reference differs')
    binary = CONTROL.parent / 'binary/capillary3d'
    if sha(binary) != control['binary_sha256']:
        raise ValueError('Archived OpenLB binary differs')
    return control, reference_path, binary


def jobs(control):
    result = {}
    for name, archived_name in NAMES.items():
        raw = control['cases'][archived_name]['raw']
        result[name] = dict(control=archived_name, steps=2*raw['steps'], dt_s=raw['dt_s']/2)
    # The original observation is one fine step before the automatic midpoint.
    steps = 2*(control['cases'][NAMES['varying']]['raw']['steps']//2)
    result['varying_midpoint'] = dict(result['varying'], steps=steps)
    return result


def request(binary, spec, control):
    raw = control['cases'][spec['control']]['raw']
    prior = read_json(CONTROL.parent / 'native' / spec['control'] / 'request.json')
    args = [raw['kind'], raw['n'], raw['nz'], spec['steps'], raw['side_m'], spec['dt_s'],
            raw['diffusivity_m2_s'], raw['speed_m_s'], raw['amplitude_kg_m3'], raw['setup']]
    return dict(command=[str(binary), *map(str, args)], spec=spec,
                binary_sha256=control['binary_sha256'], OMP_NUM_THREADS=prior['OMP_NUM_THREADS'])


def execute(base, name, spec, control):
    directory = base / 'native' / name
    directory.mkdir(parents=True)
    req = request(base / 'binary/capillary3d', spec, control)
    write_json(directory / 'request.json', req)
    print(f"OpenLB start {name}, {spec['steps']} steps", flush=True)
    start = perf_counter()
    with (directory / 'stdout.log').open('w') as stream:
        process = subprocess.run(req['command'], cwd=directory, stdout=stream, stderr=subprocess.STDOUT,
            env=dict(os.environ, OMP_NUM_THREADS=str(req['OMP_NUM_THREADS']), OMP_DYNAMIC='FALSE'), timeout=3600)
    write_json(directory / 'execution.json', dict(returncode=process.returncode, elapsed_s=perf_counter()-start))
    if process.returncode:
        raise RuntimeError(f'OpenLB failed for {name}. See {directory / "stdout.log"}')
    print(f'OpenLB complete {name}, {perf_counter()-start:.1f} s', flush=True)


def verify_native(base, name, spec, control):
    directory = base / 'native' / name
    old = CONTROL.parent / 'native' / spec['control']
    if read_json(directory / 'request.json') != request(base / 'binary/capillary3d', spec, control):
        raise ValueError('Native request differs from the prescribed comparison')
    if read_json(directory / 'execution.json')['returncode']:
        raise ValueError('Native execution failed')
    raw, mask, velocity, fields, history, coords = load(directory)
    archived, old_mask, old_velocity, old_fields, _, _ = load(old)
    changing = {'steps', 'dt_s', 'duration_s', 'tau_even', 'tau_odd', 'peak_step', 'maximum_excursion_kg_m3', 'wall_s'}
    for key in archived.keys()-changing:
        if raw[key] != archived[key]:
            raise ValueError('Unmatched physical input or method: '+key)
    dt, steps, dx = raw['dt_s'], raw['steps'], raw['dx_m']
    if (dt, steps, raw['duration_s']) != (spec['dt_s'], spec['steps'], spec['steps']*spec['dt_s']):
        raise ValueError('Incorrect time step or duration')
    expected_time = ref.TIMES[0 if name == 'varying_midpoint' else 1]
    if raw['duration_s'] != expected_time:
        raise ValueError('Snapshot does not match the archived Julia reference time')
    d_odd = (raw['tau_odd']-.5)*dx**2/(3*dt)
    d_magic = raw['magic']/(raw['tau_even']-.5)*dx**2/(3*dt)
    np.testing.assert_allclose([d_odd, d_magic], raw['diffusivity_m2_s'], rtol=1e-9, atol=0)
    np.testing.assert_array_equal(mask, old_mask)
    np.testing.assert_allclose(velocity, old_velocity, rtol=0, atol=1e-14)
    for key in ('halfway_links.json', 'halfway_link_fields.csv', 'stencil.csv'):
        if sha(directory / key) != sha(old / key):
            raise ValueError('Wall links or stencil changed: '+key)
    fluid = mask == 1
    np.testing.assert_allclose(fields['initial'][fluid], old_fields['initial'][fluid], rtol=0, atol=1e-12)
    for array in (*fields.values(), *(history[k] for k in history.dtype.names)):
        if not np.isfinite(array).all():
            raise ValueError('Non-finite field or history')
    np.testing.assert_array_equal(history['step'], np.arange(steps+1))
    np.testing.assert_allclose(history['time_s'], np.arange(steps+1)*dt, rtol=1e-14, atol=0)
    z_error = max(float(np.max(abs(f-f[:, :, :1])[fluid])) for f in fields.values())
    if z_error > 1e-12:
        raise ValueError('Axial invariance failed in liquid')
    lo, hi = history['min_concentration_kg_m3'], history['max_concentration_kg_m3']
    for field, index in [('initial', 0), ('middle', steps//2), ('peak', raw['peak_step']), ('final', steps)]:
        c = fields[field][fluid]
        np.testing.assert_allclose([c.min(), c.max()], [lo[index], hi[index]], rtol=0, atol=2e-14)
    amp = raw['amplitude_kg_m3']
    low, high = float(lo.min()), float(hi.max())
    excursion = max(.5-amp-low, high-(.5+amp), 0.)
    if abs(excursion-raw['maximum_excursion_kg_m3']) > 1e-14:
        raise ValueError('Peak excursion and history differ')
    liquid = history['fluid_inventory_kg']
    np.testing.assert_array_equal(liquid, history['kinetic_inventory_kg'])
    np.testing.assert_array_equal(history['active_wall_population_kg'], 0.)
    np.testing.assert_allclose(liquid[0], .5*dx**3*fluid.sum(), rtol=1e-11, atol=0)
    drift = float(np.max(abs(liquid/liquid[0]-1)))
    if drift > 1e-8:
        raise ValueError(f'Liquid inventory drift {drift:.3g} exceeds 1e-8')
    stats = dict(raw=raw, range_kg_m3=[low, high], maximum_excursion_kg_m3=excursion,
        peak_time_s=raw['peak_step']*dt, liquid_inventory_relative_drift=drift,
        axial_invariance_maximum_error=z_error, reconstructed_diffusivity_m2_s=d_magic,
        checks=dict(engineering_bounds=excursion<=EXCURSION_LIMIT, physical_bounds=excursion<=1e-8))
    if not amp:
        stats['checks']['constant'] = excursion <= 1e-9
    return stats, fluid, fields, history, coords


def analyze(base, specs, control, reference_path):
    constant, _, _, _, _ = verify_native(base, 'constant', specs['constant'], control)
    varying, fluid, fields, history, coords = verify_native(base, 'varying', specs['varying'], control)
    mid_stats, _, mid_fields, mid_history, _ = verify_native(base, 'varying_midpoint', specs['varying_midpoint'], control)
    np.testing.assert_array_equal(mid_history, history[:len(mid_history)])
    s, modes = ref.read_mode(reference_path)
    x, y = coords[0][:, None, None], coords[1][None, :, None]
    shell = fluid & (np.hypot(x-ref.L/2, y-ref.L/2) >= ref.R-SHELL_M)
    varying['fields'] = {}
    for label, time, field, mode in zip(('middle', 'final'), ref.TIMES,
            (mid_fields['final'], fields['final']), modes, strict=True):
        expected = np.broadcast_to(ref.profile(x, y, time, ref.A, (s, mode)), fluid.shape)
        error = (field-expected)/ref.A
        metrics = dict(time_s=time, RMS_over_amplitude=rms(error[fluid]),
            shell_RMS_over_amplitude=rms(error[shell]), interior_RMS_over_amplitude=rms(error[fluid & ~shell]),
            liquid_sd_kg_m3=float(np.std(field[fluid])), reference_sd_at_cells_kg_m3=float(np.std(expected[fluid])))
        varying['fields'][label] = metrics
        varying['checks'][label+'_profile'] = metrics['RMS_over_amplitude'] <= RMS_LIMIT
        varying['checks'][label+'_shell'] = metrics['shell_RMS_over_amplitude'] <= RMS_LIMIT
    return dict(constant=constant, varying=varying, varying_midpoint=mid_stats, midpoint_history_matches=True)


def run():
    control, reference_path, binary = control_record()
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid4().hex[:8]
    base = OUTPUT / 'runs' / stamp
    base.mkdir(parents=True)
    specs = jobs(control)
    files = [HERE / name for name in ('run.py', 'protocol.md')]
    hashes = source_manifest(files)
    for path in hashes:
        target = base / 'source_snapshot' / path
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT / path, target)
    record = dict(schema='pore2process.stationary-wall-time-step.v1', run_id=stamp, jobs=specs,
        source_sha256=hashes, control_summary=str(CONTROL.relative_to(ROOT)), control_summary_sha256=CONTROL_SHA,
        binary_sha256=control['binary_sha256'], reference=control['reference'],
        reference_path=str(reference_path.relative_to(ROOT)), new_bed_runs=0, new_DWSIM_runs=0)
    write_json(base / 'plan.json', record)
    archived_binary = base / 'binary/capillary3d'
    archived_binary.parent.mkdir()
    shutil.copy2(binary, archived_binary)
    print('Time-step comparison', base, flush=True)
    with ThreadPoolExecutor(max_workers=2) as pool:
        for future in [pool.submit(execute, base, name, spec, control) for name, spec in specs.items()]:
            future.result()
    record.update(analysis=analyze(base, specs, control, reference_path),
                  vendor=verify_vendor(), new_native_runs=len(specs))
    record['artifacts_sha256'] = {str(path.relative_to(base)): sha(path) for path in base.rglob('*') if path.is_file()}
    path = base / 'summary.json'
    write_json(path, record)
    write_json(OUTPUT / 'latest.json', dict(summary=str(path.relative_to(ROOT))))
    print('Completed time-step comparison', path, flush=True)


def check(path):
    path = path.resolve()
    saved = validate_record(path)
    control, mode, binary = control_record()
    specs = jobs(control)
    if (saved['control_summary_sha256'] != CONTROL_SHA or saved['reference'] != control['reference']
            or ROOT / saved['reference_path'] != mode or saved['jobs'] != specs
            or sha(path.parent / 'binary/capillary3d') != sha(binary)):
        raise ValueError('Archived comparison inputs differ')
    if analyze(path.parent, specs, control, mode) != saved['analysis']:
        raise ValueError('Replayed time-step comparison differs')
    result = dict(summary=str(path.relative_to(ROOT)), summary_sha256=sha(path), sources=len(saved['source_sha256']),
        artifacts=len(saved['artifacts_sha256']), native_runs=len(specs), new_native_runs_in_check=0, vendor=verify_vendor())
    write_json(OUTPUT / 'verification.json', result)
    print(result)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', type=Path)
    args = parser.parse_args()
    check(args.check) if args.check else run()
