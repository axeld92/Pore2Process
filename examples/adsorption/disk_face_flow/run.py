"""Run or replay Julia transport with prescribed balanced face flows on the disk."""
import argparse
from datetime import datetime, timezone
from pathlib import Path
import shutil
import subprocess
from time import perf_counter
from uuid import uuid4

import numpy as np

from examples.adsorption.stationary_curved_wall import reference as ref
from examples.adsorption.stock_transport_candidate.run import PARENT, PARENT_SHA, parent_record
from examples.adsorption.transport_accuracy.run import validate_record
from examples.adsorption.transport_verification.analysis import EXCURSION_LIMIT, RMS_LIMIT, SHELL_M, load
from examples.adsorption.transport_verification.reference import ROOT_J2_PRIME, bessel
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT
from pore2process_demo.provenance import source_manifest

HERE = Path(__file__).parent
OUTPUT = ROOT / 'results/adsorption/disk_face_flow'
SPACE_LIMIT = .005
TIME_LIMIT = .001
DEPTH = 4*ref.L/128
CHECK_TIMES = ('middle', 'final')


def spec(nr, na, dt, amplitude=ref.A):
    return dict(nr=nr, na=na, dt_max_s=dt, amplitude_kg_m3=amplitude)


def spatial_grid(case):
    nr, na = case['nr'], case['na']
    edges = np.linspace(0, ref.R, nr+1)
    theta = (np.arange(na)+.5)*(2*np.pi/na)
    return edges, theta


def cell_average(edges, theta, radial_mode, amplitude):
    # Two-point Gauss quadrature in r², where equal weights represent equal areas.
    a2, b2 = edges[:-1]**2, edges[1:]**2
    center, offset = (a2+b2)/2, (b2-a2)/(2*np.sqrt(3.))
    mean = (radial_mode(np.sqrt(center-offset))+radial_mode(np.sqrt(center+offset)))/2
    return .5+amplitude*np.sinc(2/len(theta))*np.real(mean[:, None]*np.exp(2j*theta))


def initial(case):
    edges, theta = spatial_grid(case)
    mode = lambda r: bessel(2, ROOT_J2_PRIME*r/ref.R)/bessel(2, ROOT_J2_PRIME)
    return cell_average(edges, theta, mode, case['amplitude_kg_m3'])


def sample(field, edges, theta, x, y):
    """Interpolate sector averages to the archived Cartesian sample centres."""
    radius = np.hypot(x, y)
    angle = np.arctan2(y, x) % (2*np.pi)
    # The origin has no preferred angle and the smooth wall has zero radial derivative.
    extended = np.vstack((np.full(len(theta), field[0].mean()), field, field[-1]))
    centers = np.r_[0., (edges[1:]+edges[:-1])/2, edges[-1]]
    ri = np.searchsorted(centers, radius, side='right')-1
    rf = (radius-centers[ri])/(centers[ri+1]-centers[ri])
    azimuth = angle/(2*np.pi/len(theta))-.5
    j = np.floor(azimuth).astype(int)
    af = azimuth-j
    lower = (1-af)*extended[ri, j % len(theta)]+af*extended[ri, (j+1) % len(theta)]
    upper = (1-af)*extended[ri+1, j % len(theta)]+af*extended[ri+1, (j+1) % len(theta)]
    return (1-rf)*lower+rf*upper


def audit_grid(rings, edges, theta):
    np.testing.assert_array_equal(rings['ring'], np.arange(1, len(edges)))
    np.testing.assert_allclose(rings['r_inner_m'], edges[:-1], rtol=1e-14, atol=0)
    np.testing.assert_allclose(rings['r_outer_m'], edges[1:], rtol=1e-14, atol=0)
    center = (edges[1:]+edges[:-1])/2
    dr = ref.R/(len(edges)-1)
    da = 2*np.pi/len(theta)
    np.testing.assert_allclose(rings['r_center_m'], center, rtol=1e-14, atol=0)
    volume = rings['cell_volume_m3']
    expected_volume = DEPTH*da*np.diff(edges**2)/2
    np.testing.assert_allclose(volume, expected_volume, rtol=1e-12, atol=0)
    np.testing.assert_allclose(rings['theta_diffusion_m3_s'], ref.D*DEPTH*dr/(center*da), rtol=1e-14, atol=0)
    np.testing.assert_allclose(rings['outer_radial_diffusion_m3_s'][:-1],
                               ref.D*DEPTH*da*edges[1:-1]/dr, rtol=1e-14, atol=0)
    q = rings['theta_flow_m3_s']
    # Two-point radial Gauss quadrature integrates the cubic velocity exactly.
    r1, r2 = center-dr/(2*np.sqrt(3)), center+dr/(2*np.sqrt(3))
    speed = lambda r: ref.OMEGA0*r*(1-(r/ref.R)**2)
    np.testing.assert_allclose(q, DEPTH*dr*(speed(r1)+speed(r2))/2, rtol=1e-12, atol=0)
    radial = rings['outer_radial_flow_m3_s']
    np.testing.assert_array_equal(radial, 0.)
    if rings['outer_radial_diffusion_m3_s'][-1] != 0:
        raise ValueError('Nonzero wall tracer conductance')
    east = np.broadcast_to(q[:, None], (len(q), len(theta)))
    net = east-np.roll(east, 1, axis=1)+(radial-np.r_[0., radial[:-1]])[:, None]
    maximum_net = float(np.max(abs(net)))
    if maximum_net > 1e-12*max(abs(q)):
        raise ValueError('Unbalanced cell carrier flow')
    return volume, dict(liquid_volume_m3=float(volume.sum()*len(theta)),
                        maximum_cell_net_flow_m3_s=maximum_net,
                        maximum_theta_face_flow_m3_s=float(max(abs(q))),
                        wall_carrier_flow_m3_s=float(radial[-1]),
                        wall_tracer_flux_kg_s=0.0)


def evaluate(directory, case, mode, samples):
    edges, theta = spatial_grid(case)
    nr, na = case['nr'], case['na']
    amp = case['amplitude_kg_m3']
    fields = {name: np.fromfile(directory / (name+'.f64'), dtype='<f8').reshape((nr, na), order='F')
              for name in ('initial', *CHECK_TIMES)}
    history = np.genfromtxt(directory / 'history.csv', delimiter=',', names=True)
    rings = np.genfromtxt(directory / 'rings.csv', delimiter=',', names=True)
    volume, grid = audit_grid(rings, edges, theta)
    for array in (*fields.values(), *(history[k] for k in history.dtype.names)):
        if not np.isfinite(array).all():
            raise ValueError('Non-finite saved concentration or history')
    np.testing.assert_array_equal(fields['initial'], initial(case))
    times = history['time_s']
    start, time_index = 0., 0
    indices = [0]
    effective_dt = []
    for end in ref.TIMES:
        count = int(np.ceil((end-start)/case['dt_max_s']))
        dt = (end-start)/count
        np.testing.assert_allclose(times[time_index+1:time_index+count+1],
                                   start+np.arange(1, count+1)*dt, rtol=1e-14, atol=0)
        time_index += count
        if times[time_index] != end:
            raise ValueError('Observation does not match the reference time')
        indices.append(time_index)
        effective_dt.append(dt)
        start = end
    if len(times) != time_index+1 or times[0] != 0.:
        raise ValueError('Incorrect integration history length')
    liquid = history['liquid_inventory_kg']
    low, high = history['min_concentration_kg_m3'], history['max_concentration_kg_m3']
    for (name, field), index in zip(fields.items(), indices, strict=True):
        np.testing.assert_allclose([field.min(), field.max()], [low[index], high[index]], rtol=0, atol=1e-14)
        np.testing.assert_allclose(np.sum(field*volume[:, None]), liquid[index], rtol=1e-12, atol=0,
                                   err_msg='Wrong liquid inventory for '+name)
    minimum, maximum = float(low.min()), float(high.max())
    excursion = max(.5-amp-minimum, maximum-(.5+amp), 0.)
    drift = float(np.max(abs(liquid/liquid[0]-1)))
    checks = dict(liquid_inventory=drift <= 1e-8, engineering_bounds=excursion <= EXCURSION_LIMIT,
                  physical_bounds=excursion <= 1e-8)
    metrics = dict(grid=grid, history_steps=time_index, effective_dt_s=effective_dt,
                   range_kg_m3=[minimum,maximum], maximum_excursion_kg_m3=excursion,
                   liquid_inventory_relative_drift=drift, checks=checks)
    if amp == 0:
        checks['uniform'] = excursion <= 1e-9
        return metrics, {}
    x, y = samples
    shell = np.hypot(x, y) >= ref.R-SHELL_M
    ring_shell = edges[:-1] >= ref.R-SHELL_M-1e-15
    s, modes = mode
    fixed = {}
    metrics['profiles'] = {}
    for name, time, radial in zip(CHECK_TIMES, ref.TIMES, modes, strict=True):
        exact = cell_average(edges, theta, lambda r: ref.sample(s, radial, r/ref.R), amp)
        err2 = np.mean(((fields[name]-exact)/amp)**2, axis=1)
        f = sample(fields[name], edges, theta, x, y)
        fixed[name] = f
        oracle = ref.profile(x+ref.L/2, y+ref.L/2, time, amp, (s, radial))
        fixed_err = (f-oracle)/amp
        errors = dict(sectors_global=float(np.sqrt(np.average(err2, weights=volume))),
                      sectors_shell=float(np.sqrt(np.average(err2[ring_shell], weights=volume[ring_shell]))),
                      fixed_global=float(np.sqrt(np.mean(fixed_err**2))),
                      fixed_shell=float(np.sqrt(np.mean(fixed_err[shell]**2))))
        metrics['profiles'][name] = dict(time_s=time, **errors)
        for key, value in errors.items():
            checks[name+'_'+key] = value <= RMS_LIMIT
    return metrics, fixed


def refinement(before, after, samples, limit):
    shell = np.hypot(*samples) >= ref.R-SHELL_M
    errors = {}
    for name in CHECK_TIMES:
        d = (after[name]-before[name])/ref.A
        errors[name] = dict(global_rms=float(np.sqrt(np.mean(d*d))), shell_rms=float(np.sqrt(np.mean(d[shell]**2))))
    return dict(errors=errors, limit=limit, passed=all(x <= limit for e in errors.values() for x in e.values()))


def passes(metrics):
    return all(value for name, value in metrics['checks'].items() if name != 'physical_bounds')


def archived():
    old, mode_path = parent_record()
    _, mask, _, _, _, q = load(PARENT.parent / 'native/shear128')
    x, y = np.meshgrid(q[0]-ref.L/2, q[1]-ref.L/2, indexing='ij')
    fluid = mask[:, :, 0] == 1
    return old['reference'], mode_path, (x[fluid], y[fluid])


def julia_command(directory, case, runtime):
    args = [case['nr'], case['na'], ref.R, DEPTH, ref.OMEGA0, ref.D, case['dt_max_s'], *ref.TIMES]
    return [*runtime, str(HERE / 'solve.jl'), str(directory), *map(str, args)]


def execute(directory, case, runtime):
    directory.mkdir(parents=True)
    cmd = julia_command(directory, case, runtime)
    write_json(directory / 'request.json', dict(command=cmd, case=case))
    initial(case).astype('<f8').ravel(order='F').tofile(directory / 'initial.f64')
    print('Julia start', directory.name, case, flush=True)
    start = perf_counter()
    with (directory / 'stdout.log').open('w') as log:
        process = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, timeout=900)
    write_json(directory / 'execution.json', dict(returncode=process.returncode, wall_s=perf_counter()-start))
    if process.returncode:
        raise RuntimeError(f'Julia failed. See {directory / "stdout.log"}')


def campaign(base, mode, samples, runtime, *, replay):
    records, fixed = {}, {}

    def case(name, settings):
        directory = base / 'cases' / name
        if not replay:
            execute(directory, settings, runtime)
        req = read_json(directory / 'request.json')
        if (req != dict(command=julia_command(directory, settings, runtime), case=settings)
                or read_json(directory / 'execution.json')['returncode']):
            raise ValueError('Execution does not match the planned case')
        metrics, fixed[name] = evaluate(directory, settings, mode, samples)
        records[name] = dict(settings=settings, metrics=metrics)
        if not replay:
            print('Julia complete', name, metrics['range_kg_m3'], metrics.get('profiles', {}), flush=True)

    case('coarse', spec(40, 256, 2e-4))
    case('medium', spec(80, 512, 2e-4))
    case('medium_half_dt', spec(80, 512, 1e-4))
    selected, same_mesh, coarser = 'medium_half_dt', 'medium', 'coarse'
    space = refinement(fixed[coarser], fixed[same_mesh], samples, SPACE_LIMIT)
    refinements = {'coarse_to_medium': space}
    if not passes(records[selected]['metrics']) or not space['passed']:
        case('fine', spec(160, 1024, 2e-4))
        case('fine_half_dt', spec(160, 1024, 1e-4))
        selected, same_mesh, coarser = 'fine_half_dt', 'fine', 'medium'
        space = refinement(fixed[coarser], fixed[same_mesh], samples, SPACE_LIMIT)
        refinements['medium_to_fine'] = space
    time = refinement(fixed[same_mesh], fixed[selected], samples, TIME_LIMIT)
    case('uniform', dict(records[selected]['settings'], amplitude_kg_m3=0.))
    qualified = (passes(records['uniform']['metrics']) and passes(records[selected]['metrics'])
                 and space['passed'] and time['passed'])
    return dict(cases=records, selected=selected, space=refinements, time=time, disk_qualified=qualified)


def run():
    qualification, mode_path, samples = archived()
    runtime = ref.runtime()
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid4().hex[:8]
    base = OUTPUT / 'runs' / stamp
    base.mkdir(parents=True)
    files = [HERE / name for name in ('run.py', 'solve.jl', 'protocol.md')]
    files += [ROOT / 'julia/adsorption/DiskFaceTransport.jl', ROOT / 'julia/Project.toml',
              ROOT / 'julia/Manifest.toml']
    hashes = source_manifest(files)
    for path in hashes:
        target = base / 'source_snapshot' / path
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT / path, target)
    shutil.copy2(mode_path, base / 'reference.csv')
    record = dict(schema='pore2process.disk-face-flow.v1', run_id=stamp, source_sha256=hashes,
        parent_summary=str(PARENT.relative_to(ROOT)), parent_summary_sha256=PARENT_SHA,
        reference=qualification, julia_command=runtime, new_bed_runs=0, new_openlb_runs=0, new_DWSIM_runs=0)
    write_json(base / 'plan.json', record)
    print('Balanced disk study', base, flush=True)
    record['analysis'] = campaign(base, ref.read_mode(base / 'reference.csv'), samples, runtime, replay=False)
    record['artifacts_sha256'] = {str(p.relative_to(base)): sha(p) for p in base.rglob('*') if p.is_file()}
    path = base / 'summary.json'
    write_json(path, record)
    write_json(OUTPUT / 'latest.json', dict(summary=str(path.relative_to(ROOT))))
    print('Completed balanced disk study', path, 'qualified', record['analysis']['disk_qualified'], flush=True)


def check(path):
    path = path.resolve()
    saved = validate_record(path)
    qualification, _, samples = archived()
    mode = path.parent / 'reference.csv'
    if (saved['parent_summary_sha256'] != PARENT_SHA or saved['reference'] != qualification
            or sha(mode) != qualification['reference_sha256']):
        raise ValueError('Comparison reference differs')
    replayed = campaign(path.parent, ref.read_mode(mode), samples, saved['julia_command'], replay=True)
    if replayed != saved['analysis']:
        raise ValueError('Replayed analysis differs')
    result = dict(summary=str(path.relative_to(ROOT)), summary_sha256=sha(path), sources=len(saved['source_sha256']),
                  artifacts=len(saved['artifacts_sha256']), disk_qualified=replayed['disk_qualified'])
    write_json(OUTPUT / 'verification.json', result)
    print(result)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', type=Path)
    args = parser.parse_args()
    check(args.check) if args.check else run()
