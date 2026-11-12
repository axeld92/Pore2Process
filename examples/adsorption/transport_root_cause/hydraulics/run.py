"""Add real, independently converged native flow grids at unchanged physical Q."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import subprocess
from time import perf_counter
from uuid import uuid4

import numpy as np

from pore2process_demo.paths import ROOT
from pore2process_demo.io import sha, write_json
from pore2process_demo.provenance import source_manifest
from pore2process_demo.runtime.openlb import build_case, verify_vendor
from examples.adsorption.bed import analyze
from examples.adsorption.transport_accuracy.run import validate_record

HERE = Path(__file__).parent
OUTPUT = ROOT / 'results/adsorption/transport_root_cause/hydraulics'
BUILD = ROOT / '.local/adsorption/transport_root_cause_flow'


def latest_parent():
    pointer = ROOT / 'results/adsorption/column_transport/hydraulics/latest.json'
    return (ROOT / json.loads(pointer.read_text())['summary']).parent


def native_checks(analysis, target):
    r = analysis['raw']
    return dict(flow=abs(analysis['superficial_velocity_m_s']/target-1) <= 1e-4,
        mass=analysis['plane_mass_relative_spread'] <= 1e-6,
        stationarity=max(analysis['stationarity'].values()) <= 1.0001e-7,
        Mach=r['max_Mach'] < .05, density=r['max_density_departure'] < 1e-4,
        porosity=abs(analysis['voxel_porosity']-.6) < .015)


def source_point(native, carrier):
    a = analyze(native)
    if not all(native_checks(a, carrier['native_target_velocity_m_s']).values()):
        raise ValueError('Refined native flow has not passed its declared screens.')
    field = np.genfromtxt(native/'field.csv', delimiter=',', names=True)
    scale = carrier['carrier_viscosity_Pa_s']/carrier['carrier_density_kg_m3']/1e-6
    maxu = float(np.max(np.sqrt(sum(field[k]**2 for k in ('ux_m_s','uy_m_s','uz_m_s')))))*scale
    return dict(field=str((native/'field.csv').resolve()), velocity_scale=scale,
        n=a['raw']['resolution'], actual_superficial_velocity_m_s=a['superficial_velocity_m_s']*scale,
        actual_flow_m3_s=a['superficial_velocity_m_s']*scale*carrier['area_m2'],
        side_m=a['raw']['unit_cell_side_m'],
        speed_m_s=a['superficial_velocity_m_s']*scale/a['voxel_porosity'],
        max_speed_m_s=maxu, voxel_porosity=a['voxel_porosity'],
        field_sha256=sha(native/'field.csv'), raw_sha256=sha(native/'raw.json'))


def refine(base, binary, n, target, gradient):
    history = []
    previous = None
    warm = None
    warm_scale = 1.
    for iteration in range(1, 9):
        dest = base/f'native_flow/d800_grid{n}/iteration{iteration}'
        dest.mkdir(parents=True)
        command = [str(binary), str(n), repr(gradient), repr(.5+.0003*n), '.0008']
        if warm is not None:
            command += [str(warm), repr(warm_scale)]
        request = dict(command=command, binary_sha256=sha(binary),
            target_superficial_velocity_m_s=target,
            initial_field_sha256=None if warm is None else sha(warm),
            initial_state='rest' if warm is None else 'previous converged native same-grid field; estimate only',
            role='Native pressure root at fixed physical flow; not an adsorption calibration.')
        write_json(dest/'request.json', request)
        started = perf_counter()
        print('Native flow start', n, iteration, 'gradient', gradient, flush=True)
        with (dest/'stdout.log').open('w') as stream:
            result = subprocess.run(command, cwd=dest, stdout=stream,
                stderr=subprocess.STDOUT, timeout=10800)
        if result.returncode:
            raise RuntimeError('Native flow failed; inspect '+str(dest/'stdout.log'))
        a = analyze(dest)
        u = a['superficial_velocity_m_s']
        error = abs(u/target-1)
        checks = native_checks(a, target)
        print('Native flow finished', n, iteration, 'flow error', error,
            'elapsed_s', perf_counter()-started, 'checks', checks, flush=True)
        item = dict(directory=str(dest.relative_to(base)), analysis=a,
            flow_relative_error=error, checks=checks, initial_field_sha256=request['initial_field_sha256'])
        history.append(item)
        write_json(base/f'grid{n}_progress.json', dict(resolution=n, iterations=history))
        if error <= 1e-4:
            if not all(checks.values()):
                raise ValueError('Native converged fixed-flow point failed screens: '+str(checks))
            return dict(accepted=item, iterations=history, checks=checks)
        next_gradient = gradient*target/u
        if previous is not None:
            old_g, old_u = previous
            if u != old_u:
                candidate = gradient+(target-u)*(gradient-old_g)/(u-old_u)
                if .5*gradient < candidate < 2*gradient:
                    next_gradient = candidate
        previous = gradient, u
        gradient = next_gradient
        warm = dest/'field.csv'
        warm_scale = target/u
    raise RuntimeError('Fixed-flow pressure root did not converge in eight native runs.')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--parent', type=Path)
    parser.add_argument('--grid', type=int, choices=[64,96], default=64)
    parser.add_argument('--previous-refinement', type=Path)
    args = parser.parse_args()
    hydraulic = (args.parent or latest_parent()).resolve()
    parent = validate_record(hydraulic/'summary.json')
    if not parent['hydraulic_checks_passed']:
        raise ValueError('The actual hydraulic parent must be qualified.')
    carrier = json.loads((hydraulic/'inputs.json').read_text())
    previous = None
    if args.previous_refinement:
        previous = validate_record(args.previous_refinement.resolve()/'summary.json')
        if previous['hydraulic_summary_sha256'] != sha(hydraulic/'summary.json'):
            raise ValueError('Previous refinement uses another physical hydraulic parent.')
        gradient = previous['native']['accepted']['analysis']['raw']['gradient_Pa_m']
    else:
        gradient = parent['cases']['d800']['native']['48']['accepted']['analysis']['raw']['gradient_Pa_m']
    rid = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid4().hex[:8]
    base = OUTPUT/'runs'/rid
    base.mkdir(parents=True)
    files = [Path(__file__), HERE/'protocol.md', *list((HERE/'openlb').iterdir()),
        ROOT/'examples/adsorption/transport_root_cause/flux_audit/population_audit.h',
        ROOT/'examples/adsorption/bed.py', ROOT/'examples/adsorption/transport_accuracy/run.py',
        ROOT/'uv.lock']
    sources = source_manifest(files)
    for name in sources:
        dest = base/'source_snapshot'/name
        dest.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(ROOT/name,dest)
    write_json(base/'inputs.json',carrier)
    build = build_case(HERE/'openlb',BUILD,base/'build')
    binary = base/'binary/capillary3d'
    binary.parent.mkdir()
    shutil.copy2(BUILD/'capillary3d',binary)
    record = dict(schema='pore2process.transport-root-hydraulic-refinement.v1',run_id=rid,
        resolution=args.grid, source_sha256=sources,
        hydraulic_summary=str((hydraulic/'summary.json').relative_to(ROOT)),
        hydraulic_summary_sha256=sha(hydraulic/'summary.json'),
        previous_refinement_summary=None if previous is None else str((args.previous_refinement.resolve()/'summary.json').relative_to(ROOT)),
        previous_refinement_sha256=None if previous is None else sha(args.previous_refinement.resolve()/'summary.json'),
        inputs=carrier,build=build,binary_sha256=sha(binary),
        geometry_role='Synthetic FCC spheres, stock indicator, same 0.8 mm diameter and nominal 0.6 porosity.',
        new_DWSIM_runs=0,new_adsorption_predictions=0)
    write_json(base/'plan.json',record)
    print('Hydraulic refinement output',base,flush=True)
    group = refine(base,binary,args.grid,carrier['native_target_velocity_m_s'],gradient)
    accepted = base/group['accepted']['directory']
    point = source_point(accepted,carrier)
    write_json(base/'source_point.json',point)
    raw = group['accepted']['analysis']['raw']
    record.update(native=group,source_point=point,new_native_flow_runs=len(group['iterations']),
        geometry=dict(resolution=args.grid,fluid_nodes=raw['fluid_nodes'],
            side_m=raw['unit_cell_side_m'],dx_m=raw['dx_m'],
            domain_volume_m3=raw['unit_cell_side_m']**3,
            physical_liquid_volume_m3=raw['fluid_nodes']*raw['dx_m']**3,
            voxel_porosity=point['voxel_porosity'],particle_diameter_m=raw['particle_diameter_m'],
            coordinate_system='Cartesian metres, node centres at (index+0.5)*dx, zero-based periodic i,j,k.'),
        checks_passed=all(group['checks'].values()),vendor=verify_vendor())
    record['artifacts_sha256']={str(p.relative_to(base)):sha(p) for p in base.rglob('*') if p.is_file()}
    write_json(base/'summary.json',record)
    write_json(OUTPUT/f'grid{args.grid}_latest.json',dict(summary=str((base/'summary.json').relative_to(ROOT))))
    print('Hydraulic refinement accepted',base,point,flush=True)


if __name__ == '__main__':
    main()
