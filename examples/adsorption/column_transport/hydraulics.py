"""Actual fixed-flow native hydraulics and actual DWSIM; no tracer substitution."""
from datetime import datetime, timezone
from concurrent.futures import ThreadPoolExecutor
import argparse
import json
import math
from pathlib import Path
import shutil
import subprocess
from time import perf_counter, sleep
from uuid import uuid4

from pore2process_demo.paths import ROOT
from pore2process_demo.io import sha, write_json
from pore2process_demo.provenance import source_manifest
from pore2process_demo.runtime.openlb import build_case, verify_vendor
from examples.adsorption.bed import analyze
from examples.adsorption.dwsim import build_adapter, solve_case

HERE = Path(__file__).parent
BUILD = ROOT / '.local/adsorption/column_transport_flow'
OUTPUT = ROOT / 'results/adsorption/column_transport'


def fixed_flow(base, name, n, diameter, target, tau):
    """Adjust driving pressure to an imposed flow; every iteration is native."""
    parent = base / 'native_flow' / name
    parent.mkdir(parents=True,exist_ok=True)
    reference = json.loads((ROOT/'results/adsorption/bed/summary.json').read_text())['cases'][f'grid{n}']
    warm = ROOT/f'results/adsorption/bed/grid{n}/field.csv'
    previous_u = reference['superficial_velocity_m_s']
    gradient = 10*(.0023/diameter)**2*target/previous_u
    history = []
    previous = None
    for iteration in range(1, 9):
        dest = parent/f'iteration{iteration}'
        reused=(dest/'raw.json').exists()
        dest.mkdir(exist_ok=True)
        command = [str(BUILD/'capillary3d'), str(n), repr(gradient), repr(tau),
                   repr(diameter), str(warm), repr(target/previous_u)]
        request=dict(command=command, binary_sha256=sha(BUILD/'capillary3d'),
            initial_field_sha256=sha(warm), target_superficial_velocity_m_s=target,
            role='Native pressure root for imposed flow; not an adsorption fit.')
        started = perf_counter()
        if reused:
            saved=json.loads((dest/'request.json').read_text())
            # A resumed copy retains its true old argv. Its warm-file bytes
            # must match even when the enclosing run directory changed.
            actual=list(saved['command']);expected=list(command)
            actual[5]=expected[5]
            if actual!=expected or any(saved[k]!=request[k] for k in ('binary_sha256','initial_field_sha256','target_superficial_velocity_m_s')):
                raise ValueError('Resumed native request/binary/initial field differs: '+str(dest))
        else:
            if (dest/'request.json').exists():
                raise ValueError('Incomplete native request must remain in its original folder: '+str(dest))
            write_json(dest/'request.json',request)
            with (dest/'stdout.log').open('w') as stream:
                result = subprocess.run(command, cwd=dest, stdout=stream, stderr=subprocess.STDOUT, timeout=1800)
            if result.returncode:
                raise RuntimeError(f'Native flow failed: {dest}/stdout.log')
        a = analyze(dest)
        u = a['superficial_velocity_m_s']
        error = abs(u/target-1)
        print(name, iteration, 'flow error', f'{error:.3g}', 'reused' if reused else f'{perf_counter()-started:.1f}s', flush=True)
        history.append(dict(directory=str(dest.relative_to(base)), analysis=a, flow_relative_error=error,reused=reused))
        if error <= 1e-4:
            checks = dict(flow=error <= 1e-4, mass=a['plane_mass_relative_spread'] <= 1e-6,
                stationarity=max(a['stationarity'].values()) <= 1.0001e-7,
                Mach=a['raw']['max_Mach'] < .05, density=a['raw']['max_density_departure'] < 1e-4,
                porosity=abs(a['voxel_porosity']-.6) < .015)
            return dict(accepted=history[-1], iterations=history, checks=checks)
        # Secant pressure root, bounded to positive moderate changes. Warm start
        # uses the immediately previous actual velocity field as an estimate.
        next_gradient = gradient*target/u
        if previous is not None:
            old_g, old_u = previous
            if u != old_u:
                candidate = gradient+(target-u)*(gradient-old_g)/(u-old_u)
                if .5*gradient < candidate < 2*gradient:
                    next_gradient = candidate
        previous = gradient, u
        gradient = next_gradient
        warm, previous_u = dest/'field.csv', u
    raise RuntimeError(f'{name}: imposed flow did not converge in eight actual native runs.')


def pump_inputs(anchor, pressure, flow):
    return dict(Temperature_K=295.15, Inlet_pressure_Pa=1e5, Outlet_pressure_Pa=1e5+pressure,
        Mass_flow_kg_s=anchor['inlet_density_kg_m3']*flow, Pump_efficiency_fraction=.75,
        Scope='Synthetic FCC bed at source Figure 7 flow and particle diameter. Native operating-point hydraulic load; no qualified new adsorption/service prediction.')


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--resume-from',type=Path)
    args=parser.parse_args()
    start = perf_counter()
    rid = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid4().hex[:8]
    base = OUTPUT/'hydraulics/runs'/rid
    base.mkdir(parents=True)
    prior=None
    if args.resume_from:
        old=args.resume_from.resolve()
        if (old/'summary.json').exists():
            raise ValueError('A completed study is immutable; start a new study without resume.')
        requests=list((old/'native_flow').glob('*/iteration*/request.json'))
        # Preserve in-flight output in the old run. Its independent group
        # waits below while other groups progress; never copy an open history.
        completed=[p for p in requests if (p.parent/'raw.json').exists()]
        for p in completed:
            target=base/'native_flow'/p.parent.relative_to(old/'native_flow')
            target.parent.mkdir(parents=True,exist_ok=True)
            shutil.copytree(p.parent,target)
        shutil.copytree(old/'water_anchor',base/'water_anchor')
        shutil.copytree(old/'source_snapshot',base/'prior_source_snapshot')
        prior=dict(directory=str(old.relative_to(ROOT)),reused_requests_at_resume=len(completed),
            in_flight_at_resume=[str(p.parent.relative_to(old)) for p in requests if p not in completed],
            source_snapshot_sha256={str(p.relative_to(old/'source_snapshot')):sha(p)
                for p in (old/'source_snapshot').rglob('*') if p.is_file()},
            native_physics_unchanged=True,reason='Resume independent native groups in two subprocesses; no dynamics, geometry, discretization or thresholds changed.')
        write_json(base/'prior_execution.json',prior)
    files = [Path(__file__), HERE/'protocol.md', *list((HERE/'openlb_flow').iterdir()),
        ROOT/'examples/adsorption/bed.py', ROOT/'examples/adsorption/dwsim.py',
        *list((ROOT/'examples/adsorption/dwsim_runner').iterdir()),
        ROOT/'src/pore2process_demo/runtime/dwsim.py']
    sources = source_manifest(files)
    for name in sources:
        target = base/'source_snapshot'/name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT/name, target)
    build = build_case(HERE/'openlb_flow', BUILD, base/'build')
    adapter = build_adapter(base/'adapter_build')
    if prior:
        anchor=json.loads((base/'water_anchor/summary.json').read_text())
    else:
        anchor = solve_case(base/'water_anchor', dict(Temperature_K=295.15, Inlet_pressure_Pa=1e5,
            Outlet_pressure_Pa=2e5, Mass_flow_kg_s=.1, Pump_efficiency_fraction=.75,
            Scope='Actual DWSIM carrier-property anchor only; not an adsorption prediction.'))
    mu, rho = anchor['inlet_viscosity_Pa_s'], anchor['inlet_density_kg_m3']
    nu = mu/rho
    area = math.pi*.011**2/4
    flow = 43.8e-6/60
    physical_u = flow/area
    target = physical_u*1e-6/nu
    source_mass = 740*area*.1
    particle_density = 740/(1-.427)
    length = source_mass/(.4*particle_density*area)
    scenarios = dict(flow_m3_s=flow, area_m2=area, superficial_velocity_m_s=physical_u,
        native_target_velocity_m_s=target, native_kinematic_viscosity_m2_s=1e-6,
        carrier_viscosity_Pa_s=mu, carrier_density_kg_m3=rho,
        nominal_porosity=.6, particle_density_kg_m3=particle_density,
        dry_mass_kg=source_mass, length_m=length,
        source_printed_bulk_density_kg_m3=740., source_printed_voidage=.427,
        source_bed_length_m=.1, inventory_role='Source density-derived mass, transferred as a scenario; not weighed per-run mass.',
        geometry_role='Synthetic FCC at .6 voidage; not experimental specimen reconstruction.')
    write_json(base/'inputs.json', scenarios)
    specs=[]
    for label,diameter in [('d800',.0008),('d450',.00045)]:
        for n in (32,48):specs.append((label,str(n),n,diameter,.5+.0003*n))
        specs.append((label,'time48',48,diameter,.5+.00015*48))
    def execute(spec):
        label,key,n,diameter,tau=spec
        name=f'{label}_'+('time48' if key=='time48' else 'grid'+key)
        if prior:
            for request in (old/'native_flow'/name).glob('iteration*/request.json'):
                dest=base/'native_flow'/name/request.parent.name
                if dest.exists():continue
                wait_start=perf_counter()
                print('Waiting for original native child',name,request.parent.name,flush=True)
                while not (request.parent/'raw.json').exists():
                    if perf_counter()-wait_start>1800:raise TimeoutError('Original native child did not finish.')
                    sleep(.5)
                # Raw JSON is written last; give native iostream destructors
                # time to close before copying the immutable completed case.
                sleep(.25)
                analyze(request.parent)
                dest.parent.mkdir(parents=True,exist_ok=True)
                shutil.copytree(request.parent,dest)
        return label,key,fixed_flow(base,name,n,diameter,target,tau)
    # Independent native processes only. DWSIM stays on this main thread after
    # all six groups finish. Saved cases/locks/Julia modules are not mutated.
    with ThreadPoolExecutor(max_workers=2) as pool:
        results=list(pool.map(execute,specs))
    groups={(label,key):value for label,key,value in results}
    cases = {}
    for label, diameter in [('d800', .0008), ('d450', .00045)]:
        grid={key:groups[(label,key)] for key in ('32','48','time48')}
        a = grid['48']['accepted']['analysis']
        coefficient = a['raw']['gradient_Pa_m']/a['superficial_velocity_m_s']
        # Same geometry and Re; rho U^2 pressure scaling between native and
        # actual carrier. This is point resistance, not a Darcy extrapolation.
        native_gradient_at_target = coefficient*target
        physical_gradient = native_gradient_at_target*(rho/1000)*(nu/1e-6)**2
        pressure = physical_gradient*length
        pump = solve_case(base/label/'dwsim', pump_inputs(anchor, pressure, flow))
        refinement = dict(grid_pressure_relative_change=abs(coefficient /
            (grid['32']['accepted']['analysis']['raw']['gradient_Pa_m']/
             grid['32']['accepted']['analysis']['superficial_velocity_m_s'])-1),
            time_pressure_relative_change=abs((grid['time48']['accepted']['analysis']['raw']['gradient_Pa_m']/
             grid['time48']['accepted']['analysis']['superficial_velocity_m_s'])/coefficient-1))
        checks = dict(native=all(all(c['checks'].values()) for c in grid.values()),
            grid=refinement['grid_pressure_relative_change'] <= .05,
            time=refinement['time_pressure_relative_change'] <= .02,
            pump=pump['pump_power_relative_error'] <= 1e-8,
            pump_mass=abs(pump['mass_balance_error_kg_s']) <= 1e-10*pump['inlet_mass_flow_kg_s'],
            pump_energy=abs(pump['energy_balance_error_W']) <= .001*pump['shaft_power_W'],
            reload=pump['save_reload_solve_checked'],
            liquid=max(pump['inlet_vapor_fraction'],pump['outlet_vapor_fraction']) < 1e-10)
        cases[label] = dict(diameter_m=diameter, native=grid, refinement=refinement, checks=checks,
            particle_Reynolds=physical_u*diameter/nu, pressure_drop_Pa=pressure,
            native_pressure_rescaling_fraction=target/a['superficial_velocity_m_s']-1,
            operating_point_apparent_permeability_m2=a['permeability_m2'],
            pressure_gradient_Pa_m=physical_gradient, dwsim=pump,
            shaft_power_W=pump['shaft_power_W'],
            specific_pump_energy_kWh_m3=pump['shaft_power_W']/flow/3.6e6,
            breakthrough_time_s=None, service_volume_m3=None, cumulative_service_energy_kWh=None,
            unavailable_reason='Native solute transport has not passed application qualification.')
        write_json(base/f'{label}.json', cases[label])
    summary = dict(schema='pore2process.column-transport-hydraulics.v1', run_id=rid,
        source_sha256=sources, inputs=scenarios, cases=cases,
        hydraulic_checks_passed=all(all(c['checks'].values()) for c in cases.values()),
        new_native_flow_runs=sum(not r['reused'] for c in cases.values() for g in c['native'].values() for r in g['iterations']),
        reused_native_flow_runs=sum(r['reused'] for c in cases.values() for g in c['native'].values() for r in g['iterations']),
        new_DWSIM_runs=2 if prior else 3,reused_DWSIM_cases=1 if prior else 0,
        prior_execution=prior,new_adsorption_predictions=0, build=build, adapter=adapter,
        vendor=verify_vendor(), elapsed_s=perf_counter()-start,
        service_prediction_qualified=False)
    summary['artifacts_sha256'] = {str(p.relative_to(base)):sha(p) for p in base.rglob('*') if p.is_file()}
    write_json(base/'summary.json', summary)
    write_json(OUTPUT/'hydraulics/latest.json', dict(summary=str((base/'summary.json').relative_to(ROOT))))
    print('Hydraulics saved', base, summary['hydraulic_checks_passed'], flush=True)


if __name__ == '__main__':
    main()
