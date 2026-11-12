"""Independent read-only structural and physical checks of half-way controls.

This file never configures a boundary, mutates populations or launches native
jobs. Its topology oracle uses the exported stock descriptor and geometry.
"""
from itertools import product
import argparse
import json
from pathlib import Path
import numpy as np

from pore2process_demo.paths import ROOT
from pore2process_demo.io import sha,write_json
from examples.adsorption.transport_accuracy.run import validate_record


def descriptor(table):
    table=np.atleast_1d(table)
    if len(table)!=27 or set(int(x) for x in table['p'])!=set(range(27)):
        raise ValueError('Stock D3Q27 export must contain each population ID once.')
    table=table[np.argsort(table['p'])]
    directions=np.column_stack([table[name] for name in ('cx','cy','cz')])
    if not np.all(np.isfinite(directions)) or not np.all(directions==np.round(directions)):
        raise ValueError('Descriptor directions must be finite integers.')
    directions=directions.astype(int);opposite=np.asarray(table['opposite'],dtype=int)
    if set(map(tuple,directions))!=set(product((-1,0,1),repeat=3)) or np.any(directions[0]):
        raise ValueError('Export differs from the declared D3Q27 direction contract.')
    if np.any(opposite<0) or np.any(opposite>=27) or not np.array_equal(opposite[opposite],np.arange(27)):
        raise ValueError('Opposite population map is not an involution.')
    if not np.array_equal(directions[opposite],-directions):
        raise ValueError('Opposite populations do not have opposite directions.')
    weight=np.asarray(table['weight'],dtype=float)
    if not np.all(np.isfinite(weight)) or np.any(weight<=0) or abs(float(np.sum(weight))-1)>1e-13:
        raise ValueError('Stock weights must be positive and normalized.')
    np.testing.assert_allclose(weight[opposite],weight,atol=1e-16,rtol=1e-14)
    return directions,opposite,weight


def mask(table):
    table=np.atleast_1d(table)
    if not len(table):raise ValueError('Empty native material export.')
    indices=np.column_stack([table[name] for name in ('i','j','k')])
    if np.any(indices<0) or not np.all(indices==np.round(indices)):
        raise ValueError('Native material indices must be nonnegative integers.')
    indices=indices.astype(int);shape=tuple(np.max(indices,axis=0)+1)
    if len(table)!=np.prod(shape) or len(set(map(tuple,indices)))!=len(table):
        raise ValueError('Native mask must contain exactly one row for every core cell.')
    materials=np.asarray(table['material'],dtype=int)
    if set(materials)!={1,2}:raise ValueError('This candidate requires declared fluid1 and solid2.')
    grid=np.zeros(shape,dtype=int);grid[tuple(indices.T)]=materials
    return grid


def expected_links(materials,directions,periodic):
    """Independently enumerate fluid→solid links, including periodic seams."""
    materials=np.asarray(materials);shape=np.array(materials.shape)
    if materials.ndim!=3 or len(periodic)!=3:raise ValueError('Three-dimensional mask/periodicity required.')
    expected=set();seams=set()
    for point in np.argwhere(materials==1):
        for p,direction in enumerate(directions):
            if p==0:continue
            neighbor=point+direction;cross=(neighbor<0)|(neighbor>=shape)
            if np.any(cross & ~np.asarray(periodic,dtype=bool)):
                raise ValueError('Fluid stencil reaches an undeclared exterior boundary.')
            neighbor%=shape
            if materials[tuple(neighbor)]==2:
                key=(*map(int,point),int(p));expected.add(key)
                if np.any(cross):seams.add(key)
    return expected,seams


def coverage(materials,directions,links,periodic):
    expected,seams=expected_links(materials,directions,periodic)
    links=np.atleast_1d(links);keys=[tuple(int(row[k]) for k in ('i','j','k','p')) for row in links]
    if len(keys)!=len(set(keys)):raise ValueError('Duplicate active boundary links.')
    active=set(keys);wrong_q=int(np.count_nonzero(np.asarray(links['q'])!=.5))
    missing=expected-active;unexpected=active-expected
    return dict(expected_fluid_solid_links=len(expected),exported_active_links=len(active),
                expected_periodic_seam_links=len(seams),covered_periodic_seam_links=len(seams&active),
                missing_links=len(missing),unexpected_links=len(unexpected),non_halfway_links=wrong_q,
                first_missing=next(iter(sorted(missing)),None),first_unexpected=next(iter(sorted(unexpected)),None),
                passed=bool(expected) and not missing and not unexpected and wrong_q==0)


def structural_directory(directory,periodic):
    directory=Path(directory)
    d=descriptor(np.genfromtxt(directory/'stencil.csv',delimiter=',',names=True))
    m=mask(np.genfromtxt(directory/'mask.csv',delimiter=',',names=True))
    links=np.genfromtxt(directory/'halfway_link_fields.csv',delimiter=',',names=True)
    result=coverage(m,d[0],links,periodic)
    native=json.loads((directory/'halfway_links.json').read_text())
    result['matches_native_counts']=(result['expected_fluid_solid_links']==native['fluid_solid_links'] and
        result['expected_periodic_seam_links']==native['periodic_seam_links'])
    result['inactive_entries_unchanged_native_check']=native['inactive_changed']==0
    result['passed']=result['passed'] and result['matches_native_counts'] and result['inactive_entries_unchanged_native_check']
    return result


def physical_directory(directory,kind):
    """Bind the actual candidate fields and commands to the declared problem."""
    directory=Path(directory);raw=json.loads((directory/'raw.json').read_text())
    request=json.loads((directory/'request.json').read_text())
    geometry=mask(np.genfromtxt(directory/'mask.csv',delimiter=',',names=True))
    initial=np.genfromtxt(directory/'initial.csv',delimiter=',',names=True)
    final=np.genfromtxt(directory/'final.csv',delimiter=',',names=True)
    history=np.genfromtxt(directory/'inventory.csv',delimiter=',',names=True)
    for table in (initial,final,history):
        if not all(bool(np.all(np.isfinite(table[name]))) for name in table.dtype.names):
            raise ValueError('Non-finite native field/history')
    h=raw['height_m']/raw['n'] if kind=='flat' else raw['dx_m']
    count=int(np.count_nonzero(geometry==1));volume=count*h**3
    inventory=history['fluid_inventory_kg'] if kind=='flat' else history['fluid_tracer_kg']
    reference_volume_error=abs(float(inventory[0])/(.5*volume)-1)
    reconstructed_D=(raw['tau_odd']-.5)*h*h/(3*raw['dt_s'])
    if abs(reconstructed_D/1e-9-1)>1e-9:raise ValueError('Candidate diffusivity changed')
    if reference_volume_error>1e-11:raise ValueError('Fluid inventory does not match exported geometry')
    fluid=initial if kind=='flat' else initial[initial['material']==1]
    amplitude=raw['initial_amplitude_kg_m3'] if kind=='flat' else .25
    expected=.5+amplitude*np.sin(2*np.pi*fluid['x_m']/raw['length_m'])
    initial_error=float(np.max(abs(fluid['concentration_kg_m3']-expected)))
    if initial_error>1e-12:raise ValueError('Initial field differs from the declared sinusoid/constant')
    result=dict(fluid_nodes=count,fluid_volume_m3=volume,
        initial_inventory_relative_volume_error=reference_volume_error,
        maximum_initial_profile_error_kg_m3=initial_error,
        physical_diffusivity_m2_s=reconstructed_D,
        maximum_relative_fluid_inventory_drift=float(np.max(abs(inventory/inventory[0]-1))),
        actual_duration_s=float(history['time_s'][-1]),
        maximum_local_balance_residual=max(json.loads(p.read_text())['max_local_fluid_balance_residual']
            for p in directory.glob('halfway_audit_*.json')),
        maximum_actual_reflection_residual=max(json.loads(p.read_text())['max_reflection_pair_residual']
            for p in directory.glob('halfway_audit_*.json')))
    if result['maximum_relative_fluid_inventory_drift']>1e-8:raise ValueError('Physical fluid mass failed')
    if max(result['maximum_local_balance_residual'],result['maximum_actual_reflection_residual'])>1e-11:
        raise ValueError('Actual local balance/reflection failed')
    if kind=='porous':
        parent=ROOT/request['parent'];parent_path=parent.parent.parent/'summary.json'
        parent_summary=validate_record(parent_path)
        old_request=json.loads((parent/'request.json').read_text())
        if request['command'][1:]!=old_request['command'][1:]:raise ValueError('Porous physical command changed')
        old_raw=json.loads((parent/'raw.json').read_text())
        for name in ('n','nx','cells','side_m','dx_m','dt_s','diffusivity_m2_s','velocity_scale','periodic','steps'):
            if raw[name]!=old_raw[name]:raise ValueError('Porous physical input changed: '+name)
        old_final=np.genfromtxt(parent/'final.csv',delimiter=',',names=True)
        np.testing.assert_array_equal(final[['i','j','k','material']],old_final[['i','j','k','material']])
        for name in ('x_m','y_m','z_m','ux_m_s','uy_m_s','uz_m_s'):
            np.testing.assert_array_equal(final[name],old_final[name])
        native=np.genfromtxt(Path(request['command'][7]),delimiter=',',names=True)
        selected=final['material']==1
        np.testing.assert_array_equal(native[['i','j','k']],final[selected][['i','j','k']])
        coordinate_error={}
        for axis in ('x','y','z'):
            coordinate_error[axis]=float(np.max(abs(native[axis+'_m']-final[selected][axis+'_m'])))
            np.testing.assert_allclose(native['u'+axis+'_m_s']*raw['velocity_scale'],
                final[selected]['u'+axis+'_m_s'],atol=1e-16,rtol=1e-11)
        if max(coordinate_error.values())>1e-15:raise ValueError('Water/scalar coordinates differ')
        result.update(parent_summary=str(parent_path.relative_to(ROOT)),parent_summary_sha256=sha(parent_path),
            parent_case_recorded=parent.name in parent_summary['cases'],
            unchanged_physical_argument_bytes=True,unchanged_geometry_and_velocity=True,
            native_field_sha256=sha(Path(request['command'][7])),
            water_scalar_maximum_coordinate_error_m=coordinate_error)
    else:
        from ..flat_channel.run import PHYSICAL
        physical=PHYSICAL
        np.testing.assert_allclose([raw['height_m'],raw['dt_s'],raw['diffusivity_m2_s']],
            [physical['height_m'],physical['dt_s'],physical['diffusivity_m2_s']],rtol=1e-14,atol=0)
        if raw['steps']!=physical['steps']:raise ValueError('Flat observation window changed')
        if raw['mean_speed_m_s'] not in (0.,physical['mean_speed_m_s']):raise ValueError('Flat flow changed')
        np.testing.assert_allclose(final['x_m'],(final['i']+.5)*h,atol=1e-16,rtol=1e-13)
        np.testing.assert_allclose(final['y_m'],(final['j']+.5)*h,atol=1e-16,rtol=1e-13)
        velocity=6*raw['mean_speed_m_s']*(final['y_m']/raw['height_m'])*(1-final['y_m']/raw['height_m'])
        np.testing.assert_allclose(final['ux_m_s'],velocity,atol=1e-14,rtol=1e-12)
        result['unchanged_flat_physical_problem']=True
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--summary',type=Path,required=True)
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args();path=args.summary.resolve();source=validate_record(path)
    cases={}
    for name in source['cases']:
        kind='porous' if name.startswith('porous') else 'flat';directory=path.parent/'native'/name
        structure=structural_directory(directory,(True,kind=='porous',True))
        if not structure['passed']:raise ValueError('Independent link coverage failed: '+name)
        cases[name]=dict(structure=structure,physics=physical_directory(directory,kind))
        print('Independent review:',name,'passed structural/input/mass checks',flush=True)
    dependencies={c['physics']['parent_summary']:c['physics']['parent_summary_sha256']
        for c in cases.values() if 'parent_summary' in c['physics']}
    review_sources=[Path(__file__),ROOT/'tests/adsorption/test_transport_half_way_review.py']
    record=dict(schema='pore2process.independent-halfway-review.v1',
        candidate_summary=str(path.relative_to(ROOT)),candidate_summary_sha256=sha(path),
        parent_summary=str(path.relative_to(ROOT)),parent_summary_sha256=sha(path),
        source_sha256={str(p.relative_to(ROOT)):sha(p) for p in review_sources},artifacts_sha256={},
        extra_summaries=[dict(path=p,sha256=d) for p,d in sorted(dependencies.items())],
        cases=cases,new_native_runs=0,new_DWSIM_runs=0,
        original_particle_transport_qualified=False,
        scope='Independent topology, input, initial-field and inventory review. Does not waive failed profile/bound screens or near-wall errors.')
    output=path.parent/'independent_review.json'
    if args.check:
        if json.loads(output.read_text())!=record:raise ValueError('Independent review does not replay')
        print('Independent review replays; zero new native/DWSIM runs',flush=True)
    else:write_json(output,record)


if __name__=='__main__':main()
