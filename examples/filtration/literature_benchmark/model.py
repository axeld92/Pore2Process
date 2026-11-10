"""Source snapshots and calls to the actual Julia numerical package."""
from examples.filtration.calibration.audit_crossflow_sources import audit
from examples.filtration.graded_support.model import load_case


def record(jl, value):
    result = {}
    for field in jl.fieldnames(jl.typeof(value)):
        name = str(field)
        v = getattr(value, name)
        result[name] = (None if v is None else str(v) if name == 'classification'
                        else bool(v) if name == 'nominal_admissible'
                        else int(v) if name.endswith('_index') else float(v))
    return result


def prepare(jl, extraction, protocol):
    """Materialize source readings; no second numerical solver in Python."""
    readings = audit(extraction)['figure_readings']['chiu2005_fig3']
    vector = jl.seval('Vector{Float64}')
    arrays = [vector([p[key]['reading_interval'][side] / scale for p in readings])
              for key, side, scale in [('x', 0, 1), ('x', 1, 1),
                                       ('critical_flux_LMH', 0, 3.6e6),
                                       ('critical_flux_LMH', 1, 3.6e6)]]
    table = extraction['garcia2008_table1']
    decisions = []
    for row in table['rows']:
        for margin in protocol['capacity_fractions']:
            point = record(jl, jl.crossflow_benchmark_point(
                row['velocity_m_s'], row['critical_flux_LMH'] / 3.6e6,
                row['pressure_loss_Pa'], row['flow_L_h'] / 3.6e6,
                table['filtration_area_m2'], protocol['process_area_m2'],
                protocol['target_flux_LMH'] / 3.6e6, 1., margin))
            bracket = record(jl, jl.critical_flux_bracket(*arrays, row['velocity_m_s'],
                protocol['target_flux_LMH'] / 3.6e6, margin))
            decisions.append({'row_id': row['row_id'], 'capacity_fraction': margin,
                'velocity_m_s': row['velocity_m_s'], 'table': point,
                'figure_monotonic_bracket': bracket,
                'cross_source_admissible': point['nominal_admissible'] and bracket['classification'] == 'admissible'})
    return decisions


def requests(decisions, protocol, source):
    for convention in protocol['flow_conventions']:
        for d in decisions:
            if d['capacity_fraction'] != protocol['nominal_capacity_fraction'] or not d['table']['nominal_admissible']:
                continue
            for support in source['supports']:
                yield {'case_id': f"{convention['id']}/row_{d['row_id']:02d}/{support['name']}",
                       'row_id': d['row_id'], 'support': support['name'],
                       'flow_convention': convention['id'], 'flow_factor': convention['factor'],
                       'validation_only': False, 'seed_options': {}}
    # Strong seed/tolerance perturbation at the highest module/permeate ratio.
    high = protocol['flow_conventions'][-1]
    row = max(d['row_id'] for d in decisions if d['table']['nominal_admissible'])
    for label, options in [('high_seed', {'seed_multiplier': .1, 'seed_temperature_offset': 2.}),
                           ('high_tolerance', {'seed_multiplier': 2., 'tolerance_multiplier': .1})]:
        yield {'case_id': 'validation/' + label, 'row_id': row, 'support': 'straight',
               'flow_convention': high['id'], 'flow_factor': high['factor'],
               'validation_only': True, 'seed_options': options}


def translate(jl, extraction, protocol, hydraulics, source, request):
    table = extraction['garcia2008_table1']
    row = next(r for r in table['rows'] if r['row_id'] == request['row_id'])
    point = record(jl, jl.crossflow_benchmark_point(
        row['velocity_m_s'], row['critical_flux_LMH'] / 3.6e6, row['pressure_loss_Pa'],
        row['flow_L_h'] / 3.6e6, table['filtration_area_m2'], protocol['process_area_m2'],
        protocol['target_flux_LMH'] / 3.6e6, request['flow_factor'], protocol['nominal_capacity_fraction']))
    if not point['nominal_admissible']:
        raise ValueError('Inadmissible source point must not reach DWSIM.')
    support = next(s for s in source['supports'] if s['name'] == request['support'])
    membrane = load_case(jl, support['body_resistance_Pa_s_m3'],
        source['supports'][0]['body_resistance_Pa_s_m3'], source['protocol'],
        source['literature_reference']['fit']['fixed_impedance_Pa_s_m'])
    if abs(membrane['flow_m3_s'] / point['permeate_reference_flow_m3_s'] - 1) > 1e-12:
        raise ValueError('Changed process flux/area requires revalidating the inherited membrane anchor.')
    process = record(jl, jl.crossflow_recycle_load(membrane['pressure_drop_Pa'],
        point['permeate_reference_flow_m3_s'], protocol['nominal_recovery'],
        point['module_to_permeate_ratio'], point['channel_resistance_Pa_s_m3'],
        hydraulics['feed_pump_efficiency_fraction'], hydraulics['circulation_pump_efficiency_fraction']))
    solute = record(jl, jl.crossflow_solute_balance(protocol['module_inlet_concentration_kg_m3'],
        process['permeate_reference_flow_m3_s'], protocol['nominal_recovery'],
        point['module_to_permeate_ratio'], protocol['assumed_rejection']))
    return {'source_point': point, 'membrane': membrane, 'process': process,
            'passive_solute_bookkeeping': solute}


def selections(cases, decisions, protocol):
    selected = []
    for convention in protocol['flow_conventions']:
        for support in sorted({c['request']['support'] for c in cases if not c['request']['validation_only']}):
            for margin in protocol['capacity_fractions']:
                for policy in ('table_nominal', 'cross_source_monotonic'):
                    permitted = {d['row_id'] for d in decisions if d['capacity_fraction'] == margin and
                        (d['table']['nominal_admissible'] if policy == 'table_nominal' else d['cross_source_admissible'])}
                    candidates = [c for c in cases if not c['request']['validation_only'] and
                        c['request']['flow_convention'] == convention['id'] and c['request']['support'] == support and
                        c['request']['row_id'] in permitted]
                    best = min(candidates, key=lambda c: c['validation']['total_shaft_power_W'], default=None)
                    selected.append({'flow_convention': convention['id'], 'support': support,
                        'capacity_fraction': margin, 'policy': policy,
                        'case_id': None if best is None else best['request']['case_id'],
                        'velocity_m_s': None if best is None else best['translation']['source_point']['velocity_m_s'],
                        'specific_shaft_energy_kWh_m3': None if best is None else best['validation']['specific_shaft_energy_kWh_m3_reference_volume']})
    return selected
