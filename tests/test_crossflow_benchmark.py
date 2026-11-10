"""Actual Julia source constraints, process coupling and rejection checks."""
import copy
import json

import pytest

from examples.filtration.literature_benchmark.model import prepare, requests, translate
from examples.filtration.literature_benchmark.check import verify
from pore2process_demo.bridge import runtime
from pore2process_demo.paths import ROOT; from pore2process_demo.io import read_json as read


def inputs():
    p = read(ROOT / 'examples/filtration/literature_benchmark/protocol.json')
    e = read(ROOT / p['extraction'])
    h = read(ROOT / 'examples/filtration/crossflow_recycle/protocol.json')
    h['recycle']['maximum_iterations'] = p['maximum_recycle_iterations']
    s = read(ROOT / 'docs/crossflow_recycle_reference.json')['source_for_Julia_replay']
    return p, e, h, s


def test_actual_julia_source_policies_preserve_transition_uncertainty():
    p, e, _, s = inputs(); jl = runtime()
    d = prepare(jl, e, p)
    assert len(d) == 12
    for margin in p['capacity_fractions']:
        group = [v for v in d if v['capacity_fraction'] == margin]
        assert not group[0]['table']['nominal_admissible']
        assert group[1]['table']['nominal_admissible'] is True
        assert group[1]['figure_monotonic_bracket']['classification'] == 'indeterminate'
        assert group[2]['cross_source_admissible'] is True
        assert group[1]['table']['source_specific_axial_energy_J_m3'] == pytest.approx(129.2134831460674)
    jobs = list(requests(d, p, s))
    assert len(jobs) == 32 and sum(not j['validation_only'] for j in jobs) == 30
    assert {j['flow_factor'] for j in jobs} == {1., 7.}
    assert all(j['row_id'] != 1 for j in jobs)
    assert len({j['case_id'] for j in jobs}) == len(jobs)


def test_actual_julia_flow_convention_and_finite_pass_solute_balance():
    p, e, h, s = inputs(); jl = runtime()
    request = {'row_id': 2, 'support': 'straight', 'flow_factor': 1.}
    a = translate(jl, e, p, h, s, request)
    b = translate(jl, e, p, h, s, {**request, 'flow_factor': 7.})
    assert a['source_point']['module_reference_flow_m3_s'] == pytest.approx(25e-3 / 3600 / .03)
    assert a['process']['channel_pressure_drop_Pa'] == pytest.approx(69.)
    assert a['source_point']['scenario_specific_axial_energy_J_m3'] == pytest.approx(575.)
    assert b['process']['circulation_shaft_power_W'] == pytest.approx(7 * a['process']['circulation_shaft_power_W'])
    assert b['process']['feed_shaft_power_W'] == pytest.approx(a['process']['feed_shaft_power_W'])
    assert b['membrane'] == a['membrane']
    assert a['passive_solute_bookkeeping']['module_inlet_concentration_kg_m3'] == .01
    assert a['passive_solute_bookkeeping']['permeate_concentration_kg_m3'] == 0
    # Overall fully retained-solute balance, independent of the module algebra.
    solute = a['passive_solute_bookkeeping']
    assert solute['fresh_concentration_kg_m3'] == pytest.approx(.5 * solute['retentate_and_purge_concentration_kg_m3'])
    assert b['passive_solute_bookkeeping']['fresh_concentration_kg_m3'] != solute['fresh_concentration_kg_m3']
    bad = copy.deepcopy(request); bad['row_id'] = 1
    with pytest.raises(ValueError, match='Inadmissible'):
        translate(jl, e, p, h, s, bad)
    changed = copy.deepcopy(p); changed['target_flux_LMH'] = 50.
    with pytest.raises(ValueError, match='Changed process flux'):
        translate(jl, e, changed, h, s, request)


def test_mutated_saved_source_fails_before_runtime_replay(tmp_path):
    from pore2process_demo.io import sha
    f = tmp_path / 'inputs.json'; f.write_text('{}\n')
    summary = {'status': 'passed', 'artifact_sha256': {'inputs.json': sha(f)}}
    (tmp_path / 'summary.json').write_text(json.dumps(summary))
    f.write_text('{"flow_factor": 1}\n')
    with pytest.raises(ValueError, match='Artifact hash mismatch'):
        verify(tmp_path / 'summary.json')


def test_actual_saved_high_circulation_dwsim_case_matches_julia():
    from examples.filtration.crossflow_recycle.model import validate
    r = read(ROOT / 'docs/crossflow_benchmark_reference.json')
    c = r['actual_high_circulation_case']
    t = translate(runtime(), r['extraction'], r['protocol'], r['hydraulics'],
                  r['source_for_Julia_replay'], c['request'])
    assert t['process']['module_to_permeate_ratio'] == pytest.approx(560 / 3)
    assert c['dwsim']['recycle']['maximum_iterations'] == 5000
    checks = validate(c['dwsim'], c['dwsim']['inputs'], t['process'],
                      r['hydraulics'], r['source_for_Julia_replay'])
    assert all(checks['checks'].values())
    assert checks['total_shaft_power_W'] == pytest.approx(t['process']['total_shaft_power_W'], rel=1e-4)


@pytest.mark.parametrize('change', ['false_convergence', 'lost_heat', 'lost_mass'])
def test_high_recycle_flags_cannot_hide_invalid_saved_states(change):
    from examples.filtration.crossflow_recycle.model import validate
    r = read(ROOT / 'docs/crossflow_benchmark_reference.json')
    c = r['actual_high_circulation_case']; actual = copy.deepcopy(c['dwsim'])
    if change == 'false_convergence':
        actual['recycle']['converged'] = False
    elif change == 'lost_heat':
        actual['powers']['cooling_removed_W'] = actual['powers']['cooling_energy_stream_W'] = 0
    else:
        actual['streams']['Return']['mass_flow_kg_s'] *= .99
    with pytest.raises((ValueError, RuntimeError)):
        validate(actual, actual['inputs'], c['translation']['process'],
                 r['hydraulics'], r['source_for_Julia_replay'])
