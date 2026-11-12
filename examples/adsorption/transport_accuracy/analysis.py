"""Read native results; distinguish local drift, liquid inventory and RTD error."""
import json
from pathlib import Path

import numpy as np

from examples.adsorption.column_transport.run import analyze as strict_analysis
from examples.adsorption.packed_bed_tracer.bridge import moments


def spatial_diagnostic(path):
    a = np.genfromtxt(Path(path) / 'final.csv', delimiter=',', names=True)
    raw = json.loads((Path(path) / 'raw.json').read_text())
    shape = (raw['nx'], raw['n'], raw['n'])
    c = a['concentration_kg_m3'].reshape(shape)
    fluid = a['material'].reshape(shape) != 2
    near = np.zeros(shape, dtype=bool)
    for axis in range(3):
        for shift in (-1, 1):
            near |= ~np.roll(fluid, shift, axis=axis)
    near &= fluid
    values = c[fluid]
    return dict(physical_fluid_mean_kg_m3=float(np.mean(values)),
                physical_fluid_RMS_departure=float(np.sqrt(np.mean((values - 1.) ** 2))),
                fraction_fluid_nodes_next_to_solid=float(np.mean(near[fluid])),
                mean_absolute_departure_next_to_solid=float(np.mean(abs(c[near] - 1.))) if np.any(near) else None,
                mean_absolute_departure_away_from_solid=float(np.mean(abs(c[fluid & ~near] - 1.))) if np.any(fluid & ~near) else None,
                proximity_definition='Six immediate Cartesian neighbours; periodic roll; diagnostic of voxel walls, not surface distance.')


def analyze(path):
    path = Path(path)
    r = strict_analysis(path, references=False)
    raw = r['raw']
    a = np.genfromtxt(path / 'inventory.csv', delimiter=',', names=True)
    f = np.genfromtxt(path / 'flux.csv', delimiter=',', names=True)
    constant = raw['initial_concentration'] == 1.
    quarter = a['time_s'] >= .75 * a['time_s'][-1]
    fluid0 = float(a['fluid_tracer_kg'][0])
    Q = raw['flow_m3_s']
    V = raw['fluid_volume_m3']
    mass = a['fluid_tracer_kg'] + a['wall_population_tracer_kg']
    checks = dict(finite=r['checks']['finite'], kinetic_mass=r['checks']['mass'],
                  concentration_nonnegative=r['minimum_concentration_kg_m3'] >= -1e-8,
                  concentration_upper_bound=r['maximum_concentration_kg_m3'] <= 1.02)
    liquid_change = None
    drift = None
    diagnostic = None
    m = None
    observed_moments = None
    inlet_moments = None
    moments_error = None
    plateau_error = None
    if constant:
        liquid_change = float(np.max(abs(a['fluid_tracer_kg'] / fluid0 - 1.)))
        drift = max(float(np.ptp(a['min_concentration_kg_m3'][quarter])),
                    float(np.ptp(a['max_concentration_kg_m3'][quarter])))
        checks.update(local_departure=r['uniform_maximum_departure'] <= .02,
                      liquid_inventory=liquid_change <= .001, bounded_late_drift=drift <= .002)
        diagnostic = spatial_diagnostic(path)
    else:
        try:
            m = moments(f['time_s'].tolist(), f['outlet_flux_kg_s'].tolist())
            observed_moments = dict(m)
            if 'inlet_ramp_s' in raw:
                feed = np.genfromtxt(path / 'feed.csv', delimiter=',', names=True)
                np.testing.assert_allclose(feed['time_s'], f['time_s'], rtol=1e-12, atol=1e-14)
                inlet_moments = moments(feed['time_s'].tolist(), feed['prescribed_concentration_kg_m3'].tolist())
                T = raw['inlet_ramp_s']
                checks['inlet_moment_audit'] = (
                    abs(inlet_moments['mean_s'] / (T / 2) - 1) <= 1e-5 and
                    abs(inlet_moments['variance_s2'] / (T*T/28) - 1) <= 1e-4)
                m = remove_input_moments(observed_moments, inlet_moments)
        except Exception as error:
            # Preserve the solver output and the reason; do not replace the tail.
            moments_error = str(error)
            m = None
        if m is None:
            checks.update(positive_RTD_moments=False, absolute_plateau=False,
                          completed_front=False, first_moment_tail=False, second_moment_tail=False)
        else:
            plateau_error = abs(m['steady_exit_flux_kg_s'] / Q - 1.)
            checks.update(positive_RTD_moments=True, absolute_plateau=plateau_error <= .01,
                          completed_front=m['final_fraction'] >= .999,
                          first_moment_tail=m['tail_area_fraction'] <= .01,
                          second_moment_tail=m['tail_second_moment_fraction'] <= .01,
                          outlet_flux_bounds=m['minimum_fraction'] >= -1e-8 and m['maximum_fraction'] <= 1.02)
    final_wall_change = float(a['wall_population_tracer_kg'][-1] - a['wall_population_tracer_kg'][0])
    extrema_locations = None
    if 'min_i' in a.dtype.names:
        extrema_locations = {}
        for role, quantity, select in [('minimum', 'min_concentration_kg_m3', np.argmin),
                                       ('maximum', 'max_concentration_kg_m3', np.argmax)]:
            row = a[int(select(a[quantity]))]
            prefix = 'min' if role == 'minimum' else 'max'
            extrema_locations[role] = dict(time_s=float(row['time_s']),
                concentration_kg_m3=float(row[quantity]),
                zero_based_indices=[int(row[prefix + '_' + k]) for k in ('i', 'j', 'k')],
                role='Location of recorded liquid extremum; does not identify its cause.')
    return dict(strict_diagnostic=r, engineering_checks=checks,
                engineering_screen_passed=all(checks.values()),
                actual_duration_s=float(a['time_s'][-1]),
                physical_void_time_s=V / Q,
                observed_transits=float(a['time_s'][-1] / (V / Q)),
                maximum_average_liquid_concentration_relative_change=liquid_change,
                final_quarter_extrema_drift=drift, final_spatial_diagnostic=diagnostic,
                final_wall_population_inventory_change_kg=final_wall_change,
                final_wall_population_change_over_physical_liquid_plateau_inventory=final_wall_change / V,
                final_physical_liquid_inventory_kg=float(a['fluid_tracer_kg'][-1]),
                extrema_locations=extrema_locations,
                maximum_global_kinetic_inventory_kg=float(np.max(mass)),
                RTD=m, observed_outlet_moments=observed_moments, imposed_inlet_moments=inlet_moments,
                RTD_unavailable_reason=moments_error,
                absolute_outlet_plateau_relative_error=plateau_error,
                note='Kinetic inventory includes temporary wall populations; physical liquid volume excludes all solid nodes. A finite-window RTD is not a bound on unobserved tails.')


def remove_input_moments(response, inlet):
    """Independent convolution identity, without fitting or tail extrapolation."""
    mean = response['mean_s'] - inlet['mean_s']
    variance = response['variance_s2'] - inlet['variance_s2']
    if mean <= 0 or variance <= 0:
        raise ValueError('Input-adjusted finite-window transport moments are nonpositive.')
    result = dict(response)
    result.update(mean_s=mean, variance_s2=variance, equivalent_sections=mean*mean/variance,
        tail_area_fraction=response['tail_area_fraction'] * response['mean_s'] / mean,
        tail_second_moment_fraction=response['tail_second_moment_fraction'] *
            (response['variance_s2'] + response['mean_s']**2) / (variance + mean**2))
    return result


def refinement_metrics(a, b):
    if a['RTD'] is None or b['RTD'] is None:
        return dict(available=False, passed=False)
    delta = {k: abs(b['RTD'][k] / a['RTD'][k] - 1.) for k in ('mean_s', 'variance_s2')}
    return dict(available=True, relative_changes=delta,
                passed=delta['mean_s'] <= .05 and delta['variance_s2'] <= .1)


def service_gate(cases, required_case_ids):
    """Engineering tolerance never silently becomes the historical strict gate."""
    missing = sorted(set(required_case_ids) - set(cases))
    failed = {k: [name for name, passed in c['engineering_checks'].items() if not passed]
              for k, c in cases.items() if not c['engineering_screen_passed']}
    return dict(eligible_for_closure_assessment=not missing and not failed,
                missing_cases=missing, failed_native_screens=failed,
                closure_qualified=False, service_prediction_qualified=False,
                service_time_s=None, service_volume_m3=None, energy_per_on_spec_volume_Wh_m3=None,
                reason='A native screen is a prerequisite, not a qualified closure or a service prediction.')
