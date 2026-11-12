"""Compare saved stationary-disk fields and reflection inventory with the archived controls."""
import json

import numpy as np

from examples.adsorption.stationary_curved_wall import reference as ref
from examples.adsorption.stationary_curved_wall.analysis import geometry, rms
from examples.adsorption.transport_verification.analysis import EXCURSION_LIMIT, RMS_LIMIT, SHELL_M, load

SETUPS = {
    'direct_fullway': (False, False),
    'regularized_halfway': (True, True),
    'direct_halfway': (False, True),
}


def analyze(directory, archived_directory, radial):
    raw, mask, velocity, fields, history, coordinates = load(directory)
    archived, old_mask, old_velocity, old_fields, _, _ = load(archived_directory)
    spec = json.loads((directory / 'request.json').read_text())['spec']
    regularized, halfway = SETUPS[spec['setup']]
    amp, steps, dt, dx = raw['amplitude_kg_m3'], raw['steps'], raw['dt_s'], raw['dx_m']
    if (raw['setup'], raw['n'], amp) != (spec['setup'], spec['n'], spec['amplitude']):
        raise ValueError('Native case differs from its request')
    matched = ('kind', 'descriptor', 'n', 'nz', 'steps', 'dt_s', 'dx_m', 'side_m', 'radius_m',
               'speed_m_s', 'diffusivity_m2_s', 'amplitude_kg_m3', 'magic', 'tau_even', 'tau_odd',
               'active_wall_links')
    for key in matched:
        if raw[key] != archived[key]:
            raise ValueError('Native case differs from the archived control in ' + key)
    expected_collision = ('stock CombinedAdvectionDiffusionRLBdynamics + ' if regularized else 'stock ')
    expected_wall = 'stock BouzidiPostProcessor q=0.5' if halfway else 'stock BounceBack'
    if raw['collision'] != expected_collision + 'AdvectionDiffusionTRTdynamics' or raw['wall'] != expected_wall:
        raise ValueError('Unexpected stock collision or reflection')
    for array in (velocity, *fields.values(), *(history[k] for k in history.dtype.names)):
        if not np.isfinite(array).all():
            raise ValueError('Non-finite native data')
    np.testing.assert_array_equal(mask, old_mask)
    np.testing.assert_allclose(velocity, old_velocity, rtol=0, atol=1e-14)
    fluid = mask == 1
    np.testing.assert_allclose(fields['initial'][fluid], old_fields['initial'][fluid], rtol=0, atol=1e-12)
    np.testing.assert_array_equal(history['step'], np.arange(steps+1))
    np.testing.assert_allclose(history['time_s'], np.arange(steps+1)*dt, rtol=1e-14, atol=0)
    # Half-way solid cells hold stream copies, not a physical concentration.
    z_error = max(float(np.max(abs(f-f[:, :, :1])[fluid])) for f in fields.values())
    if z_error > 1e-12:
        raise ValueError('Disk field is not axially invariant')
    for name, index in [('initial', 0), ('middle', steps//2), ('peak', raw['peak_step']), ('final', steps)]:
        field = fields[name][fluid]
        np.testing.assert_allclose([field.min(), field.max()],
            [history['min_concentration_kg_m3'][index], history['max_concentration_kg_m3'][index]],
            rtol=0, atol=2e-14)
    low = float(history['min_concentration_kg_m3'].min())
    high = float(history['max_concentration_kg_m3'].max())
    excursion = max(.5-amp-low, high-(.5+amp), 0.)
    peak_index = raw['peak_step']
    if abs(excursion-raw['maximum_excursion_kg_m3']) > 1e-14:
        raise ValueError('Peak excursion does not match history')
    adjacent, links, initial_wall_mass, link_velocity = geometry(directory, fluid, velocity, fields['initial'], dx)
    if links != raw['active_wall_links']:
        raise ValueError('Wall link count differs')
    liquid = history['fluid_inventory_kg']
    wall = history['active_wall_population_kg']
    closed = history['kinetic_inventory_kg']
    np.testing.assert_allclose(closed, liquid+wall, rtol=1e-14, atol=0)
    if halfway:
        coverage = json.loads((directory / 'halfway_links.json').read_text())
        if coverage['q'] != .5 or coverage['fluid_solid_links'] != links or any(
                coverage[k] for k in ('missing_links', 'unexpected_active_links', 'inactive_changed')):
            raise ValueError('Half-way link coverage failed')
        if np.any(wall):
            raise ValueError('Half-way inventory counts solid copies')
    elif abs(initial_wall_mass/wall[0]-1) > 1e-11:
        raise ValueError('Full-way initial storage disagrees with reflected link inventory')
    if abs(liquid[0]/(.5*np.count_nonzero(fluid)*dx**3)-1) > 1e-11:
        raise ValueError('Initial liquid inventory disagrees with the disk geometry')
    drift = float(np.max(abs(closed/closed[0]-1)))
    if drift > 1e-8:
        raise ValueError(f'Closed inventory drift {drift:.3g} exceeds 1e-8')
    x = coordinates[0][:, None, None]
    y = coordinates[1][None, :, None]
    shell = fluid & (np.hypot(x-ref.L/2, y-ref.L/2) >= ref.R-SHELL_M)
    checks = dict(engineering_bounds=excursion <= EXCURSION_LIMIT, physical_bounds=excursion <= 1e-8)
    result = dict(raw=raw, range_kg_m3=[low, high], maximum_excursion_kg_m3=excursion,
        peak_time_s=peak_index*dt, closed_inventory_relative_drift=drift,
        liquid_inventory_relative_drift=float(np.max(abs(liquid/liquid[0]-1))),
        axial_invariance_maximum_error=z_error, checks=checks, **link_velocity)
    if not amp:
        checks['constant'] = excursion <= 1e-9
        return result
    radial_coordinates, modes = radial
    result['fields'] = {}
    for key, time, mode in zip(('middle', 'final'), ref.TIMES, modes, strict=True):
        expected = ref.profile(x, y, time, amp, (radial_coordinates, mode))
        error = (fields[key]-expected)/amp
        stats = dict(time_s=time, RMS_over_amplitude=rms(error[fluid]),
            shell_RMS_over_amplitude=rms(error[shell]), interior_RMS_over_amplitude=rms(error[fluid & ~shell]),
            link_RMS_over_amplitude=rms(error[adjacent]),
            liquid_sd_kg_m3=float(np.std(fields[key][fluid])),
            reference_sd_at_cells_kg_m3=float(np.std(np.broadcast_to(expected, mask.shape)[fluid])))
        result['fields'][key] = stats
        checks[key+'_profile'] = stats['RMS_over_amplitude'] <= RMS_LIMIT
        checks[key+'_shell'] = stats['shell_RMS_over_amplitude'] <= RMS_LIMIT
    return result


def passes(constant, varying):
    # Physical bounds are tracked separately from the engineering acceptance.
    return all(passed for case in (constant, varying) for key, passed in case['checks'].items()
               if key != 'physical_bounds')
