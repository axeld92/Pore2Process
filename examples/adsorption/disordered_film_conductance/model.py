"""Turn first-contact pulses into independent liquid-film estimates."""
import numpy as np

from examples.adsorption.disordered_early_bypass.model import archive_lags


def spherical_area_by_half(centres_x, diameter, length):
    """Area of the spherical zones on each side of the periodic mid-plane."""
    if length <= diameter or diameter <= 0:
        raise ValueError("Invalid spherical packing dimensions")
    radius = diameter / 2
    images = np.mod(centres_x, length)[:, None] + length * np.array([-1, 0, 1])
    return np.array([
        2 * np.pi * radius * np.maximum(
            0., np.minimum(upper, images + radius) - np.maximum(lower, images - radius)
        ).sum()
        for lower, upper in ((0., length / 2), (length / 2, length))
    ])


def observations(raw):
    arrivals, hits = np.split(raw, 2, axis=1)
    starts = arrivals - archive_lags(arrivals)
    if np.isinf(hits).any() or np.any(hits < starts - 1e-12) or np.any(hits > arrivals + 1e-12):
        raise ValueError("Invalid solid-wall contact times")
    first_hit = np.min(np.where(np.isnan(hits), np.inf, hits), axis=1)
    half, end = arrivals[:, 3], arrivals[:, 7]
    through_half = np.minimum(first_hit, half)
    through_end = np.minimum(first_hit, end)
    captured_first = np.isfinite(hits[:, :4]).any(axis=1)
    captured_whole = np.isfinite(first_hit)
    counts = np.column_stack((captured_first, captured_whole & ~captured_first, captured_whole))
    times = np.column_stack((through_half, through_end - through_half, through_end))
    if np.any(times < 0):
        raise ValueError("Negative liquid occupation")
    return counts.astype(float), times, end


def film_rate(captures, durations, void_volume):
    """Uptake per volume-averaged liquid concentration, in m³/s."""
    if np.any(np.asarray(durations) <= 0) or np.any(np.asarray(captures) < 0) or np.any(void_volume <= 0):
        raise ValueError("Invalid capture, occupation time or liquid volume")
    return void_volume * captures / durations


def stats(raw, void_volume, surface_area, seed):
    captures, times, inert = observations(raw)
    rate = film_rate(captures.sum(axis=0), times.sum(axis=0), void_volume)
    # Whole-marker resampling preserves the relationship between capture and
    # remaining time in all three domains. It quantifies sampling error only.
    rng = np.random.default_rng(seed)
    samples = np.column_stack((captures, times))
    boot = np.empty((1200, 3))
    for draw in boot:
        total = samples[rng.integers(len(raw), size=len(raw))].sum(axis=0)
        draw[:] = film_rate(total[:3], total[3:], void_volume) / surface_area
    return dict(markers=len(raw), captured=captures.sum(axis=0).astype(int).tolist(),
                mean_capped_times_s=times.mean(axis=0).tolist(), inert_transit_s=float(inert.mean()),
                archive_film_conductance_m3_s=rate.tolist(),
                archive_film_coefficient_m_s=(rate / surface_area).tolist(),
                sampling_95pct_m_s=np.quantile(boot, [.025, .975], axis=0).tolist())


def effective_coefficient(k, direct_fraction, extra_resistance):
    """Parallel direct and slow areas with a fixed extra resistance in series."""
    return k * (direct_fraction + (1 - direct_fraction) / (1 + k * extra_resistance))


def required_coefficient(target, area, direct_fraction, extra_resistance):
    """Invert the increasing clean-film conductance curve by bisection."""
    if target <= 0 or area <= 0 or not 0 < direct_fraction <= 1 or extra_resistance < 0:
        raise ValueError("Invalid film threshold or access resistance")
    low, high = 0., target / (area * direct_fraction)
    for _ in range(55):
        mid = (low + high) / 2
        if area * effective_coefficient(mid, direct_fraction, extra_resistance) < target:
            low = mid
        else:
            high = mid
    return (low + high) / 2
