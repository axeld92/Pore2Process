"""Nonnegative transit-time distribution and its finite-pulse outlet curve."""
import numpy as np


def paths(cumulative, dt, pulse_s, transit_s):
    """Fit from one dose-normalized cumulative curve using its known square input."""
    intervals = round(pulse_s / dt)
    if intervals < 1 or not np.isclose(intervals * dt, pulse_s):
        raise ValueError("The inlet pulse must span a whole number of time steps")
    mass = np.diff(cumulative)
    if np.any(mass < -1e-12) or cumulative[-1] < 0.999:
        raise ValueError("Calibration requires nonnegative, near-complete recovery")
    mass = np.maximum(mass, 0)
    mass /= mass.sum()
    inlet = np.ones(intervals) / intervals
    # A is convolution by inlet. A transpose is correlation. The fixed
    # first-difference penalty stabilizes modes averaged away by the pulse.
    penalty = 0.002
    weights = mass[:len(mass) - intervals + 1].copy()
    extrapolate = weights.copy()
    acceleration = 1.
    for _ in range(2000):
        residual = np.convolve(extrapolate, inlet) - mass
        gradient = np.correlate(residual, inlet, mode="valid")
        slope = penalty * np.diff(extrapolate)
        gradient[:-1] -= slope
        gradient[1:] += slope
        updated = np.maximum(extrapolate - gradient / (1 + 4 * penalty), 0)
        next_acceleration = (1 + np.sqrt(1 + 4 * acceleration**2)) / 2
        extrapolate = updated + (acceleration - 1) / next_acceleration * (updated - weights)
        weights, acceleration = updated, next_acceleration
    weights /= weights.sum()
    t = np.arange(len(weights)) * dt
    # Correct liquid inventory on the slow paths. Taper at the measured
    # window edge to keep every delay within the calibrated support.
    adjustment = np.clip(t / transit_s - 1, 0, 1) * np.clip((t[-1] - t) / transit_s, 0, 1)
    mean_error = transit_s - np.dot(t, weights)
    amplitude = mean_error / np.dot(adjustment, weights) if mean_error else 0.
    delays = t + amplitude * adjustment
    if np.any(np.diff(delays) <= 0):
        raise ValueError("Inventory correction changed the transit-time ordering")
    return delays, weights, float(amplitude)


def pulse_cdf(times, delays, weights, pulse_s):
    """Cumulative exit mass for a unit dose across parallel delay paths."""
    mass = np.r_[0., np.cumsum(weights)]
    first = np.r_[0., np.cumsum(weights * delays)]
    # Integrated impulse CDF, evaluated at both edges of the input box.
    def integral(t):
        j = np.searchsorted(delays, t, side="right")
        return t * mass[j] - first[j]
    return np.clip((integral(times) - integral(times - pulse_s)) / pulse_s, 0., 1.)


def moments(delays, weights, pulse_s):
    mean = np.dot(delays, weights)
    variance = np.dot((delays - mean)**2, weights)
    return float(mean + pulse_s / 2), float(variance + pulse_s**2 / 12)
