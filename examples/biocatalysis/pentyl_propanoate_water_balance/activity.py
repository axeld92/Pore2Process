"""Composition projections anchored to water activities and measured tie lines."""

from dataclasses import dataclass

import numpy as np

from examples.biocatalysis.pentyl_propanoate_transport.phase.run import ternary_split


def liquid_moles(conversion, water_increment_per_mol_acid=0.0):
    """Moles of acid, pentanol, water and ester per mole of acid charged."""
    water = conversion + water_increment_per_mol_acid
    if not np.isfinite(water) or not 0 <= conversion < 1 or water < 0:
        raise ValueError("Use conversion below one and finite, nonnegative liquid water.")
    return np.array([1 - conversion, 2 - conversion, water, conversion])


def measured_pentanol_curve(observations):
    """Hot VLE surrogate with the dry limit and the reported gas corrections."""
    _, liquid_pentanol, _, gamma_water, _ = np.array(observations).T
    water = np.r_[0.0, 1 - liquid_pentanol]
    activity = np.r_[0.0, water[1:] * gamma_water]
    if not (np.all(np.diff(water) > 0) and np.all(np.diff(activity) > 0)
            and 0 < water[-1] < activity[-1] < 1):
        raise ValueError("The tabulated VLE branch must be ordered and physically admissible.")
    return water, activity


def binary_activity(moles, curve):
    water_fraction = moles[2] / moles.sum()
    water, activity = curve
    if water_fraction > water[-1]:
        return None
    return float(np.interp(water_fraction, water, activity))


@dataclass(frozen=True)
class TernaryProjection:
    endpoints: np.ndarray

    @classmethod
    def from_tie_lines(cls, rows):
        pairs = np.array(rows, dtype=float).reshape(-1, 2, 3)
        pairs /= pairs.sum(axis=2, keepdims=True)
        organic = pairs[:, 0, :]
        acid_ratios = organic[:, 0] / (1 - organic[:, 2])
        if not (np.all(np.diff(acid_ratios) > 0) and np.all(pairs[:, 1, 2] > organic[:, 2])):
            raise ValueError("Ordered organic compositions and water-rich partners are required.")
        return cls(pairs)

    def state(self, moles, ester_as_pentanol, curve):
        acid, alcohol, water, ester = moles
        if ester_as_pentanol:
            alcohol += ester
        feed = np.array([acid, alcohol, water])
        # Intersect the organic binodal with a ray at fixed acid:alcohol.
        acid_ratio = acid / (acid + alcohol)
        organic = self.endpoints[:, 0, :]
        offset = organic[:, 0] - acid_ratio * (1 - organic[:, 2])
        if offset[0] > 0 or offset[-1] < 0:
            return None
        right = max(1, int(np.flatnonzero(offset >= 0)[0]))
        blend = -offset[right - 1] / (offset[right] - offset[right - 1])
        saturated_organic, saturated_aqueous = (
            self.endpoints[right - 1] + blend * (self.endpoints[right] - self.endpoints[right - 1])
        )
        relative_water = (water / feed.sum()) / saturated_organic[2]
        if relative_water <= 1 + 1e-12:
            relative_water = min(relative_water, 1.0)
            if curve is not None:
                binary_water, binary_aw = curve
                relative_water = np.interp(relative_water * binary_water[-1], binary_water, binary_aw) / binary_aw[-1]
            # Water-rich branch as an approximately ideal aqueous reference.
            return float(saturated_aqueous[2] * relative_water), False
        phases = ternary_split(feed, self.endpoints.reshape(-1, 6))
        if phases is None:
            return None
        _, aqueous = phases
        return float(aqueous[2] / aqueous.sum()), True


def projected_activity(moles, weighted_projections, ester_as_pentanol, curve):
    states = [(weight, projection.state(moles, ester_as_pentanol, curve))
              for weight, projection in weighted_projections]
    if any(state is None for _, state in states):
        return None
    return sum(weight * state[0] for weight, state in states)
