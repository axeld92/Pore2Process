module ParticleTransport

export effectiveness, modulus_for, relative_effectiveness

"""First-order efficiency of a spherical bead with the same total intrinsic rate.

`modulus` uses the full bead volume for the intrinsic rate density.
`core` is the inactive inner radius divided by the bead radius. The rate
density in the active shell is increased to preserve the catalyst inventory.
"""
function effectiveness(modulus, core=0.0)
    isfinite(modulus) && modulus >= 0 && 0 <= core < 1 ||
        throw(ArgumentError("Use a finite nonnegative modulus and an inactive core in [0, 1)."))
    volume_factor = 1 + core + core^2
    if modulus < 1e-2
        # Small-modulus expansion of the volume-averaged shell concentration.
        curvature = (1-core) * (1 + 3core + 6core^2 + 5core^3) / (15volume_factor^2)
        return 1 - curvature * modulus^2
    end
    shell_modulus = modulus / sqrt((1-core) * volume_factor)
    t = tanh(shell_modulus * (1-core))
    surface_gradient = (1 + core*shell_modulus*t) / (core + t/shell_modulus) - 1
    3 * surface_gradient / modulus^2
end

"""Invert the uniform-sphere expression using the tabulated reference efficiency."""
function modulus_for(efficiency)
    isfinite(efficiency) && 0 < efficiency <= 1 ||
        throw(ArgumentError("Reference effectiveness must be in (0, 1]."))
    efficiency == 1 && return 0.0
    lower, upper = 0.0, 3 / efficiency
    for _ in 1:60
        mid = (lower + upper) / 2
        if effectiveness(mid) > efficiency
            lower = mid
        else
            upper = mid
        end
    end
    (lower + upper) / 2
end

"""Return a relative correction to an apparent rate from the reference bead state."""
function relative_effectiveness(modulus, diffusivity_ratio, diameter_ratio; core=0.0)
    isfinite(diffusivity_ratio) && diffusivity_ratio > 0 &&
        isfinite(diameter_ratio) && diameter_ratio > 0 ||
        throw(ArgumentError("Ratios of effective diffusivity and bead diameter must be positive."))
    # At fixed enzyme per particle the envelope's intrinsic rate density is
    # proportional to diameter^-3. Consequently modulus² scales as 1/(D * d).
    changed = modulus / sqrt(diffusivity_ratio * diameter_ratio)
    effectiveness(changed, core) / effectiveness(modulus, core)
end

end
