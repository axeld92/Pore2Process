module BeadFilm

import ..LorenzoPackedBed
import ..LorenzoPackedBed: rate, rate_terms

export FilmKinetics, film_coefficients, surface_progress, inlet_balance

"""Liquid-side film coefficient, using the superficial-velocity Reynolds number."""
function film_coefficients(u_m_s, diameter_m, rho_kg_m3, mu_Pa_s, diffusivity_m2_s,
                           void_fraction, correlation)
    properties = (u_m_s, diameter_m, rho_kg_m3, mu_Pa_s, diffusivity_m2_s)
    all(x -> isfinite(x) && x > 0, properties) && 0 < void_fraction < 1 ||
        throw(ArgumentError("Use positive liquid properties, speed, diameter and a physical void fraction."))
    reynolds = rho_kg_m3 * u_m_s * diameter_m / mu_Pa_s
    schmidt = mu_Pa_s / (rho_kg_m3 * diffusivity_m2_s)
    sherwood = if correlation == "wilson_geankoplis"
        1.09 / void_fraction * cbrt(reynolds * schmidt)
    elseif correlation == "corrected_printed_B2"
        1.17 * reynolds^0.585 * cbrt(schmidt)
    elseif correlation == "isolated_sphere"
        2.0
    else
        throw(ArgumentError("Unknown film correlation: $correlation"))
    end
    (; reynolds, schmidt, sherwood, coefficient_m_s=sherwood * diffusivity_m2_s / diameter_m)
end

"""Scalar surface progress for a reversible reaction with one common film mobility."""
function surface_progress(reaction, bulk, equilibrium, capacity)
    lower, upper = minmax(bulk, equilibrium)
    # At bulk progress the residual is minus the signed rate. At equilibrium
    # it is the film term alone. The bracket also handles reverse reaction.
    for _ in 1:44
        midpoint = (lower + upper) / 2
        if capacity * (midpoint - bulk) > reaction(midpoint)
            upper = midpoint
        else
            lower = midpoint
        end
    end
    (lower + upper) / 2
end

struct FilmKinetics
    kinetics::LorenzoPackedBed.PingPong
    capacity_mol_L_min::Float64
    eta_internal::Float64
    function FilmKinetics(kinetics, capacity_mol_L_min, eta_internal)
        isfinite(capacity_mol_L_min) && capacity_mol_L_min > 0 && 0 < eta_internal <= 1 ||
            throw(ArgumentError("Use positive film capacity and internal effectiveness in (0, 1]."))
        new(kinetics, capacity_mol_L_min, eta_internal)
    end
end

rate_terms(model::FilmKinetics, args...) = rate_terms(model.kinetics, args...)

function surface_state(model, terms, acid0, alcohol0, catalyst_g_L, bulk)
    ratio = alcohol0 / acid0
    equilibrium = 2ratio / (1 + ratio + sqrt((ratio - 1)^2 + 4ratio * terms.k2 / terms.k1))
    reaction(x) = model.eta_internal * rate(model.kinetics, terms, acid0, alcohol0, catalyst_g_L, x)
    surface = surface_progress(reaction, bulk, equilibrium, model.capacity_mol_L_min)
    (; surface, reaction=reaction(surface))
end

function rate(model::FilmKinetics, terms, acid0, alcohol0, catalyst_g_L, bulk)
    # The original axial solver multiplies the returned rate by eta_internal.
    # The same factor enters the local film balance exactly once.
    surface_state(model, terms, acid0, alcohol0, catalyst_g_L, bulk).reaction / model.eta_internal
end

function inlet_balance(model::FilmKinetics, temperature_C, catalyst_g_L; acid0=3.3, alcohol0=6.6)
    terms = rate_terms(model, temperature_C, acid0, alcohol0, false)
    state = surface_state(model, terms, acid0, alcohol0, catalyst_g_L, 0.0)
    uncoated = model.eta_internal * rate(model.kinetics, terms, acid0, alcohol0, catalyst_g_L, 0.0)
    (; surface_acid_relative_to_bulk=1 - state.surface,
       rate_relative_to_no_film=state.reaction / uncoated,
       film_balance_error_mol_L_min=abs(state.reaction - model.capacity_mol_L_min * state.surface))
end

end
