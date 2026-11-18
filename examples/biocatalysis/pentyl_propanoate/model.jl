module PentylPropanoate

using LinearAlgebra: I, Tridiagonal, norm

export Kinetics, coefficients, fit_batch, fit_at_temperature, batch_conversion, bed_prediction,
       equilibrium_conversion, rate_constant, equilibrium_constant, axial_peclet

const ACID_MOLAR_MASS_KG_MOL = 0.0740785
const REFERENCE_T_K = 333.15
const TEMPERATURE_SCALE_K = REFERENCE_T_K^2 / 20

"""Apparent reversible kinetics per kilogram of immobilised catalyst, for 2:1 feed."""
struct Kinetics
    log_k_reference::Float64
    log_k_slope::Float64
    log_K_reference::Float64
    log_K_slope::Float64
end

coefficients(kinetics::Kinetics) = [
    kinetics.log_k_reference, kinetics.log_k_slope,
    kinetics.log_K_reference, kinetics.log_K_slope,
]

temperature_coordinate(temperature_C) =
    TEMPERATURE_SCALE_K * (inv(REFERENCE_T_K) - inv(temperature_C + 273.15))

rate_constant(kinetics, temperature_C) =
    exp(kinetics.log_k_reference + kinetics.log_k_slope * temperature_coordinate(temperature_C))

equilibrium_constant(kinetics, temperature_C) =
    exp(kinetics.log_K_reference + kinetics.log_K_slope * temperature_coordinate(temperature_C))

equilibrium_conversion(kinetics, temperature_C) =
    4 / (3 + sqrt(1 + 8 / equilibrium_constant(kinetics, temperature_C)))

function conversion(kinetics, temperature_C, exposure)
    rate = rate_constant(kinetics, temperature_C) / 9
    gap = sqrt(1 + 8 / equilibrium_constant(kinetics, temperature_C))
    change = -expm1(-gap * rate * exposure)
    2 * change / (gap + (3 - gap) * change / 2)
end

"""Acid conversion from a dry, product-free 2:1 feed after `time_min` in batch.

The catalyst loading is mass percent relative to the initial acid mass.
Stoichiometry conserves total moles and creates one mole each of water and ester.
"""
function batch_conversion(kinetics, temperature_C, loading_pct, time_min)
    exposure = loading_pct / 100 * ACID_MOLAR_MASS_KG_MOL * (60 * time_min)
    conversion(kinetics, temperature_C, exposure)
end

function batch_exposures(loading_pct, time_min, observed)
    length(loading_pct) == length(time_min) == length(observed) ||
        throw(ArgumentError("Batch columns must have equal lengths."))
    length(observed) >= 4 && all(x -> isfinite(x) && 0 <= x <= 1, observed) ||
        throw(ArgumentError("Supply at least four finite batch conversions in [0, 1]."))
    all(w -> isfinite(w) && w > 0, loading_pct) &&
        all(t -> isfinite(t) && t >= 0, time_min) ||
        throw(ArgumentError("Catalyst loading and time must be physical."))
    loading_pct .* time_min .* (0.6 * ACID_MOLAR_MASS_KG_MOL)
end

function least_squares(predict, observed, initial)
    p = copy(initial)
    residual = predict(p) .- observed
    damping = 1e-3
    for iteration in 1:100
        jacobian = zeros(length(observed), length(p))
        for j in eachindex(p)
            shift = zeros(length(p))
            shift[j] = cbrt(eps(Float64)) * max(1.0, abs(p[j]))
            jacobian[:, j] = (predict(p + shift) - predict(p - shift)) / (2 * shift[j])
        end
        gradient = jacobian' * residual
        if norm(gradient, Inf) <= 1e-9
            return (; parameters=p, predicted=residual .+ observed,
                    iterations=iteration - 1, gradient_inf=norm(gradient, Inf))
        end
        step = -(jacobian' * jacobian + damping * I) \ gradient
        trial = predict(p + step) .- observed
        if norm(trial) < norm(residual)
            p += step
            residual = trial
            damping *= 0.3
        else
            damping *= 10
        end
    end
    error("The batch-only least-squares fit did not converge.")
end

"""Fit one catalyst-mass-specific rate with an apparent Kx from the matching batch."""
function fit_at_temperature(temperature_C, loading_pct, time_min, observed, apparent_K)
    isfinite(temperature_C) && temperature_C > -273.15 &&
        isfinite(apparent_K) && apparent_K > 0 ||
        throw(ArgumentError("Use a physical temperature and a positive apparent Kx."))
    exposures = batch_exposures(loading_pct, time_min, observed)
    kinetics(p) = Kinetics(only(p), 0, log(apparent_K), 0)
    predict(p) = conversion.(Ref(kinetics(p)), temperature_C, exposures)
    fit = least_squares(predict, observed, [log(0.12)])
    (; kinetics=kinetics(fit.parameters), fit.predicted, fit.iterations, fit.gradient_inf)
end

"""Fit a smooth temperature dependence of the rate and apparent Kx to batch only."""
function fit_batch(temperature_C, loading_pct, time_min, observed)
    length(temperature_C) == length(observed) &&
        all(t -> isfinite(t) && t > -273.15, temperature_C) ||
        throw(ArgumentError("Supply a physical temperature for every batch observation."))
    exposures = batch_exposures(loading_pct, time_min, observed)
    predict(p) = conversion.(Ref(Kinetics(p...)), temperature_C, exposures)
    fit = least_squares(predict, observed, [log(0.12), 0.5, log(5.5), 0.3])
    (; kinetics=Kinetics(fit.parameters...), fit.predicted, fit.iterations, fit.gradient_inf)
end

# Danckwerts inlet and zero-gradient outlet flux conditions. This is the
# dimensionless pulse-free first-passage variance for the axial dispersion model.
function relative_variance(peclet)
    peclet < 1e-3 && return 1 - peclet / 3 + peclet^2 / 12 - peclet^3 / 60
    z = inv(peclet)
    2z * (1 + z * expm1(-peclet))
end

"""Moment-match a dimensionless tracer variance to the Danckwerts axial Péclet number."""
function axial_peclet(variance)
    isfinite(variance) && 0 < variance < 1 ||
        throw(ArgumentError("The pulse-free relative variance must be between zero and one."))
    lo, hi = 0.0, 2 / variance
    for _ in 1:60
        mid = (lo + hi) / 2
        if relative_variance(mid) > variance
            lo = mid
        else
            hi = mid
        end
    end
    (lo + hi) / 2
end

"""Predict the outlet acid conversion at known inlet acid molar flow and catalyst mass.

With finite `peclet`, solve the conservative axial dispersion balance over a
uniform catalyst inventory. `peclet=Inf` gives the exact plug-flow solution.
"""
function bed_prediction(kinetics, temperature_C, catalyst_mass_kg, acid_feed_mol_s;
                        peclet=Inf, cells=256)
    isfinite(temperature_C) && temperature_C > -273.15 &&
        isfinite(catalyst_mass_kg) && catalyst_mass_kg >= 0 &&
        isfinite(acid_feed_mol_s) && acid_feed_mol_s > 0 && peclet > 0 && cells >= 2 ||
        throw(ArgumentError("Use a physical temperature, inventory, acid flow, Péclet number and cell count."))
    exposure = catalyst_mass_kg / acid_feed_mol_s
    if isinf(peclet)
        return (; conversion=conversion(kinetics, temperature_C, exposure),
                balance_error=0.0)
    end
    equilibrium = equilibrium_conversion(kinetics, temperature_C)
    reaction = exposure * rate_constant(kinetics, temperature_C) / (9 * cells)
    quadratic = 1 - 1 / equilibrium_constant(kinetics, temperature_C)
    # Exponential flux between cell centres is monotone for every cell Péclet number.
    downstream = inv(expm1(peclet / cells))
    upstream = 1 + downstream
    transport = fill(upstream + downstream, cells)
    transport[1] = upstream
    transport[end] = 1 + downstream
    lower, upper = fill(-upstream, cells - 1), fill(-downstream, cells - 1)
    source(x) = reaction * (2 - 3x + quadratic * x^2)
    x = [conversion(kinetics, temperature_C, exposure * (i - 0.5) / cells) for i in 1:cells]
    tolerance = max(2e-11, 8eps(Float64) * cells * maximum(transport))
    for _ in 1:30
        internal_flux = @views x[1:end-1] .+ downstream .* (x[1:end-1] .- x[2:end])
        released = source.(x)
        residual = diff(vcat(0.0, internal_flux, x[end])) - released
        balance_error = abs(x[end] - sum(released))
        if norm(residual, 1) <= tolerance && balance_error <= 2e-11
            return (; conversion=x[end], balance_error)
        end
        slope = reaction .* (-3 .+ 2quadratic .* x)
        delta = Tridiagonal(lower, transport - slope, upper) \ residual
        x = clamp.(x - delta, 0.0, equilibrium)
    end
    error("The axial-dispersion balance did not converge.")
end

end
