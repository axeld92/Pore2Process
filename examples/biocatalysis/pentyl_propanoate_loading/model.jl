module CatalystLoading

import ..PentylPropanoate

export rate_scale, progress_conversion, batch_conversion, fit_loading, predict_batch

const RATE_BOUNDS = (1e-4, 10.0)
const EQUILIBRIUM_BOUNDS = (1e-3, 1e3)

"""Mass-specific rate relative to the positive reference concentration."""
function rate_scale(law, shape, relative_concentration)
    u = relative_concentration
    isfinite(u) && u > 0 && isfinite(shape) ||
        throw(ArgumentError("Use a positive finite relative concentration and shape."))
    if law == "proportional"
        return 1.0
    elseif law == "power"
        0 <= shape <= 3 || throw(ArgumentError("Power exponent must be in [0, 3]."))
        return u^(shape-1)
    elseif law == "crowding"
        0 <= shape <= 1 || throw(ArgumentError("Crowding strength must be in [0, 1]."))
        return 1/(1 + shape*(u-1))
    end
    throw(ArgumentError("Unknown loading law: $law"))
end

"""
Exact dry-feed solution for X from a reversible mole-fraction rate law.
Exposure is kg catalyst*s/mol initial acid, k is mol/(kg*s), and
R is the initial alcohol/acid mole ratio. At R=2 this is the
preceding PentylPropanoate model's solution.
"""
function progress_conversion(k, Kx, R, exposure)
    all(isfinite, (k, Kx, R, exposure)) && k > 0 && Kx > 0 && R > 0 && exposure >= 0 ||
        throw(ArgumentError("Use positive rate, Kx and feed ratio, and nonnegative exposure."))
    gap = sqrt((R-1)^2 + 4R/Kx)
    change = -expm1(-gap*k*exposure/(1+R)^2)
    2R*change / (2gap + (1+R-gap)*change)
end

function batch_conversion(law, shape, k_ref, Kx, loading_pct, time_min, R;
                          reference_loading_pct=10.0)
    isfinite(loading_pct) && loading_pct > 0 && isfinite(time_min) && time_min >= 0 ||
        throw(ArgumentError("Use a positive loading and nonnegative elapsed time."))
    exposure = loading_pct/100 * PentylPropanoate.ACID_MOLAR_MASS_KG_MOL * 60time_min
    k = k_ref * rate_scale(law, shape, loading_pct/reference_loading_pct)
    progress_conversion(k, Kx, R, exposure)
end

function predict_batch(law, shape, k_ref, Kx, loadings, times, ratios;
                       reference_loading_pct=10.0)
    length(loadings) == length(times) == length(ratios) ||
        throw(ArgumentError("Batch input columns must have equal lengths."))
    [batch_conversion(law, shape, k_ref, Kx, w, t, R; reference_loading_pct)
     for (w,t,R) in zip(loadings, times, ratios)]
end

"""A deterministic bounded scalar profile; endpoints are candidate minima."""
function profile_minimum(evaluate, lower, upper, points)
    shapes = collect(range(lower, upper; length=points))
    profile = evaluate.(shapes)
    best = profile[argmin([f.weighted_sse for f in profile])]
    golden = (sqrt(5)-1)/2
    for i in eachindex(shapes)
        (i == 1 || profile[i].weighted_sse <= profile[i-1].weighted_sse) &&
            (i == length(shapes) || profile[i].weighted_sse <= profile[i+1].weighted_sse) || continue
        lo, hi = shapes[max(1,i-1)], shapes[min(length(shapes),i+1)]
        a = hi-golden*(hi-lo)
        b = lo+golden*(hi-lo)
        fa, fb = evaluate(a), evaluate(b)
        for _ in 1:60
            hi-lo < 1e-7 && break
            if fa.weighted_sse < fb.weighted_sse
                hi, b, fb = b, a, fa
                a = hi-golden*(hi-lo)
                fa = evaluate(a)
            else
                lo, a, fa = a, b, fb
                b = lo+golden*(hi-lo)
                fb = evaluate(b)
            end
        end
        for result in (fa, fb, evaluate((lo+hi)/2))
            result.weighted_sse < best.weighted_sse && (best = result)
        end
    end
    best
end

"""Guard ill-identified extreme-shape trials with a bounded numerical profile.

The wide numerical bounds include the limiting high-K region and prevent
the nuisance fit from diverging at shapes with little support. Selected
fits must be checked independently for boundary activity and refinement.
"""
function bounded_conditional_fit(law, shape, scaled, ratios, observed, sqrtweights, fixed_K)
    function at_equilibrium(log_K)
        K = exp(log_K)
        function at_rate(log_k)
            k = exp(log_k)
            cost = sum(w^2*(progress_conversion(k,K,R,s)-x)^2
                       for (w,R,s,x) in zip(sqrtweights,ratios,scaled,observed))
            (; law, shape, k_ref=k, Kx=K, weighted_sse=cost)
        end
        profile_minimum(at_rate, log(first(RATE_BOUNDS)), log(last(RATE_BOUNDS)), 17)
    end
    fit = fixed_K > 0 ? at_equilibrium(log(fixed_K)) :
          profile_minimum(at_equilibrium, log(first(EQUILIBRIUM_BOUNDS)),
                          log(last(EQUILIBRIUM_BOUNDS)), 17)
    merge(fit, (; gradient_inf=NaN, solver="bounded_profile"))
end

"""Fit the common k and Kx at a given concentration-law shape."""
function conditional_fit(law, shape, reference_loading, exposures, loadings, ratios,
                         observed, sqrtweights, fixed_K)
    scaled = [s*rate_scale(law, shape, w/reference_loading)
              for (s,w) in zip(exposures, loadings)]
    fixed = fixed_K > 0
    weighted_y = sqrtweights .* observed
    k(p) = exp(first(p))
    equilibrium(p) = fixed ? fixed_K : exp(p[2])
    predict(p) = sqrtweights .* [progress_conversion(k(p), equilibrium(p), R, s)
                                for (R,s) in zip(ratios, scaled)]
    # Reuse the previous weighted nonlinear least-squares routine; the
    # generic alternate starts only safeguard the two nuisance parameters.
    for (k0,K0) in ((0.12,4.0), (0.06,2.0), (0.24,10.0))
        initial = fixed ? [log(k0)] : log.([k0,K0])
        try
            fit = PentylPropanoate.least_squares(predict, weighted_y, initial)
            p = fit.parameters
            if first(RATE_BOUNDS) <= k(p) <= last(RATE_BOUNDS) &&
                first(EQUILIBRIUM_BOUNDS) <= equilibrium(p) <= last(EQUILIBRIUM_BOUNDS)
                return (; law, shape, k_ref=k(p), Kx=equilibrium(p),
                        weighted_sse=sum(abs2, fit.predicted-weighted_y),
                        gradient_inf=fit.gradient_inf, solver="least_squares")
            end
        catch err
            # Do not suppress a programming error; only retry numerical
            # convergence or a nonphysical trial during nonlinear fitting.
            if !(err isa ArgumentError || (err isa ErrorException &&
                    occursin("least-squares fit did not converge", err.msg)))
                rethrow()
            end
        end
    end
    bounded_conditional_fit(law, shape, scaled, ratios, observed, sqrtweights, fixed_K)
end

"""Profile the one extra loading parameter, retaining the nested boundary."""
function fit_loading(law, loading_pct, time_min, ratios, observed, weights;
                     reference_loading_pct=10.0, fixed_K=0.0, profile_points=25)
    n = length(observed)
    n >= 4 && n == length(loading_pct) == length(time_min) == length(ratios) == length(weights) &&
        all(x -> isfinite(x) && 0 <= x <= 1, observed) &&
        all(x -> isfinite(x) && x > 0, loading_pct) &&
        all(x -> isfinite(x) && x > 0, time_min) &&
        all(x -> isfinite(x) && x > 0, ratios) &&
        all(x -> isfinite(x) && x > 0, weights) &&
        isfinite(reference_loading_pct) && reference_loading_pct > 0 &&
        isfinite(fixed_K) && fixed_K >= 0 && profile_points >= 9 ||
        throw(ArgumentError("Use matching finite, physical post-start batch inputs."))
    # Materialise bridge arrays once. Every conditional fit is native
    # Julia, including the bounded profile for difficult trial shapes.
    loading_pct, time_min, ratios, observed, weights =
        (Float64.(collect(v)) for v in (loading_pct, time_min, ratios, observed, weights))
    exposures = loading_pct .* time_min .* (0.6 * PentylPropanoate.ACID_MOLAR_MASS_KG_MOL)
    sqrtweights = sqrt.(weights)
    fit_at(shape) = conditional_fit(law, shape, reference_loading_pct, exposures,
        loading_pct, ratios, observed, sqrtweights, fixed_K)
    law == "proportional" && return fit_at(1.0)
    hi = law == "power" ? 3.0 : law == "crowding" ? 1.0 :
        throw(ArgumentError("Unknown loading law: $law"))

    # First locate the minimum in a broad deterministic shape scan.
    # Refine a bracket around each local grid minimum, including a
    # possible optimum near an endpoint. Always compare exact null.
    null = fit_at(law == "power" ? 1.0 : 0.0)
    best = profile_minimum(fit_at, 0.0, hi, profile_points)
    best.weighted_sse < null.weighted_sse-1e-13 ? best : null
end

end
