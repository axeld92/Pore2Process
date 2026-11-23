module SharedLoadingEvidence

import ..CatalystLoading

export fit_shared, predict_series, fit_crossing_law, log_total_rate

"""Keep an unidentifiable apparent equilibrium visible rather than inventing a plateau."""
function conditional_with_limit(law,shape,reference,E,C,R,X,sqrtW,fixed_K)
    fit = CatalystLoading.conditional_fit(law,shape,reference,E,C,R,X,sqrtW,fixed_K)
    if fixed_K > 0 || fit.Kx < .999*last(CatalystLoading.EQUILIBRIUM_BOUNDS)
        return merge(fit,(; limit_probe_change=0.))
    end
    scaled = E .* [CatalystLoading.rate_scale(law,shape,c/reference) for c in C]
    probe(K) = CatalystLoading.bounded_conditional_fit(law,shape,scaled,R,X,sqrtW,K)
    # A fitted K approaching the wide numerical ceiling may mean that
    # the slow training curves do not identify a common plateau. Treat
    # the near-irreversible limit as a mathematical diagnostic, not a
    # measured K; verify that expanding the cap does not affect progress.
    coarse,fine = probe(1e6),probe(1e9)
    predict(f) = [CatalystLoading.progress_conversion(f.k_ref,f.Kx,r,e) for (r,e) in zip(R,scaled)]
    change = maximum(abs.(predict(coarse)-predict(fine)))
    change <= 1e-5 || error("Apparent equilibrium limit has not converged.")
    fine.weighted_sse <= fit.weighted_sse + 1e-10 ?
        merge(fine,(; solver="upper_equilibrium_limit",limit_probe_change=change)) :
        merge(fit,(; limit_probe_change=change))
end

"""Joint apparent rate-law shape, independent kinetic/equilibrium nuisance fits by feed."""
function fit_shared(law, concentrations, exposures, ratios, observed, weights, groups;
                    reference_concentration, fixed_shape=NaN, profile_points=25,
                    fixed_equilibria=Float64[])
    n = length(observed)
    n >= 4 && all(v -> length(v) == n, (concentrations, exposures, ratios, weights, groups)) &&
        all(v -> all(x -> isfinite(x) && x > 0, v), (concentrations, exposures, ratios, weights)) &&
        all(x -> isfinite(x) && 0 <= x <= 1, observed) && isfinite(reference_concentration) &&
        reference_concentration > 0 && profile_points >= 9 ||
        throw(ArgumentError("Expected matching physical post-start observations and positive reference concentration."))
    C, E, R, X, W = (Float64.(collect(v)) for v in (concentrations, exposures, ratios, observed, weights))
    g = Int.(collect(groups))
    series = sort(unique(g))
    indices = [findall(==(key), g) for key in series]
    all(i -> length(i) >= 4, indices) ||
        throw(ArgumentError("Each feed needs at least four post-start observations."))
    prescribed = isempty(fixed_equilibria) ? zeros(length(series)) : Float64.(collect(fixed_equilibria))
    length(prescribed) == length(series) && all(x->isfinite(x) && x>=0,prescribed) ||
        throw(ArgumentError("Provide one finite, nonnegative fixed K per sorted feed or none."))
    function at_shape(shape)
        fits = [conditional_with_limit(law, shape, reference_concentration,
            E[i], C[i], R[i], X[i], sqrt.(W[i]), K) for (i,K) in zip(indices,prescribed)]
        (; law, shape, series, fits, weighted_sse=sum(f.weighted_sse for f in fits))
    end
    law == "proportional" && return at_shape(1.0)
    isfinite(fixed_shape) && return at_shape(fixed_shape)
    hi = law == "power" ? 3.0 : law == "crowding" ? 1.0 :
        throw(ArgumentError("Unknown loading law $law"))
    null = at_shape(law == "power" ? 1.0 : 0.0)
    best = CatalystLoading.profile_minimum(at_shape, 0., hi, profile_points)
    best.weighted_sse < null.weighted_sse-1e-13 ? best : null
end

function predict_series(law, shape, k_ref, Kx, concentrations, exposures, ratios;
                        reference_concentration)
    length(concentrations) == length(exposures) == length(ratios) ||
        throw(ArgumentError("Prediction columns must have equal lengths."))
    [CatalystLoading.progress_conversion(k_ref*CatalystLoading.rate_scale(
        law, shape, C/reference_concentration), Kx, R, E)
     for (C,E,R) in zip(concentrations,exposures,ratios)]
end

"""Log of total rate relative to the positive concentration anchor."""
log_total_rate(law, shape, u) = log(u) + log(CatalystLoading.rate_scale(law,shape,u))

"""
Profile progress-time measurements with independent log-time intercepts.
One block = one study and matched conversion level. Different chemistries
need not share either the kinetics or the conversion level.
"""
function fit_crossing_law(law, relative_concentrations, times, weights, blocks;
                          profile_points=25)
    n = length(times)
    n >= 2 && all(v -> length(v) == n, (relative_concentrations,weights,blocks)) &&
        all(v -> all(x -> isfinite(x) && x > 0, v), (relative_concentrations,times,weights)) &&
        profile_points >= 9 ||
        throw(ArgumentError("Expected positive finite matched times, concentrations and weights."))
    u,t,w = (Float64.(collect(v)) for v in (relative_concentrations,times,weights))
    b = Int.(collect(blocks))
    ids = sort(unique(b))
    idx = [findall(==(key), b) for key in ids]
    all(i -> length(unique(u[i])) >= 2, idx) ||
        throw(ArgumentError("Each training progress block needs two distinct concentrations."))
    function at_shape(shape)
        logtime = log.(t) .+ [log_total_rate(law,shape,v) for v in u]
        intercepts = [sum(w[i].*logtime[i])/sum(w[i]) for i in idx]
        sse = sum(sum(w[i].*(logtime[i].-a).^2) for (i,a) in zip(idx,intercepts))
        (; law, shape, weighted_sse=sse, blocks=ids, intercepts)
    end
    law == "proportional" && return at_shape(1.0)
    hi = law == "power" ? 3.0 : law == "crowding" ? 1.0 :
        throw(ArgumentError("Unknown loading law $law"))
    null = at_shape(law == "power" ? 1.0 : 0.0)
    best = CatalystLoading.profile_minimum(at_shape,0.,hi,profile_points)
    best.weighted_sse < null.weighted_sse-1e-13 ? best : null
end

end
