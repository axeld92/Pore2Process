"""Carman–Kozeny packed-primary-particle reference with an empirical multiplier.

Particle diameter is in m; porosity and multiplier are dimensionless; k is m².
This reduced closure does not identify a connected pore geometry. The caller
must check applicability above percolation and for the preparation family.
"""
function support_permeability(particle_diameter_m, porosity; multiplier=1.0)
    d = positive("primary particle diameter_m", particle_diameter_m)
    e = positive("porosity fraction", porosity)
    e < 1 || throw(ArgumentError("Porosity fraction must be strictly between 0 and 1."))
    b = positive("empirical permeability multiplier", multiplier)
    # Logs avoid premature underflow/overflow in d² and (1-e)².
    positive("predicted permeability_m2", exp(log(b)+2log(d)+3log(e)-log(180.)-2log1p(-e)))
end

"""Fit ONLY the supplied calibration vectors, minimising equal-weight squared
log residuals. The analytic minimiser is log(b)=mean(log(k)-log(k_CK)).
No uncertainty model or confidence interval is implied.
"""
function fit_support_multiplier(diameter_m::AbstractVector, porosity::AbstractVector,
                                measured_k_m2::AbstractVector)
    n = length(measured_k_m2)
    n > 0 && length(diameter_m) == length(porosity) == n ||
        throw(ArgumentError("Calibration diameter, porosity and measured-k vectors must have equal, nonzero length."))
    baseline = support_permeability.(diameter_m, porosity)
    measured = [positive("measured permeability_m2", x) for x in measured_k_m2]
    log_ratios = log.(measured) .- log.(baseline)
    log_b = sum(log_ratios)/n
    b = positive("fitted multiplier", exp(log_b))
    predictions = support_permeability.(diameter_m, porosity; multiplier=b)
    (; multiplier=b, log_multiplier=log_b, calibration_count=n,
       calibration_log_rmse=sqrt(sum(abs2, log_ratios .- log_b)/n), predictions_m2=predictions)
end

"""Prediction diagnostics for explicit vectors (no fitting). Relative errors use
the measured value as denominator. All errors are dimensionless.
"""
function permeability_errors(predicted_k_m2::AbstractVector, measured_k_m2::AbstractVector)
    n = length(measured_k_m2)
    n > 0 && length(predicted_k_m2) == n ||
        throw(ArgumentError("Prediction and measured-k vectors must have equal, nonzero length."))
    predicted = [positive("predicted permeability_m2", x) for x in predicted_k_m2]
    measured = [positive("measured permeability_m2", x) for x in measured_k_m2]
    relative = predicted ./ measured .- 1
    all(isfinite,relative) || throw(ArgumentError("Relative permeability error overflows; use less extreme inputs."))
    logs = log.(predicted) .- log.(measured)
    (; count=n, mean_absolute_relative_error=sum(abs,relative)/n,
       max_absolute_relative_error=maximum(abs,relative), mean_relative_error=sum(relative)/n,
       log_rmse=sqrt(sum(abs2,logs)/n))
end
