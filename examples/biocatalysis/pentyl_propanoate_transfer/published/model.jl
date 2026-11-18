module LorenzoPackedBed

using LinearAlgebra: Tridiagonal, norm

export PingPong, batch_prediction, packed_bed, disappearance_rate

"""Parameters in the order of Table S3, using its tabulated numeric basis."""
struct PingPong
    k1ref::Float64
    k1exp::Float64
    k2ref::Float64
    k2exp::Float64
    k3::Float64
    k4::Float64
    k5::Float64
    k6::Float64
    k7::Float64
    k8::Float64
    k9::Float64
    k10::Float64
    Kdref::Float64
    Kd1::Float64
    Kd2::Float64
    Kdexp::Float64
    Tref::Float64
end

PingPong(values::AbstractVector) = PingPong(values...)

function rate_terms(p, temperature_C, acid0, alcohol0, temperature_dependent_ester_term)
    inverse_shift = inv(temperature_C + 273.15) - inv(p.Tref)
    k1 = p.k1ref * exp(-p.k1exp * inverse_shift)
    k2 = p.k2ref * exp(-p.k2exp * inverse_shift)
    kd = p.Kdref * acid0^p.Kd1 / (alcohol0 / acid0)^p.Kd2 * exp(-p.Kdexp * inverse_shift)
    single_ester_ratio = temperature_dependent_ester_term ? k2 / k1 : p.k2ref / p.k1ref
    (; k1, k2, kd, single_ester = single_ester_ratio * p.k3 * p.k4 / p.k8,
       acid_ester = k2 * p.k4 * p.k6 / (k1 * p.k8),
       alcohol_water = p.k7 * p.k8 / p.k4,
       water_ester = k2 * p.k4 * (p.k9 + kd) / (k1 * p.k8) -
                     k2 * (p.k5 + kd) / k1 + p.k8 * p.k10 / p.k4)
end

function rate(p, terms, acid0, alcohol0, catalyst_g_L, conversion)
    # Eq. S1 does not give the bed species profiles. Close them with a
    # dry feed, single stoichiometric reaction and constant liquid flow.
    ester = water = acid0 * conversion
    acid = acid0 - ester
    alcohol = alcohol0 - ester
    den = p.k3 * acid + p.k4 * alcohol + (p.k5 + terms.kd) * acid * alcohol +
          p.k6 * acid^2 + p.k7 * alcohol^2 + p.k8 * water + terms.single_ester * ester +
          (p.k9 + terms.kd) * acid * water + p.k10 * alcohol * ester +
          terms.acid_ester * acid * ester + terms.alcohol_water * alcohol * water +
          terms.water_ester * water * ester
    catalyst_g_L * (terms.k1 * acid * alcohol - terms.k2 * water * ester) / den
end

function validate_rate_inputs(temperature_C, acid0, alcohol0, catalyst_g_L)
    all(isfinite, (temperature_C, acid0, alcohol0, catalyst_g_L)) &&
        temperature_C > -273.15 && alcohol0 >= acid0 > 0 && catalyst_g_L >= 0 ||
        throw(ArgumentError("Use a physical temperature, alcohol-rich feed and catalyst concentration."))
end

"""Positive acid disappearance from Eq. S1, with g/L and minute as a conditional basis."""
function disappearance_rate(p, temperature_C, conversion, catalyst_g_L;
                            acid0=3.3, alcohol0=6.6, temperature_dependent_ester_term=false)
    validate_rate_inputs(temperature_C, acid0, alcohol0, catalyst_g_L)
    0 <= conversion <= 1 || throw(ArgumentError("Conversion must lie in [0, 1]."))
    terms = rate_terms(p, temperature_C, acid0, alcohol0, temperature_dependent_ester_term)
    rate(p, terms, acid0, alcohol0, catalyst_g_L, conversion)
end

"""Integrate the printed rate law in batch as a check on its numerical basis."""
function batch_prediction(p, temperature_C, time_min, catalyst_g_L;
                          acid0=3.3, alcohol0=6.6, max_step_min=0.05)
    validate_rate_inputs(temperature_C, acid0, alcohol0, catalyst_g_L)
    isfinite(time_min) && isfinite(max_step_min) && time_min >= 0 && max_step_min > 0 ||
        throw(ArgumentError("Use a nonnegative time and a finite positive integration step."))
    terms = rate_terms(p, temperature_C, acid0, alcohol0, false)
    dxdt(x) = rate(p, terms, acid0, alcohol0, catalyst_g_L, x) / acid0
    steps = max(1, ceil(Int, time_min / max_step_min))
    dt = time_min / steps
    x = 0.0
    for _ in 1:steps
        a = dxdt(x)
        b = dxdt(x + 0.5dt * a)
        c = dxdt(x + 0.5dt * b)
        d = dxdt(x + dt * c)
        x += dt * (a + 2b + 2c + d) / 6
    end
    x
end

"""Solve the acid-sink interpretation of Eq. 2, using the inlet in Eq. 4 by default."""
function packed_bed(p, temperature_C, residence_min, catalyst_g_L, effectiveness, peclet;
                    acid0=3.3, alcohol0=6.6, cells=1024, inlet_flux=false,
                    temperature_dependent_ester_term=false)
    validate_rate_inputs(temperature_C, acid0, alcohol0, catalyst_g_L)
    all(isfinite, (residence_min, effectiveness, peclet)) && residence_min > 0 &&
        0 <= effectiveness <= 1 && peclet > 0 && cells >= 2 ||
        throw(ArgumentError("Use physical residence time, effectiveness, Péclet number and grid size."))
    terms = rate_terms(p, temperature_C, acid0, alcohol0, temperature_dependent_ester_term)
    factor = residence_min * effectiveness / (acid0 * cells)
    abstract_rate(x) = rate(p, terms, acid0, alcohol0, catalyst_g_L, x)
    cell_source(x) = factor * abstract_rate(x)
    # The numerator crosses zero once in the physical region. Never step
    # past that equilibrium while solving the acid-consumption problem.
    lo, hi = 0.0, 1.0
    for _ in 1:52
        mid = (lo + hi) / 2
        if terms.k1 * (1 - mid) * (alcohol0 / acid0 - mid) > terms.k2 * mid^2
            lo = mid
        else
            hi = mid
        end
    end
    equilibrium = (lo + hi) / 2
    # Exponentially fitted face flux J = X - (1/Pe)dX/d(z/L).
    down = inv(expm1(peclet / cells))
    inlet = inlet_flux ? 0.0 : inv(expm1(peclet / (2cells)))
    diagonal = fill(1 + 2down, cells)
    diagonal[1] = 1 + down + inlet
    diagonal[end] = 1 + down
    lower, upper = fill(-(1 + down), cells - 1), fill(-down, cells - 1)
    function residual(x)
        internal = @views x[1:end-1] .+ down .* (x[1:end-1] .- x[2:end])
        diff(vcat(-inlet * x[1], internal, x[end])) - cell_source.(x)
    end
    x = zeros(cells)
    for _ in 1:40
        r = residual(x)
        if norm(r, 1) <= 5e-10
            source = sum(cell_source, x)
            return (; conversion=x[end], inlet_backdiffusion_fraction=inlet * x[1],
                    integrated_source_fraction=source, balance_error=abs(x[end] + inlet * x[1] - source))
        end
        δ = 1e-6
        slope(y) = (cell_source(y + δ) - cell_source(y - δ)) / (2δ)
        step = Tridiagonal(lower, diagonal - slope.(x), upper) \ r
        for halve in 0:20
            proposal = clamp.(x - step / 2.0^halve, 0.0, equilibrium)
            if norm(residual(proposal), 1) < norm(r, 1)
                x = proposal
                break
            end
        end
    end
    error("The finite-volume reaction and flux balance did not converge.")
end

end
