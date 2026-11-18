module EffectiveEquilibrium

using ..FilmParticle
using ..PublishedFreundlich

struct Blend
    weight::Float64
    function Blend(weight)
        0 < weight < 1 || throw(ArgumentError("A blend requires a weight between zero and one"))
        new(weight)
    end
end

function FilmParticle.isotherm(qmax, b, law::Blend)
    langmuir, dlangmuir = FilmParticle.isotherm(qmax, b, nothing)
    freundlich, dfreundlich = FilmParticle.isotherm(qmax, b, PublishedFreundlich.TABLE1)
    q(c) = (1-law.weight)*langmuir(c) + law.weight*freundlich(c)
    slope(c) = (1-law.weight)*dlangmuir(c) + law.weight*dfreundlich(c)
    q, slope
end

"""Solve C + beta*q(C) = rhs for an intermediate mixture of the source laws."""
@inline function surface_concentration(rhs, beta, qmax, b, langmuir_upper, law)
    rhs > 0 && beta > 0 || return rhs
    source = PublishedFreundlich.TABLE1
    alpha = beta*law.weight*source.K
    scaled_capacity = beta*(1-law.weight)*qmax
    n = inv(source.exponent)
    lo = 0.
    hi = min(langmuir_upper^source.exponent, rhs/alpha)
    x = hi
    for _ in 1:40
        power = x^(n-1)
        c = x*power
        c == 0 && return 0.
        bc = b*c
        denom = 1+bc
        residual = muladd(alpha, x, c + scaled_capacity*bc/denom) - rhs
        abs(residual) <= 16eps(rhs) && return c
        if residual > 0
            hi = x
        else
            lo = x
        end
        derivative = muladd(n*power, 1 + scaled_capacity*b/denom^2, alpha)
        candidate = x-residual/derivative
        x = lo < candidate < hi ? candidate : (lo+hi)/2
    end
    error("Mixed-equilibrium surface concentration root failed for rhs=$rhs and beta=$beta")
end

end
