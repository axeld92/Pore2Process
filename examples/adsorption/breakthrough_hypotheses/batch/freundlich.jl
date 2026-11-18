module PublishedFreundlich

using ..FilmParticle

struct Isotherm
    K::Float64
    exponent::Float64
end

# Alhamed Table 1 gives q in mg/g and C in mg/L.
const TABLE1 = Isotherm(24.6 * 1000^(1/2.74) / 1000, 1/2.74)

@inline loading(c, law) = law.K * c^law.exponent

function FilmParticle.isotherm(_, _, law::Isotherm)
    qe(c) = loading(c, law)
    # Keep the derivative finite at zero for the existing Newton bookkeeping.
    # The equilibrium value itself is evaluated at the exact concentration.
    slope(c) = law.K * law.exponent * max(c, floatmin(Float64))^(law.exponent - 1)
    qe, slope
end

"""Invert C + beta*q(C) = rhs without differentiating the power at C=0."""
@inline function surface_concentration(rhs, beta, law)
    rhs > 0 && beta > 0 || return rhs
    # Setting x = C^exponent gives an increasing convex function of x.
    alpha = beta * law.K
    n = inv(law.exponent)
    x = min(rhs^law.exponent, rhs/alpha)
    for _ in 1:18
        power = x^(n-1)
        c = x*power
        # Newton starts at an upper bound. If its mapped concentration
        # already underflows, the solution is zero at this precision.
        c == 0 && return 0.
        residual = muladd(alpha, x, c) - rhs
        abs(residual) <= 8eps(rhs) && return c
        x -= residual/(alpha+n*power)
    end
    error("Freundlich surface concentration root failed for rhs=$rhs and beta=$beta")
end

end
