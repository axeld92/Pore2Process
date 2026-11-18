using Test

include("model.jl")
using .LorenzoPackedBed

# k4 times alcohol is the only material denominator term for this test.
# The tiny positive k8 keeps the original rate-law ratios defined.
first_order = PingPong([0.7, 0, 0, 0, 0, 1, 0, 0, 0, 1e-100, 0, 0, 0, 0, 0, 0, 313.15])

function exact_first_order_outlet(peclet, damkohler, flux_inlet)
    λ = (peclet - sqrt(peclet^2 + 4peclet * damkohler)) / 2
    μ = peclet - λ
    b_over_a = -λ / μ * exp(λ)
    correction = b_over_a * exp(-μ)
    a = flux_inlet ? inv(1 - λ / peclet + correction * (1 - μ / peclet)) : inv(1 + correction)
    1 - a * (exp(λ) + b_over_a)
end

@testset "Rate and finite-volume balances" begin
    dose = 1.2
    residence = 2.0
    eta = 0.8
    peclet = 19.5
    @test disappearance_rate(first_order, 40, 0.2, dose) ≈ 0.7 * dose * 3.3 * 0.8
    @test batch_prediction(first_order, 40, residence, dose) ≈ 1 - exp(-0.7dose * residence) atol=1e-8
    @test packed_bed(first_order, 40, residence, 0, eta, peclet).conversion == 0.0
    for flux in (false, true)
        expected = exact_first_order_outlet(peclet, 0.7dose * residence * eta, flux)
        coarse = packed_bed(first_order, 40, residence, dose, eta, peclet; inlet_flux=flux)
        outcome = packed_bed(first_order, 40, residence, dose, eta, peclet; inlet_flux=flux, cells=2048)
        @test outcome.conversion ≈ expected atol=2.5e-7
        @test abs(outcome.conversion - expected) < 0.3 * abs(coarse.conversion - expected)
        @test outcome.balance_error < 1e-9
        @test flux == iszero(outcome.inlet_backdiffusion_fraction)
    end
end
