using Test
include("model.jl")
using .ParticleTransport

@testset "Sphere and active shell with fixed total intrinsic rate" begin
    @test effectiveness(0) == 1
    @test effectiveness(1e-6, 0.7) ≈ 1
    @test effectiveness(3) ≈ (3 / tanh(3) - 1) / 3
    @test effectiveness(modulus_for(0.869)) ≈ 0.869
    @test relative_effectiveness(2, 1, 1; core=0.7) ≈ 1

    # Independently integrate the concentration and cumulative radial flux
    # across the shell and normalize by the computed outer concentration.
    source(x, y) = (derivative=y[2] / x^2, release=4 / (1-0.7^3) * x^2 * y[1])
    x, y = 0.7, [1.0, 0.0]
    step = (1-x) / 1000
    for _ in 1:1000
        a = collect(source(x, y))
        b = collect(source(x + step/2, y + step*a/2))
        c = collect(source(x + step/2, y + step*b/2))
        d = collect(source(x + step, y + step*c))
        y += step * (a + 2b + 2c + d) / 6
        x += step
    end
    @test effectiveness(2, 0.7) ≈ 3 * y[2] / (4 * y[1]) atol=1e-12
    @test effectiveness(2, 0.7) > effectiveness(2)
    @test relative_effectiveness(2, 0.55, 0.83; core=0.7) >
        relative_effectiveness(2, 0.55, 0.83)

    lower_limit = sqrt(0.55 * 0.83)
    for phi in (0.1, 1.0, 4.0, 10.0, 100.0, 10_000.0)
        @test lower_limit < relative_effectiveness(phi, 0.55, 0.83) < 1
    end
    @test_throws ArgumentError effectiveness(3, 1.0)
    @test_throws ArgumentError relative_effectiveness(2, 0, 1)
end
