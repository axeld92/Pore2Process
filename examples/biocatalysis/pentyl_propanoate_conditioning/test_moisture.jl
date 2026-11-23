using Test
include("moisture.jl")
using .BeadMoisture

@testset "Sphere with fixed surface: independent diffusion solution" begin
    R, eps, tau, D = 296.5e-6, 0.40, 5.3, 8e-10
    exact(t) = 6/pi^2 * sum(exp(-n^2*pi^2*D*t/(tau*R^2))/n^2 for n in 1:512)
    result = relax_water(R, eps, tau, 592., D, Inf, 0., 1., 0., 1.;
                         cells=128, dt=0.25, duration=600., sample_step=60.)
    reverse = relax_water(R, eps, tau, 592., D, Inf, 0., 1., 1., 0.;
                          cells=128, dt=0.25, duration=600., sample_step=60.)
    @test result.remaining[1] ≈ 1
    @test maximum(abs.(result.remaining[2:end]-exact.(result.times_s[2:end]))) <= 0.001
    @test result.remaining ≈ reverse.remaining atol=1e-9
    @test result.remaining ≈ result.pore_remaining atol=1e-10
    @test result.max_relative_balance_error <= 1e-7
    @test all(result.remaining_shell .<= result.remaining .+ 1e-10)
end

@testset "Reversible storage and film balances" begin
    solve(qmax, c0, cb, cells, dt) = relax_water(
        296.5e-6, 0.22, 5.3, 592., 8e-10, 1e-5, qmax, 0.42, c0, cb;
        cells, dt, duration=1800., sample_step=60.)
    plain = solve(0., 0., 0.42, 64, 2.)
    for (c0, cb) in ((0., 0.42), (0.42, 0.))
        coarse = solve(0.045, c0, cb, 64, 2.)
        fine = solve(0.045, c0, cb, 128, 1.)
        @test maximum(abs.(coarse.remaining-fine.remaining)) <= 0.005
        @test fine.max_relative_balance_error <= 1e-7
        @test all(diff(fine.remaining) .<= 1e-10)
        @test all(fine.remaining .>= -1e-10)
        @test all(fine.remaining_shell .<= fine.remaining .+ 1e-10)
        @test fine.remaining[end] > plain.remaining[end]
    end
    @test_throws ArgumentError relax_water(0., 0.40, 5.3, 592., 8e-10, 1e-5,
                                           0.025, 0.42, 0., 1.)
end
