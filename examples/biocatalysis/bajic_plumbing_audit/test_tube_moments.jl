using Test
using LinearAlgebra
include("tube_moments.jl")
using .TubeMoments: spectrum, moments

const radius=0.25e-3
const diffusivity=7.1e-10
const scale=radius^2/diffusivity

@testset "Conservation and uniform plug" begin
    p=spectrum(32;profile=:plug)
    for t in (0.,2.,10.,100.)
        m=moments(p,t,radius,diffusivity)
        @test m.mean_s ≈ t rtol=1e-13
        @test m.variance_s2 <= 1e-22
    end
end

@testset "Circular tube independent analytical limits" begin
    p=spectrum(256)
    @test sum(p.weight) ≈ 1 rtol=1e-13
    # Flux-weighted stationary injection gives volume/flow exactly.
    @test p.mean_factor ≈ 1 rtol=1e-13
    @test p.slope_dimensionless ≈ 1/24 rtol=1e-3
    @test moments(p,0.,radius,diffusivity).variance_s2 == 0.
    # Short-distance limit at fixed annular grid: covariance(0)*m².
    tiny=1e-10
    meaninv=p.mean_factor
    direct=dot(p.weight,(1 ./p.relative_velocity .-meaninv).^2)
    @test moments(p,tiny,radius,diffusivity).variance_s2/tiny^2 ≈ direct rtol=1e-3
    # Long-distance slope tends to the closed-form continuum result.
    later=moments(p,2000.,radius,diffusivity)
    earlier=moments(p,1900.,radius,diffusivity)
    @test (later.variance_s2-earlier.variance_s2)/100 ≈ scale/24 rtol=1e-3
    # Finite segment/extension: covariance is nonnegative; slope never
    # exceeds the same spectrum's infinite-length limit.
    for base in (0.,10.,100.), added in (2.,10.5)
        first=moments(p,base,radius,diffusivity)
        last=moments(p,base+added,radius,diffusivity)
        @test 0 < last.variance_s2-first.variance_s2 < last.asymptotic_slope_s*added
    end
end

@testset "Radial convergence" begin
    a=spectrum(256); b=spectrum(512)
    for t in (.5,2.,10.5,100.)
        va=moments(a,t,radius,diffusivity).variance_s2
        vb=moments(b,t,radius,diffusivity).variance_s2
        @test isapprox(va,vb;rtol=.01)
    end
    @test_throws ArgumentError moments(b,-1.,radius,diffusivity)
end
