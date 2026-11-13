using Test, LinearAlgebra
include("model.jl")
using .DisorderedExchange: calibrate, propagate

@testset "Inert passage, wall transfer, and held-out propagation" begin
    n = 100
    localtime = collect(1.:n)
    arrivals = hcat(localtime,2localtime)
    hits = fill(NaN,n,2)
    model = calibrate(arrivals,hits)
    @test model.counts == 10 .* Matrix{Int}(I,10,10)
    @test model.wallfree == model.transition
    @test sum(model.flux) ≈ 1
    inert = propagate(model,8)
    @test all(inert[:,1:2] .≈ 1.)
    @test all(isapprox.(inert[:,3:4],0.;atol=1e-15))
    @test inert[1,5] ≈ 3*sum(localtime)/n
    @test inert[:,7] ≈ ones(6)

    # Half of the slowest decile reaches the wall in each period.
    hits[96:100,2] = arrivals[96:100,1]
    model = calibrate(arrivals,hits)
    @test model.wallfree[10,10] == .5
    @test model.states[10,4] == .05
    result = propagate(model,8)
    @test result[:,1] ≈ [.9+.05*.5^i for i in 1:6]
    @test result[:,2] ≈ ones(6)
    @test result[1,3] ≈ .025/.95
    @test result[:,5:6] == inert[:,5:6]

    # A directed exchange checks that donor rows feed receiver columns.
    reordered = hcat(localtime,localtime+circshift(localtime,10))
    model = calibrate(reordered,fill(NaN,n,2))
    @test model.transition[1,10] == 1.
    @test model.transition[2,1] == 1.
    @test propagate(model,3)[1,7] ≈ 7/9
    @test_throws ArgumentError calibrate(arrivals[:,1:1],hits[:,1:1])
end
