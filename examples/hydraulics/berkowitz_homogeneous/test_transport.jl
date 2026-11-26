using Test, Statistics
include("first_passage.jl")
using .GlassPassage
using .GlassPassage.Track: Flow

@testset "Periodic flux, solid reflection and first passage" begin
    before = Int32[2 1; 0 0; 0 0]
    after = copy(before)
    rates = Float64[1 1; 0 0; 0 0]
    grid = Flow(before,after,copy(rates),copy(rates))
    planes = [3,7]
    times = GlassPassage.arrivals(grid,[1],planes,1.,0.,.37,9.,32,10)
    @test maximum(abs,times.-planes)<1e-12
    tracer = zeros(1)
    x,p,axial,next = GlassPassage.diffuse(grid,(.2,.3,.9),1,0,3,.8,0.,[2],tracer,1)
    @test x[3] ≈ .3
    @test (p,axial,next) == (1,0,1)
    x,p,axial,next = GlassPassage.diffuse(grid,x,p,axial,1,-.6,0.,[2],tracer,next)
    @test (p,axial,next) == (2,-1,1)
    @test x[1] ≈ .6
    # Infinite homogeneous tube with constant drift and molecular diffusion.
    # The exact first-passage moments are L/u and 2DL/u^3.
    D, L = .15, 24
    times = vec(GlassPassage.arrivals(grid,[1],[L],1.,D,.05,70.,4096,130))
    @test all(isfinite,times)
    @test abs(mean(times)/L-1)<.01
    @test abs(var(times)/(2D*L)-1)<.08
end
