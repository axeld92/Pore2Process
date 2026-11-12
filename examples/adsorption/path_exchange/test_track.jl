using Random, Statistics, Test
include("track.jl")
using .Track: Flow, advect, crossing, diffusion_axis, drift, trace

@testset "Marker travel and reflection" begin
    @test crossing(.25,1.,1.,0.) == (.75,1)
    @test crossing(.25,-1.,-1.,0.) == (.25,-1)
    @test crossing(.25,1.,0.,-1.) == (Inf,0)
    dt,side = crossing(.25,1.,2.,1.)
    @test side==1
    @test dt ≈ log(2/1.25)
    @test drift(.25,1.,2.,1.,dt) ≈ 1.
    @test drift(.25,1.,0.,-1.,log(2.)) ≈ .625
    @test first(crossing(.75,1.,1e-17,-1.)) ≈ log(.25/1e-17)
    before = Int32[0 1; 0 0; 0 0]
    after = Int32[2 0; 0 0; 0 0]
    normal = [2. 2.; 0. 0.; 0. 0.]
    flow = Flow(before,after,normal,normal)
    @test trace(flow,(0.,.5,.5),1,Xoshiro(42),1.,0.,.1,2.) == (.5,.5)
    @test advect(flow,(0.,.5,.5),1,NaN,0.,.75)[1:3] == ((.5,.5,.5),2,.5)
    @test diffusion_axis(flow,(.25,.75,.5),1,NaN,2,.5,.125)[1:2] == ((.25,.75,.5),1)
    @test diffusion_axis(flow,(.75,.5,.5),1,NaN,1,.5,.125) == ((.25,.5,.5),2,.125)
    @test diffusion_axis(flow,(.25,.5,.5),2,.125,1,-.5,.25) == ((.75,.5,.5),1,.125)
    rng = Xoshiro(100)
    reflected = [diffusion_axis(flow,(.5,rand(rng),.5),1,NaN,2,.4*randn(rng),.1)[1][2] for _ in 1:50_000]
    @test abs(mean(reflected)-.5)<.005
    @test abs(var(reflected)-1/12)<.002
end
