using Test
include("track.jl")
include("../wall_uptake/track_uptake.jl")
using .TrackRepeats: Flow, sample, tile

@testset "Repeat connections, first passage and wall clocks" begin
    before = Int32[0 1;0 0;0 0]
    after = Int32[2 0;0 0;0 0]
    v = [2. 2.;0. 0.;0. 0.]
    flow = Flow(before,after,v,v)
    tiled = tile(flow,4)
    @test tiled[1].before[1,:]==[0,1,2,3]
    @test tiled[1].after[1,:]==[2,3,4,0]
    inert = sample(flow,1.,0.,.01,4.,10,12,4,1)
    @test all(inert[1:4,:] .≈ [.5,1.,1.5,2.])
    @test all(isnan,inert[5:8,:])

    original = TrackUptake.sample(TrackUptake.Flow(before,after,v,v),1.,.04,.01,3.,1000,42)
    observed = sample(flow,1.,.04,.01,3.,1000,42,2,1)
    @test observed[1,:]==original[1,:]
    @test observed[2,:].-observed[1,:]==original[2,:]
    @test isequal(observed[3,:],original[3,:])
    @test isequal(observed[4,:].-observed[1,:],original[4,:])
end
