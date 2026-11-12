using Test
include("track_uptake.jl")
using .TrackUptake: FirstHit, sample
using .TrackUptake.TrackContacts: Flow, FaceEvent, sample_contacts

@testset "First absorbing contact and inert reference" begin
    ob = FirstHit()
    ob.time = .02
    push!(ob,FaceEvent(1,4,1.))
    ob.time = .04
    push!(ob,FaceEvent(1,4,1.))
    @test ob.before==.02
    @test isnan(ob.after)
    push!(ob,FaceEvent(-1,4,1.))
    @test ob.after==.04
    before = Int32[0 1;0 0;0 0]
    after = Int32[2 0;0 0;0 0]
    normal = [2. 2.;0. 0.;0. 0.]
    flow = Flow(before,after,normal,normal)
    ref,_ = sample_contacts(flow,1.,.04,.01,3.,1000,42)
    hit = sample(flow,1.,.04,.01,3.,1000,42)
    @test hit[1:2,:]==ref[1:2,:]
    @test isfinite.(hit[3:4,:])==(ref[3:4,:].>0)
    for i in 1:2
        good = isfinite.(hit[i+2,:])
        @test any(good)
        @test all(0 .<= hit[i+2,good] .<= hit[i,good])
    end
end
