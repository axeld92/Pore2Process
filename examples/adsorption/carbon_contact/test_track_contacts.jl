using Test
include("track_contacts.jl")
using .TrackContacts: Contacts, FaceEvent, diffusion_contacts, sample_contacts
using .TrackContacts.Track: Flow, sample

@testset "Carbon face classification and chronological contact" begin
    before = Int32[0 1; 0 0; 0 0]
    after = Int32[2 0; 0 0; 0 0]
    normal = [2. 2.; 0. 0.; 0. 0.]
    flow = Flow(before,after,normal,normal)
    events = FaceEvent[]
    x,p,first,stats = diffusion_contacts(flow,(.25,.75,.5),1,NaN,2,.5,.125,2.,Contacts(),events,7)
    @test (x,p) == ((.25,.75,.5),1)
    @test stats.hits_before==1 && stats.hits_after==0 && stats.reflection_m_before==1.
    @test events==[FaceEvent(7,4,1.)]
    _,_,_,stats = diffusion_contacts(flow,(.25,.75,.5),1,NaN,1,-.5,.125,2.,stats,events,7)
    @test stats.hits_before==1 && stats.external_hits==1 && length(events)==1
    _,_,_,stats = diffusion_contacts(flow,(.75,.75,.5),2,.1,1,.5,.125,2.,stats,events,7)
    @test stats.hits_after==0 && stats.external_hits==2 && length(events)==1
    x,p,first,stats = diffusion_contacts(flow,(.75,.5,.5),1,NaN,1,.5,.125,2.,stats,events,7)
    @test (x,p,first) == ((.25,.5,.5),2,.125)
    x,p,first,stats = diffusion_contacts(flow,(.25,.5,.5),2,first,1,-.5,.25,2.,stats,events,7)
    @test (x,p,first) == ((.75,.5,.5),1,.125)
    _,_,_,stats = diffusion_contacts(flow,x,p,first,3,-1.,.25,2.,stats,events,7)
    @test stats.hits_after==1 && stats.backtracked_hits==1 && stats.reflection_m_after==2.
    @test events[end]==FaceEvent(-7,5,2.)
    _,_,_,stats = diffusion_contacts(flow,(.25,.5,.75),2,first,3,.5,.25,2.,stats,events,7)
    @test stats.hits_after==2 && stats.backtracked_hits==1 && stats.reflection_m_after==3.
    @test events[end]==FaceEvent(-7,6,1.)
    @test sizeof(FaceEvent)==16
    @test collect(reinterpret(UInt32,[FaceEvent(7,4,1.)]))==[7,4,0,0x3ff00000]
end

@testset "Observer preserves archived marker dynamics" begin
    before = Int32[0 1; 0 0; 0 0]
    after = Int32[2 0; 0 0; 0 0]
    normal = [2. 2.; 0. 0.; 0. 0.]
    flow = Flow(before,after,normal,normal)
    for d in (0., .04)
        expected = sample(flow,1.,d,.01,3.,1000,42)
        data,events = sample_contacts(flow,1.,d,.01,3.,1000,42)
        @test isequal(expected,data[1:2,:])
        @test sum(length,events)==sum(data[3:4,:])
        positive = sum(event.reflection_m for group in events for event in group if event.marker_phase>0;init=0.)
        negative = sum(event.reflection_m for group in events for event in group if event.marker_phase<0;init=0.)
        @test positive ≈ sum(data[5,:])
        @test negative ≈ sum(data[6,:])
    end
end
