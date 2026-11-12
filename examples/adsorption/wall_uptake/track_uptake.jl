module TrackUptake

using Random
include(joinpath(@__DIR__,"../carbon_contact/track_contacts.jl"))
using .TrackContacts: Contacts, FaceEvent, Flow, advect, diffusion_contacts

mutable struct FirstHit
    time::Float64
    before::Float64
    after::Float64
end

FirstHit() = FirstHit(NaN,NaN,NaN)

@inline function Base.push!(observer::FirstHit,event::FaceEvent)
    if event.marker_phase>0
        isnan(observer.before) && (observer.before=observer.time)
    else
        isnan(observer.after) && (observer.after=observer.time)
    end
    observer
end

"""Use the existing reflector with a first-hit clock instead of a face log."""
@inline function trace(flow,x,p,rng,dx,D,dt,stop,marker)
    first = NaN
    observer = FirstHit()
    stats = Contacts()
    sigma = sqrt(2D*dt)/dx
    half = dt/2
    for step in 1:round(Int,stop/dt)
        t = (step-1)*dt
        x,p,first,exit = advect(flow,x,p,first,t,half)
        isfinite(exit) && return first,exit-first,observer.before,observer.after-first
        observer.time = t+half
        for axis in 1:3
            x,p,first,stats = diffusion_contacts(flow,x,p,first,axis,sigma*randn(rng),
                t+half,dx,stats,observer,marker)
        end
        x,p,first,exit = advect(flow,x,p,first,t+half,half)
        isfinite(exit) && return first,exit-first,observer.before,observer.after-first
    end
    first,NaN,observer.before,observer.after-first
end

function sample(flow,dx,D,dt,stop,count,seed)
    inlet = findall(i->flow.before[1,i]==0 && flow.low[1,i]>0,1:flow.first_count)
    isempty(inlet) && error("No forward inlet")
    weights = cumsum(flow.low[1,inlet])
    result = Matrix{Float64}(undef,4,count)
    Threads.@threads for i in 1:count
        rng = Xoshiro(seed+i)
        p = inlet[searchsortedfirst(weights,rand(rng)*weights[end])]
        x = (0.,rand(rng),rand(rng))
        result[:,i] .= trace(flow,x,p,rng,dx,D,dt,stop,i)
    end
    result
end

function main(args)
    length(args)==8 || error("Use compact output.f64 dx D dt stop count seed")
    source,output = args[1:2]
    dx,D,dt,stop = parse.(Float64,args[3:6])
    count,seed = parse.(Int,args[7:8])
    Base.ENDIAN_BOM==0x04030201 || error("Expected little-endian storage")
    dx>0 && D>0 && 0<dt<=stop && count>0 || error("Invalid wall-uptake settings")
    read_grid(name,::Type{T}) where T = reshape(reinterpret(T,read(joinpath(source,name))),3,:)
    before = read_grid("before.i32",Int32)
    after = read_grid("after.i32",Int32)
    low = read_grid("incoming.f64",Float64) ./ dx^3
    high = read_grid("outgoing.f64",Float64) ./ dx^3
    flow = Flow(before,after,low,high)
    data = sample(flow,dx,D,dt,stop,count,seed)
    open(io->write(io,data),output,"w")
    println("Completed ",sum(isfinite,view(data,2,:))," pairs, first- and second-period captures ",
        sum(isfinite,view(data,3,:))," / ",sum(isfinite,view(data,4,:)))
end

if abspath(PROGRAM_FILE)==@__FILE__
    main(ARGS)
end

end
