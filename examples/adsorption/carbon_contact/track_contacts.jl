module TrackContacts

using Random
include(joinpath(@__DIR__, "../path_exchange/track.jl"))
using .Track: Flow, advect, destination

struct Contacts
    hits_before::Int
    hits_after::Int
    reflection_m_before::Float64
    reflection_m_after::Float64
    external_hits::Int
    backtracked_hits::Int
end

Contacts() = Contacts(0, 0, 0., 0., 0, 0)

struct FaceEvent
    # Sign records chronological phase: + before first axial seam crossing, - after.
    marker_phase::Int32
    # 1-based face of the period-wrapped original liquid voxel.
    face_id::UInt32
    reflection_m::Float64
end

@inline function external_end(flow, p, axis, side)
    axis == 1 && (side < 0 ? p <= flow.first_count && flow.low[1,p] > 0 :
                           p > flow.first_count && flow.high[1,p] > 0)
end

@inline function record_reflection(stats, events, flow, p, axis, side, first, marker, distance)
    if external_end(flow,p,axis,side)
        return Contacts(stats.hits_before,stats.hits_after,stats.reflection_m_before,
                        stats.reflection_m_after,stats.external_hits+1,stats.backtracked_hits)
    end
    wrapped = p <= flow.first_count ? p : p-flow.first_count
    face = 6*(wrapped-1) + 2*(axis-1) + (side>0 ? 2 : 1)
    if isnan(first)
        push!(events,FaceEvent(marker,face,distance))
        return Contacts(stats.hits_before+1,stats.hits_after,stats.reflection_m_before+distance,
                        stats.reflection_m_after,stats.external_hits,stats.backtracked_hits)
    end
    push!(events,FaceEvent(-marker,face,distance))
    Contacts(stats.hits_before,stats.hits_after+1,stats.reflection_m_before,
             stats.reflection_m_after+distance,stats.external_hits,
             stats.backtracked_hits+(p<=flow.first_count))
end

# Same coordinate folding as Track.diffusion_axis, with an observer on reflections.
@inline function diffusion_contacts(flow, x, p, first, axis, displacement, t, dx, stats, events, marker)
    while true
        target = x[axis]+displacement
        if 0<=target<=1
            return Base.setindex(x,target,axis),p,first,stats
        end
        side = target>1 ? 1 : -1
        wall = side>0 ? 1. : 0.
        @inbounds next_p = destination(flow,p,axis,side)
        remainder = target-wall
        if next_p==0
            # The reflection changes the unreflected displacement by twice its normal overshoot.
            stats = record_reflection(stats,events,flow,p,axis,side,first,marker,2*abs(remainder)*dx)
            x = Base.setindex(x,wall,axis)
            displacement = -remainder
        else
            x = Base.setindex(x,1-wall,axis)
            displacement = remainder
            p = Int(next_p)
            if isnan(first) && p>flow.first_count
                first = t
            end
        end
    end
end

@inline function trace_contacts(flow, x, p, rng, dx, diffusivity, dt, stop, events, marker)
    first = NaN
    stats = Contacts()
    if diffusivity==0
        _,_,first,exit = advect(flow,x,p,first,0.,stop)
        return first,exit-first,stats
    end
    half = dt/2
    sigma = sqrt(2*diffusivity*dt)/dx
    for step in 1:round(Int,stop/dt)
        t = (step-1)*dt
        x,p,first,exit = advect(flow,x,p,first,t,half)
        isfinite(exit) && return first,exit-first,stats
        for axis in 1:3
            x,p,first,stats = diffusion_contacts(flow,x,p,first,axis,sigma*randn(rng),t+half,dx,stats,events,marker)
        end
        x,p,first,exit = advect(flow,x,p,first,t+half,half)
        isfinite(exit) && return first,exit-first,stats
    end
    first,NaN,stats
end

function sample_contacts(flow, dx, diffusivity, dt, stop, count, seed)
    inlet = findall(i->flow.before[1,i]==0 && flow.low[1,i]>0,1:flow.first_count)
    isempty(inlet) && error("No forward inlet")
    weights = cumsum(flow.low[1,inlet])
    result = Matrix{Float64}(undef,8,count)
    chunks = min(Threads.nthreads(),count)
    events = [FaceEvent[] for _ in 1:chunks]
    Threads.@threads for chunk in 1:chunks
        lo = (chunk-1)*count÷chunks+1
        hi = chunk*count÷chunks
        for i in lo:hi
            rng = Xoshiro(seed+i)
            p = inlet[searchsortedfirst(weights,rand(rng)*weights[end])]
            x = (0.,rand(rng),rand(rng))
            first,second,stats = trace_contacts(flow,x,p,rng,dx,diffusivity,dt,stop,events[chunk],i)
            result[:,i] .= (first,second,stats.hits_before,stats.hits_after,stats.reflection_m_before,
                            stats.reflection_m_after,stats.external_hits,stats.backtracked_hits)
        end
    end
    result,events
end

function main(args)
    length(args)==8 || error("Usage: track_contacts.jl compact output.f64 dx D dt stop count seed")
    source, output = args[1:2]
    dx,d,dt,stop = parse.(Float64,args[3:6])
    count,seed = parse.(Int,args[7:8])
    all(isfinite,(dx,d,dt,stop)) && dx>0 && d>=0 && 0<dt<=stop ||
        error("Invalid geometry, diffusivity or time")
    0<count<=typemax(Int32) || error("Marker count exceeds event identifier range")
    isapprox(stop/dt,round(stop/dt);rtol=1e-12) || error("Stop time must span whole steps")
    Base.ENDIAN_BOM==0x04030201 && sizeof(FaceEvent)==16 || error("Expected a packed little-endian output")
    read_grid(name,::Type{T}) where T = reshape(reinterpret(T,read(joinpath(source,name))),3,:)
    before = read_grid("before.i32",Int32)
    after = read_grid("after.i32",Int32)
    low = read_grid("incoming.f64",Float64) ./ dx^3
    high = read_grid("outgoing.f64",Float64) ./ dx^3
    size(before)==size(after)==size(low)==size(high) && iseven(size(before,2)) ||
        error("Expected two repeated voxel periods")
    grid = Flow(before,after,low,high)
    6*grid.first_count<=typemax(UInt32) || error("Wall face identifier exceeds event range")
    result,events = sample_contacts(grid,dx,d,dt,stop,count,seed)
    open(output,"w") do io
        write(io,result)
    end
    open(output*".faces","w") do io
        for chunk in events
            write(io,reinterpret(UInt8,chunk))
        end
    end
    println("Completed ",sum(isfinite,view(result,2,:))," / ",count,
            " pairs with ",sum(length,events)," solid reflections at D=",d,", dt=",dt)
end

if abspath(PROGRAM_FILE)==@__FILE__
    main(ARGS)
end

end
