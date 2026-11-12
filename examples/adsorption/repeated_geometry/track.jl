module TrackRepeats

using Random
include(joinpath(@__DIR__, "../carbon_contact/track_contacts.jl"))
using .TrackContacts: Contacts, FaceEvent, Flow, advect, diffusion_contacts

"""Reuse the two-cell face connections to join copies of the same TRT field."""
function tile(two::Flow, periods)
    periods >= 2 || throw(ArgumentError("At least two periods are required"))
    n = two.first_count
    size(two.before,2) == 2n || error("Expected two archived periods")
    periods*n <= typemax(Int32) || error("Too many tiled voxels")
    if periods == 2
        return [two]
    end
    before = repeat(two.before[:,1:n], 1, periods)
    after = repeat(two.after[:,1:n], 1, periods)
    low = repeat(two.low[:,1:n], 1, periods)
    high = repeat(two.high[:,1:n], 1, periods)
    for k in 1:periods, j in 1:n, a in 1:3
        offset = (k-1)*n
        b = k == 1 ? two.before[a,j] : two.before[a,j+n]
        e = two.after[a,j]
        before[a,j+offset] = b == 0 ? 0 : b+max(k-2,0)*n
        after[a,j+offset] = k == periods && e > n ? 0 : e == 0 ? 0 : e+offset
    end
    slope = high-low
    [Flow(before,after,low,high,slope,k*n) for k in 1:periods-1]
end

mutable struct WallClock{T}
    hits::T
    period::Int
    time::Float64
end

@inline function Base.push!(clock::WallClock, event::FaceEvent)
    phase = min(length(clock.hits), clock.period+(event.marker_phase<0))
    isnan(clock.hits[phase]) && (clock.hits[phase] = clock.time)
    clock
end

@inline function mark_seam!(arrivals, period, seam)
    if period < length(arrivals) && isfinite(seam)
        arrivals[period] = seam
        return period+1
    end
    period
end

@inline function trace!(arrivals,hits,flows,x,p,rng,sigma,dx,dt,steps,marker,period_voxels)
    period = 1
    last = length(arrivals)
    clock = WallClock(hits,1,NaN)
    stats = Contacts()
    half = dt/2
    for step in 1:steps
        t = (step-1)*dt
        flow = flows[min(period,last-1)]
        x,p,seam,exit = advect(flow,x,p,period==last ? Inf : NaN,t,half)
        period = mark_seam!(arrivals,period,seam)
        if isfinite(exit)
            period==last || error("Skipped an internal seam during advection")
            arrivals[last] = exit
            return
        end
        clock.time = t+half
        for axis in 1:3
            displacement = sigma*randn(rng)
            axis==1 && abs(displacement)>=period_voxels && error("Diffusion can skip an axial seam")
            clock.period = period
            flow = flows[min(period,last-1)]
            x,p,seam,stats = diffusion_contacts(flow,x,p,period==last ? Inf : NaN,axis,
                displacement,t+half,dx,stats,clock,marker)
            period = mark_seam!(arrivals,period,seam)
        end
        flow = flows[min(period,last-1)]
        x,p,seam,exit = advect(flow,x,p,period==last ? Inf : NaN,t+half,half)
        period = mark_seam!(arrivals,period,seam)
        if isfinite(exit)
            period==last || error("Skipped an internal seam during advection")
            arrivals[last] = exit
            return
        end
    end
end

function sample(two,dx,D,dt,stop,count,seed,periods,period_voxels)
    flows = tile(two,periods)
    n = two.first_count
    inlet = findall(i->two.before[1,i]==0 && two.low[1,i]>0,1:n)
    isempty(inlet) && error("No forward inlet")
    maximum(abs,view(two.low,1,:))*dt/2 < period_voxels &&
        maximum(abs,view(two.high,1,:))*dt/2 < period_voxels ||
        error("One advection step can skip an axial seam")
    weights = cumsum(two.low[1,inlet])
    data = fill(NaN,2periods,count)
    sigma = sqrt(2D*dt)/dx
    steps = round(Int,stop/dt)
    Threads.@threads for i in 1:count
        rng = Xoshiro(seed+i)
        p = inlet[searchsortedfirst(weights,rand(rng)*weights[end])]
        x = (0.,rand(rng),rand(rng))
        trace!(view(data,1:periods,i),view(data,periods+1:2periods,i),
            flows,x,p,rng,sigma,dx,dt,steps,i,period_voxels)
    end
    data
end

function main(args)
    length(args)==10 || error("Use compact output dx D dt stop count seed periods period_voxels")
    source,output = args[1:2]
    dx,D,dt,stop = parse.(Float64,args[3:6])
    count,seed,periods,period_voxels = parse.(Int,args[7:10])
    Base.ENDIAN_BOM==0x04030201 || error("Expected little-endian storage")
    all(isfinite,(dx,D,dt,stop)) && dx>0 && D>0 && 0<dt<=stop &&
        count>0 && periods>=2 && period_voxels>0 || error("Invalid tracer settings")
    read_grid(name,::Type{T}) where T = reshape(reinterpret(T,read(joinpath(source,name))),3,:)
    before = read_grid("before.i32",Int32)
    after = read_grid("after.i32",Int32)
    low = read_grid("incoming.f64",Float64) ./ dx^3
    high = read_grid("outgoing.f64",Float64) ./ dx^3
    n = size(low,2)÷2
    size(before)==size(after)==size(low)==size(high) && iseven(size(low,2)) ||
        error("Expected a two-period compact grid")
    low[:,1:n]==low[:,n+1:end] && high[:,1:n]==high[:,n+1:end] ||
        error("The archived velocity is not repeated across periods")
    data = sample(Flow(before,after,low,high),dx,D,dt,stop,count,seed,periods,period_voxels)
    open(io->write(io,data),output,"w")
    println("Completed ",Base.count(isfinite,view(data,periods,:))," of ",size(data,2),
        " trajectories through ",periods," cells")
end

if abspath(PROGRAM_FILE)==@__FILE__
    main(ARGS)
end

end
