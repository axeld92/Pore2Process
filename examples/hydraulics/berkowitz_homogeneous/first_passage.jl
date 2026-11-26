module GlassPassage

using Random
include(joinpath(@__DIR__, "../../adsorption/path_exchange/track.jl"))
using .Track: Flow, crossing, drift, destination

"""Advance through the periodic flow field, retaining the unwrapped axial index.

Each first-passage plane is a voxel face. After crossing the last plane the
marker is absorbed. Axial backtracking upstream of the injection face is free.
"""
@inline function advance(flow, x, p, axial, t, duration, levels, times, next)
    remaining = duration
    while remaining > 0
        @inbounds hits = ntuple(a->crossing(x[a],flow.low[a,p],flow.high[a,p],flow.slope[a,p]),3)
        a = hits[1][1]<=hits[2][1] ? 1 : 2
        a = hits[3][1]<hits[a][1] ? 3 : a
        hit, side = hits[a]
        elapsed = min(hit,remaining)
        @inbounds x = ntuple(b->drift(x[b],flow.low[b,p],flow.high[b,p],flow.slope[b,p],elapsed),3)
        hit > remaining && break
        p = Int(destination(flow,p,a,side))
        p > 0 || error("A nonzero advective flux crossed solid")
        x = Base.setindex(x,side>0 ? 0. : 1.,a)
        t += elapsed
        remaining -= elapsed
        if a == 1
            axial += side
            if axial == levels[next]
                times[next] = t
                next += 1
                next > length(levels) && break
            end
        end
    end
    x,p,axial,next
end

@inline function diffuse(flow, x, p, axial, a, motion, t, levels, times, next)
    side = motion >= 0 ? 1 : -1
    remaining = abs(motion)
    while remaining > 0
        distance = side > 0 ? 1-x[a] : x[a]
        if remaining < distance
            x = Base.setindex(x,clamp(x[a]+side*remaining,0.,1.),a)
            break
        end
        remaining -= distance
        neighbour = Int(destination(flow,p,a,side))
        if neighbour == 0
            x = Base.setindex(x,side > 0 ? 1. : 0.,a)
            side = -side
        else
            p = neighbour
            x = Base.setindex(x,side > 0 ? 0. : 1.,a)
            if a == 1
                axial += side
                if axial == levels[next]
                    times[next] = t
                    next += 1
                    next > length(levels) && break
                end
            end
        end
    end
    x,p,axial,next
end

function arrivals(flow, inlet, levels, h, D, dt, stop, markers, seed)
    all(>(0),diff([0;levels])) || error("Exit planes must increase downstream")
    all(>(0),flow.low[1,inlet]) || error("Injection requires forward-flow faces")
    cdf = cumsum(flow.low[1,inlet])
    data = fill(NaN,length(levels),markers)
    sigma = sqrt(2D*dt)/h
    half = dt/2
    Threads.@threads for marker in 1:markers
        rng = Xoshiro(seed+marker)
        p = Int(inlet[searchsortedfirst(cdf,rand(rng)*last(cdf))])
        x = (0.,rand(rng),rand(rng))
        axial, next = 0,1
        times = view(data,:,marker)
        for step in 1:ceil(Int,stop/dt)
            t = (step-1)*dt
            x,p,axial,next = advance(flow,x,p,axial,t,half,levels,times,next)
            next>length(levels) && break
            for k in 1:3
                a = isodd(step) ? k : 4-k
                x,p,axial,next = diffuse(flow,x,p,axial,a,sigma*randn(rng),t+half,levels,times,next)
                next>length(levels) && break
            end
            next>length(levels) && break
            x,p,axial,next = advance(flow,x,p,axial,t+half,half,levels,times,next)
            next>length(levels) && break
        end
    end
    data
end

function main(args)
    length(args)==9 || error("Use faces output dx D dt stop marker_count seed levels_csv")
    folder, path = args[1:2]
    h,D,dt,stop = parse.(Float64,args[3:6])
    markers, seed = parse.(Int,args[7:8])
    levels = parse.(Int,split(args[9],','))
    load(name,T) = collect(reinterpret(T,read(joinpath(folder,name))))
    grid(name,T) = reshape(load(name,T),3,:)
    flow = Flow(grid("before.i32",Int32),grid("after.i32",Int32),
                grid("low.f64",Float64),grid("high.f64",Float64))
    inlet = load("inlet.i32",Int32)
    data = arrivals(flow,inlet,levels,h,D,dt,stop,markers,seed)
    write(path,data)
    println("Arrivals at each plane ",vec(sum(isfinite,data;dims=2))," / ",markers)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
end
