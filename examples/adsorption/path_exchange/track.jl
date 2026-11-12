module Track

using Random

struct Flow
    before::Matrix{Int32}
    after::Matrix{Int32}
    low::Matrix{Float64}
    high::Matrix{Float64}
    slope::Matrix{Float64}
    first_count::Int
end

Flow(before, after, low, high) = Flow(before, after, low, high, high-low, size(low,2)÷2)

# x is in voxel widths. The ODE inside one voxel is x' = low + slope*x.
@inline function crossing(x, low, high, slope)
    speed = muladd(x,high,(1-x)*low)
    if speed>0 && high>0
        distance = 1-x
        face_speed = high
    elseif speed<0 && low<0
        distance = -x
        face_speed = low
    else
        return Inf,0
    end
    if slope==0
        dt = distance/speed
    else
        change = slope*distance/speed
        log_ratio = abs(change)<=.5 ? log1p(change) : log(abs(face_speed))-log(abs(speed))
        dt = log_ratio/slope
    end
    max(dt,0.),speed>0 ? 1 : -1
end

@inline function drift(x, low, high, slope, dt)
    factor = slope==0 ? dt : expm1(slope*dt)/slope
    clamp(muladd(muladd(x,high,(1-x)*low),factor,x),0.,1.)
end

@inline function destination(flow, p, a, side)
    side>0 ? flow.after[a,p] : flow.before[a,p]
end

@inline function advect(flow, x, p, first, t, dt)
    stop = t+dt
    while t<stop
        @inbounds begin
            times = ntuple(a->crossing(x[a],flow.low[a,p],flow.high[a,p],flow.slope[a,p]),3)
            a = times[1][1]<=times[2][1] ? 1 : 2
            a = times[3][1]<times[a][1] ? 3 : a
            hit,side = times[a]
            remaining = stop-t
            advance = min(hit,remaining)
            next_x = ntuple(i->drift(x[i],flow.low[i,p],flow.high[i,p],flow.slope[i,p],advance),3)
            if hit>remaining
                return next_x,p,first,NaN
            end
            t += hit
            next_p = destination(flow,p,a,side)
            if next_p==0
                (a==1 && side>0 && p>flow.first_count) || error("Invalid advective exit")
                return next_x,p,first,t
            end
            x = Base.setindex(next_x,side>0 ? 0. : 1.,a)
            p = Int(next_p)
            if isnan(first) && p>flow.first_count
                first = t
            end
        end
    end
    x,p,first,NaN
end

@inline function diffusion_axis(flow, x, p, first, a, displacement, t)
    while true
        target = x[a]+displacement
        if 0<=target<=1
            return Base.setindex(x,target,a),p,first
        end
        side = target>1 ? 1 : -1
        wall = side>0 ? 1. : 0.
        @inbounds next_p = destination(flow,p,a,side)
        remainder = target-wall
        if next_p==0
            x = Base.setindex(x,wall,a)
            displacement = -remainder
        else
            x = Base.setindex(x,1-wall,a)
            displacement = remainder
            p = Int(next_p)
            if isnan(first) && p>flow.first_count
                first = t
            end
        end
    end
end

@inline function trace(flow, x, p, rng, dx, diffusivity, dt, stop)
    first = NaN
    if diffusivity==0
        _,_,first,exit = advect(flow,x,p,first,0.,stop)
        return first,exit-first
    end
    half = dt/2
    sigma = sqrt(2*diffusivity*dt)/dx
    for step in 1:round(Int,stop/dt)
        t = (step-1)*dt
        x,p,first,exit = advect(flow,x,p,first,t,half)
        isfinite(exit) && return first,exit-first
        for a in 1:3
            x,p,first = diffusion_axis(flow,x,p,first,a,sigma*randn(rng),t+half)
        end
        x,p,first,exit = advect(flow,x,p,first,t+half,half)
        isfinite(exit) && return first,exit-first
    end
    first,NaN
end

function sample(flow, dx, diffusivity, dt, stop, count, seed)
    inlet = findall(i->flow.before[1,i]==0 && flow.low[1,i]>0,1:flow.first_count)
    isempty(inlet) && error("No forward inlet")
    weights = cumsum(flow.low[1,inlet])
    times = Matrix{Float64}(undef,2,count)
    Threads.@threads for i in 1:count
        rng = Xoshiro(seed+i)
        p = inlet[searchsortedfirst(weights,rand(rng)*weights[end])]
        x = (0.,rand(rng),rand(rng))
        first,second = trace(flow,x,p,rng,dx,diffusivity,dt,stop)
        times[1,i] = first
        times[2,i] = second
    end
    times
end

function main(args)
    source, output = args[1:2]
    dx,d,dt,stop = parse.(Float64,args[3:6])
    count,seed = parse.(Int,args[7:8])
    read_grid(name,::Type{T}) where T = reshape(reinterpret(T,read(joinpath(source,name))),3,:)
    before = read_grid("before.i32",Int32)
    after = read_grid("after.i32",Int32)
    low = read_grid("incoming.f64",Float64) ./ dx^3
    high = read_grid("outgoing.f64",Float64) ./ dx^3
    grid = Flow(before,after,low,high)
    times = sample(grid,dx,d,dt,stop,count,seed)
    open(output,"w") do io
        write(io,times)
    end
    println("Completed ",sum(isfinite,view(times,2,:))," / ",count," pairs at D=",d,", dt=",dt)
end

if abspath(PROGRAM_FILE)==@__FILE__
    main(ARGS)
end

end
