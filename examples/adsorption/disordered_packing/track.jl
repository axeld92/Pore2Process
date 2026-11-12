include(joinpath(@__DIR__,"../repeated_geometry/track.jl"))
using .TrackRepeats: Flow, trace!
using Random

function sample(directory,dx,D,dt,stop,count,seed,period_voxels,counts)
    read_grid(name,T) = reshape(reinterpret(T,read(joinpath(directory,name))),3,:)
    before = read_grid("before.i32",Int32)
    after = read_grid("after.i32",Int32)
    low = read_grid("incoming.f64",Float64)./dx^3
    high = read_grid("outgoing.f64",Float64)./dx^3
    size(before)==size(after)==size(low)==size(high) && size(low,2)==last(counts) &&
        all(>(0),diff([0;counts])) || error("Invalid eight-cell compact field")
    slope = high-low
    flows = [Flow(before,after,low,high,slope,k) for k in counts[1:end-1]]
    maximum(abs,view(low,1,:))*dt/2 < period_voxels &&
        maximum(abs,view(high,1,:))*dt/2 < period_voxels ||
        error("One advection step can skip an axial seam")
    inlet = findall(i->before[1,i]==0 && low[1,i]>0,1:counts[1])
    isempty(inlet) && error("No forward inlet")
    weights = cumsum(low[1,inlet])
    periods = length(counts)
    data = fill(NaN,2periods,count)
    steps = round(Int,stop/dt)
    sigma = sqrt(2D*dt)/dx
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
    length(args)==10 || error("Use compact output dx D dt stop count seed period_voxels end_counts_csv")
    dx,D,dt,stop = parse.(Float64,args[3:6])
    count,seed,width = parse.(Int,args[7:9])
    ends = parse.(Int,split(args[10],','))
    all(isfinite,(dx,D,dt,stop)) && dx>0 && D>0 && 0<dt<=stop &&
        count>0 && width>0 && length(ends)==8 || error("Invalid tracer settings")
    data = sample(args[1],dx,D,dt,stop,count,seed,width,ends)
    open(io->write(io,data),args[2],"w")
    println("Inert outlet recovery ",Base.count(isfinite,view(data,8,:))," of ",size(data,2))
end

if abspath(PROGRAM_FILE)==@__FILE__
    main(ARGS)
end
