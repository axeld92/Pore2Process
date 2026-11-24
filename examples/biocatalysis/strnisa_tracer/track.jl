# Same exact-within-voxel advection and reflected molecular moves
# as the previously qualified porous tracer. The x seam is open.
include(joinpath(@__DIR__,"../../adsorption/path_exchange/track.jl"))
using .Track: Flow, sample

function passage(before,after,low,high,first_count,dx,d,dt,stop,count,seed)
    size(before)==size(after)==size(low)==size(high) ||
        error("Inconsistent compact flow arrays")
    0<first_count<size(low,2) || error("Invalid interior partition")
    d>=0 && dx>0 && dt>0 && stop>dt && count>0 ||
        throw(ArgumentError("Invalid tracer input"))
    12sqrt(2d*dt)/dx<1 || error("Brownian increment outside previously qualified bound")
    grid=Flow(before,after,low/dx^3,high/dx^3,(high-low)/dx^3,first_count)
    pair=sample(grid,dx,d,dt,stop,count,seed)
    vec(pair[1,:]+pair[2,:])
end

function main(args)
    length(args)==9 || error("Use compact output midpoint_cell_count dx D dt stop count seed")
    folder,output=args[1:2]
    first_count=parse(Int,args[3])
    dx,d,dt,stop=parse.(Float64,args[4:7])
    count,seed=parse.(Int,args[8:9])
    Base.ENDIAN_BOM==0x04030201 || error("Expected little-endian binary arrays")
    load(file,T)=reshape(reinterpret(T,read(joinpath(folder,file))),3,:)
    before=load("before.i32",Int32); after=load("after.i32",Int32)
    low=load("incoming.f64",Float64); high=load("outgoing.f64",Float64)
    times=passage(before,after,low,high,first_count,dx,d,dt,stop,count,seed)
    write(output,times)
    println("Recovered ",sum(isfinite,times)," / ",count," tracer markers")
end
if abspath(PROGRAM_FILE)==@__FILE__
    main(ARGS)
end
