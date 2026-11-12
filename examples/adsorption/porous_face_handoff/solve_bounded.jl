include(joinpath(@__DIR__,"../../../julia/adsorption/BoundedFaceTransport.jl"))
using .BoundedFaceTransport

function input(directory,name,::Type{T}) where T
    path = joinpath(directory,name)
    n, remainder = divrem(filesize(path),sizeof(T))
    remainder==0 || error("Invalid binary length")
    open(io->read!(io,Vector{T}(undef,n)),path,"r")
end

function solve_case(args)
    length(args)==8 || error("Use source_dir output_dir dx D t_mid t_final uniform_or_varying time_factor")
    source, directory = args[1:2]
    dx,d,t1,t2 = parse.(Float64,args[3:6])
    left = input(source,"left.i32",Int32)
    right = input(source,"right.i32",Int32)
    axis = input(source,"axis.u8",UInt8)
    flow = input(source,"balanced.f64",Float64)
    initial = input(source,"varying_initial.f64",Float64)
    kind = args[7]
    kind in ("uniform","varying","half_step") || error("Unknown tracer")
    if kind=="uniform"
        fill!(initial,.5)
    end
    result = evolve_bounded(left,right,axis,flow,initial,dx,d,(t1,t2);time_factor=parse(Float64,args[8]))
    for (name, c) in zip(("middle","final"),result.snapshots)
        write(joinpath(directory,kind*"_"*name*".f64"),c)
    end
    open(joinpath(directory,kind*"_history.csv"),"w") do io
        println(io,"time_s,mean_kg_m3,min_kg_m3,max_kg_m3,integrated_dissipation")
        for row in result.history
            println(io,join(row,','))
        end
    end
end

solve_case(ARGS)
