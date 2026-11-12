include(joinpath(@__DIR__,"../../../julia/adsorption/BoundedFaceTransport.jl"))
include(joinpath(@__DIR__,"../../../julia/adsorption/ThreadedBoundedFaceTransport.jl"))

const Serial = BoundedFaceTransport
const Parallel = ThreadedBoundedFaceTransport

function load_vector(directory,name,::Type{T}) where T
    path = joinpath(directory,name)
    count,remainder = divrem(filesize(path),sizeof(T))
    remainder==0 || error("Invalid binary length")
    open(io->read!(io,Vector{T}(undef,count)),path,"r")
end

function save(directory,result)
    for (name,c) in zip(("middle","final"),result.snapshots)
        write(joinpath(directory,"varying_"*name*".f64"),c)
    end
    open(joinpath(directory,"varying_history.csv"),"w") do io
        println(io,"time_s,mean_kg_m3,min_kg_m3,max_kg_m3,integrated_dissipation")
        for row in result.history
            println(io,join(row,','))
        end
    end
end

function verify(grid,initial,archived,dt)
    fields = (initial,load_vector(archived,"varying_middle.f64",Float64),
              load_vector(archived,"varying_final.f64",Float64))
    for (name,c) in zip(("initial","middle","final"),fields)
        serial = copy(c)
        work = Serial.workspace(c,length(grid.rate))
        a,b = similar(c),similar(c)
        dissipation = 0.
        for _ in 1:4
            d1 = Serial.limited_step!(a,serial,grid,work,dt,grid.d)
            d2 = Serial.limited_step!(b,a,grid,work,dt,grid.d)
            @. b = serial+(b-serial)/4
            d3 = Serial.limited_step!(a,b,grid,work,dt,grid.d)
            @. serial += 2*(a-serial)/3
            dissipation += dt*(d1+d2+4d3)/6
        end
        # Set the interval using the exact step size, so time is not a variable.
        threaded = Parallel.evolve_threaded(grid,c,(4dt,),(4,))
        difference = maximum(abs,only(threaded.snapshots)-serial)
        relative = abs(threaded.history[end][5]/dissipation-1)
        println(join((name,difference,relative),','))
        difference==0 && relative<=1e-10 || error("Threaded equivalence failed")
    end
end

function main(args)
    length(args)==10 || error("Use mode face_dir output_or_archived_dir dx D t1 t2 steps1 steps2 threads")
    mode,source,output = args[1:3]
    dx,d,t1,t2 = parse.(Float64,args[4:7])
    n1,n2,threads = parse.(Int,args[8:10])
    Base.ENDIAN_BOM == 0x04030201 || error("The saved field contract requires a little-endian host")
    Threads.nthreads()==threads || error("Unexpected Julia thread count")
    left = load_vector(source,"left.i32",Int32)
    right = load_vector(source,"right.i32",Int32)
    axis = load_vector(source,"axis.u8",UInt8)
    flow = load_vector(source,"balanced.f64",Float64)
    initial = load_vector(source,"varying_initial.f64",Float64)
    grid = Parallel.prepare(left,right,axis,flow,length(initial),dx,d)
    if mode=="verify"
        println("field,max_concentration_difference_kg_m3,relative_dissipation_difference")
        verify(grid,initial,output,t1/n1)
    elseif mode=="solve"
        println("Tracer grid ready, cap ",grid.cap," s, threads ",threads)
        flush(stdout)
        save(output,Parallel.evolve_threaded(grid,initial,(t1,t2),(n1,n2);progress=true))
    else
        error("Unknown mode")
    end
end

main(ARGS)
