for module_name in ("BoundedFaceTransport","ThreadedBoundedFaceTransport","OpenBoundedFaceTransport")
    include(joinpath(@__DIR__,"../../../julia/adsorption",module_name*".jl"))
end

const Parallel = ThreadedBoundedFaceTransport
const OpenFlow = OpenBoundedFaceTransport

function load_vector(directory,name,::Type{T}) where T
    path = joinpath(directory,name)
    count,remainder = divrem(filesize(path),sizeof(T))
    remainder==0 || error("Invalid binary length")
    open(io->read!(io,Vector{T}(undef,count)),path,"r")
end

function main(args)
    length(args)==8 || error("Use face_dir output_dir dx diffusivity dt pulse_steps steps threads")
    source,directory = args[1:2]
    dx,d,dt = parse.(Float64,args[3:5])
    pulse_steps,steps,threads = parse.(Int,args[6:8])
    Base.ENDIAN_BOM == 0x04030201 || error("The saved field contract requires a little-endian host")
    Threads.nthreads()==threads || error("Unexpected Julia thread count")
    left = load_vector(source,"left.i32",Int32)
    right = load_vector(source,"right.i32",Int32)
    axis = load_vector(source,"axis.u8",UInt8)
    flow = load_vector(source,"balanced.f64",Float64)
    boundary = load_vector(source,"boundary_rate.f64",Float64)
    grid = Parallel.prepare(left,right,axis,flow,length(boundary),dx,d)
    open(joinpath(directory,"history.csv"),"w") do io
        println(io,"time_s,liquid_kg,min_kg_m3,max_kg_m3,injected_kg,recovered_kg,first_moment_kg_s,second_moment_kg_s2,interval_exit_kg")
        flush(io)
        r = OpenFlow.evolve_open(grid,zeros(length(boundary)),boundary,dx,dt,pulse_steps,steps;progress=true,stream=io)
        write(joinpath(directory,"final.f64"),r.concentration)
    end
end

if abspath(PROGRAM_FILE)==@__FILE__
    main(ARGS)
end
