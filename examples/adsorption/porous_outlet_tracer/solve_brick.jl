include("solve.jl")
include(joinpath(@__DIR__,"../../../julia/adsorption/BrickOpenTracer.jl"))

const Brick = BrickOpenTracer

function verify_grid(grid,legacy,axis,boundary)
    grid.boundary==boundary || error("Boundary rates differ from the archived open faces")
    count(!iszero,grid.forward)==length(legacy.rate) || error("Face count differs")
    @inbounds for e in eachindex(legacy.rate)
        i,j,a = legacy.left[e],legacy.right[e],axis[e]
        grid.forward[a,i]==j && grid.backward[a,j]==i || error("Fine face indexing differs")
        p,s = divrem(i-1,grid.width)
        grid.rate[grid.offset[a,s+1],a,p+1]==legacy.rate[e] || error("Fine face rate differs")
    end
end

function verify(grid,source,periodic,late,dx,d,dt)
    left = load_vector(source,"left.i32",Int32)
    right = load_vector(source,"right.i32",Int32)
    axis = load_vector(source,"axis.u8",UInt8)
    flow = load_vector(source,"balanced.f64",Float64)
    boundary = load_vector(source,"boundary_rate.f64",Float64)
    legacy = Parallel.prepare(left,right,axis,flow,length(boundary),dx,d)
    verify_grid(grid,legacy,axis,boundary)
    states = (("clean",zeros(grid.count),4),
              ("smooth",load_vector(periodic,"varying_initial.f64",Float64),4),
              ("late",load_vector(late,"final.f64",Float64),0))
    println("field,max_concentration_difference_kg_m3,max_relative_outlet_difference")
    for (name,c,pulse_steps) in states
        a = OpenFlow.evolve_open(legacy,c,boundary,dx,dt,pulse_steps,4)
        b = Brick.evolve(grid,c,dt,pulse_steps,4)
        err = maximum(abs,a.concentration-b.concentration)
        scale = max(abs(a.history[end][6]),eps(Float64)*grid.volume)
        mass_error = maximum(abs(x[6]-y[6]) for (x,y) in zip(a.history,b.history))/scale
        println(join((name,err,mass_error),','))
        err==0 && mass_error<=1e-12 || error("Outlet kernel equivalence failed")
    end
end

function driver(args)
    length(args)>=8 || error("Use mode compact_dir output_dir ratio dx D dt steps [open_faces periodic late]")
    mode,source,output = args[1:3]
    ratio = parse(Int,args[4])
    dx,d,dt = parse.(Float64,args[5:7])
    steps = parse(Int,args[8])
    Base.ENDIAN_BOM == 0x04030201 || error("The saved field contract requires little-endian")
    before = reshape(load_vector(source,"before.i32",Int32),3,:)
    after = reshape(load_vector(source,"after.i32",Int32),3,:)
    inflow = reshape(load_vector(source,"incoming.f64",Float64),3,:)
    outflow = reshape(load_vector(source,"outgoing.f64",Float64),3,:)
    grid = Brick.prepare(before,after,inflow,outflow,ratio,dx,d)
    if mode=="verify"
        verify(grid,args[9],args[10],args[11],dx,d,dt)
    elseif mode=="solve"
        println("N384 grid ready, depletion number ",dt*grid.max_depletion,", threads ",Threads.nthreads())
        flush(stdout)
        open(joinpath(output,"history.csv"),"w") do io
            println(io,"time_s,liquid_kg,min_kg_m3,max_kg_m3,injected_kg,recovered_kg,first_moment_kg_s,second_moment_kg_s2,interval_exit_kg")
            flush(io)
            result = Brick.evolve(grid,zeros(grid.count),dt,round(Int,.004/dt),steps;stream=io,progress=true)
            write(joinpath(output,"final.f64"),result.concentration)
        end
    else
        error("Unknown mode")
    end
end

driver(ARGS)
