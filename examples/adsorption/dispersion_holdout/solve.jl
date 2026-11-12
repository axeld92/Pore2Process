include(joinpath(@__DIR__,"../porous_outlet_tracer/solve.jl"))
include(joinpath(@__DIR__,"../../../julia/adsorption/BrickOpenTracer.jl"))
using .BrickOpenTracer: prepare, evolve

struct ProfileRejected <: Exception end

mutable struct HoldoutMonitor
    file::IO
    predicted::Vector{Float64}
    dose::Float64
    sample::Int
    sum_squared::Float64
end

Base.flush(m::HoldoutMonitor) = flush(m.file)

function Base.println(m::HoldoutMonitor,row::String)
    println(m.file,row)
    m.sample += 1
    recovered = parse(Float64,split(row,',')[6])/m.dose
    m.sum_squared += (recovered-m.predicted[m.sample])^2
    bound = sqrt(m.sum_squared/length(m.predicted))
    if recovered>=.1 && bound>.03
        flush(m)
        println("Irreversible profile rejection at sample ",m.sample,", RMSE lower bound ",bound)
        flush(stdout)
        throw(ProfileRejected())
    end
end

function driver(args)
    source,output,model = args[1:3]
    dx,d,dt,dose = parse.(Float64,args[4:7])
    steps = parse(Int,args[8])
    before = reshape(load_vector(source,"before.i32",Int32),3,:)
    after = reshape(load_vector(source,"after.i32",Int32),3,:)
    inflow = reshape(load_vector(source,"incoming.f64",Float64),3,:)
    outflow = reshape(load_vector(source,"outgoing.f64",Float64),3,:)
    grid = prepare(before,after,inflow,outflow,3,dx,d)
    println("Two-cell grid ready, cells ",grid.count,", depletion number ",dt*grid.max_depletion)
    flush(stdout)
    predicted = open(model) do io
        readline(io)
        [parse(Float64,split(row,',')[6])/dose for row in eachline(io)]
    end
    length(predicted)==steps+1 || error("The frozen model must span the planned window")
    status = open(joinpath(output,"history.csv"),"w") do io
        println(io,"time_s,liquid_kg,min_kg_m3,max_kg_m3,injected_kg,recovered_kg,first_moment_kg_s,second_moment_kg_s2,interval_exit_kg")
        flush(io)
        monitor = HoldoutMonitor(io,predicted,dose,0,0.)
        try
            result = evolve(grid,zeros(grid.count),dt,round(Int,.004/dt),steps;stream=monitor,progress=true)
            write(joinpath(output,"final.f64"),result.concentration)
            "complete"
        catch error
            error isa ProfileRejected || rethrow()
            "profile_rejected"
        end
    end
    open(joinpath(output,"status.txt"),"w") do io
        println(io,status)
    end
end

driver(ARGS)
