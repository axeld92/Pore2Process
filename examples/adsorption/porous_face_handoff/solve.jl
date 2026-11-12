include(joinpath(@__DIR__,"../../../julia/adsorption/PorousFaceTransport.jl"))
using .PorousFaceTransport
using LinearAlgebra

function read_vector(directory, name, ::Type{T}) where T
    path = joinpath(directory,name)
    sizeof(T)*div(filesize(path),sizeof(T)) == filesize(path) || error("Invalid binary length")
    data = Vector{T}(undef,filesize(path)÷sizeof(T))
    open(io->read!(io,data),path,"r")
end

function save(directory, prefix, result)
    for (name, c) in zip(("middle","final"),result.snapshots)
        write(joinpath(directory,prefix*"_"*name*".f64"),c)
    end
    open(joinpath(directory,prefix*"_history.csv"),"w") do io
        println(io,"time_s,mean_kg_m3,min_kg_m3,max_kg_m3,integrated_dissipation")
        for row in result.history
            println(io,join(row,','))
        end
    end
end

function main(args)
    length(args)==9 || error("Use directory N dx D Qx Qy Qz t_mid t_final")
    directory = args[1]
    n = parse(Int,args[2])
    dx,d,qx,qy,qz,t_mid,t_final = parse.(Float64,args[3:end])
    Base.ENDIAN_BOM == 0x04030201 || error("The binary contract requires little-endian")
    BLAS.set_num_threads(1)
    left = read_vector(directory,"left.i32",Int32)
    right = read_vector(directory,"right.i32",Int32)
    axis = read_vector(directory,"axis.u8",UInt8)
    raw = read_vector(directory,"unbalanced.f64",Float64)
    initial = read_vector(directory,"varying_initial.f64",Float64)
    count = length(initial)
    balanced = balance_faces(left,right,axis,raw,count,n*[qx,qy,qz])
    write(joinpath(directory,"balanced.f64"),balanced)
    println("Balanced ",length(balanced)," faces for ",count," pore cells");flush(stdout)
    operator = face_operator(left,right,balanced,count,dx,d)
    times = (t_mid,t_final)
    uniform = evolve(operator,fill(.5,count),times,left,right,d/dx^2)
    save(directory,"uniform",uniform)
    varying = evolve(operator,initial,times,left,right,d/dx^2)
    save(directory,"varying",varying)
    if n == 64
        for (name, c) in zip(("middle","final"),exponential_reference(operator,initial,times))
            write(joinpath(directory,"exponential_"*name*".f64"),c)
        end
    elseif n == 96
        save(directory,"half_step",evolve(operator,initial,times,left,right,d/dx^2;time_factor=.5))
    end
end

main(ARGS)
