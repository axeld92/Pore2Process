include(joinpath(@__DIR__,"../../../julia/adsorption/PorousFaceTransport.jl"))
using .PorousFaceTransport: balance_faces
using LinearAlgebra

function main(args)
    length(args)==5 || error("Use directory liquid_cell_count sum_Qx sum_Qy sum_Qz")
    dir = args[1]
    cells = parse(Int,args[2])
    directional_sums = parse.(Float64,args[3:end])
    BLAS.set_num_threads(1)
    load(name,T) = reinterpret(T,read(joinpath(dir,name)))
    left = load("left.i32",Int32)
    right = load("right.i32",Int32)
    axis = load("axis.u8",UInt8)
    raw = load("unbalanced.f64",Float64)
    q = balance_faces(left,right,axis,raw,cells,directional_sums)
    write(joinpath(dir,"balanced.f64"),q)
end

main(ARGS)
