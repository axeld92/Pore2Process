# Reuse the project pore-graph projection, with only one periodic axis.
include(joinpath(@__DIR__, "../../../../julia/adsorption/PorousFaceTransport.jl"))
using .PorousFaceTransport: graph_laplacian, project, divergence
using LinearAlgebra

function main(args)
    length(args) == 3 || error("Use directory liquid_cell_count sum_axial_flow_m3_s")
    dir, count_arg, flow_arg = args
    count = parse(Int, count_arg)
    total = parse(Float64, flow_arg)
    BLAS.set_num_threads(1)
    load(name, T) = reinterpret(T, read(joinpath(dir, name)))
    left = load("left.i32", Int32)
    right = load("right.i32", Int32)
    axial = load("axis.u8", UInt8) .== 1
    raw = load("unbalanced.f64", Float64)
    scale = maximum(abs, raw)
    scale > 0 && total > 0 && any(axial) || error("Expected nonzero axial flow")
    matrix = graph_laplacian(left, right, count)
    inverse_diagonal = 1 ./ diag(matrix)
    balance(flow) = project(matrix, inverse_diagonal, left, right, flow)
    q = balance(raw / scale)
    constraint = balance(Float64.(axial))
    sum_constraint = sum(view(constraint, axial))
    sum_constraint > 0 || error("Pore topology does not conduct axially")
    q += constraint * ((total / scale - sum(view(q, axial))) / sum_constraint)
    q *= scale
    maximum(abs, divergence(left, right, q, count)) <= 1e-10 * maximum(abs, q) ||
        error("Projected flow does not close locally")
    write(joinpath(dir, "balanced.f64"), q)
end

main(ARGS)
