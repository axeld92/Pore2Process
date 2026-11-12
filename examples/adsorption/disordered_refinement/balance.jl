module RefinedFaceProjection

using LinearAlgebra, SparseArrays
include(joinpath(@__DIR__,"../../../julia/adsorption/PorousFaceTransport.jl"))
using .PorousFaceTransport: divergence, graph_laplacian

struct AggregationLevel{C,M}
    matrix::M
    invdiag::Vector{Float64}
    groups::Vector{Int32}
    coarse::C
    coarse_count::Int
end

struct CoarseSolve{F}
    factor::F
end

function aggregate(coords)
    shape = ntuple(d->maximum(view(coords,d,:))÷2+1,3)
    ids = zeros(Int32,shape)
    groups = Vector{Int32}(undef,size(coords,2))
    coarse = NTuple{3,Int32}[]
    @inbounds for j in eachindex(groups)
        ijk = ntuple(d->coords[d,j]÷2+1,3)
        id = ids[ijk...]
        if id==0
            push!(coarse,ijk .- 1)
            id = ids[ijk...] = length(coarse)
        end
        groups[j] = id
    end
    centres = reshape(reinterpret(Int32,coarse),3,:)
    groups, centres
end

function hierarchy(A,coords)
    n = size(A,1)
    if n<=1200
        return CoarseSolve(cholesky(Symmetric(A)))
    end
    groups, centres = aggregate(coords)
    coarse_count = size(centres,2)
    P = sparse(Int32.(1:n),groups,ones(n),n,coarse_count)
    Ac = sparse(transpose(P)*(A*P))
    AggregationLevel(A,1 ./ Vector(diag(A)),groups,hierarchy(Ac,centres),coarse_count)
end

struct CycleWork{W}
    product::Vector{Float64}
    coarse_rhs::Vector{Float64}
    coarse_solution::Vector{Float64}
    coarse_work::W
end

work(::CoarseSolve) = nothing
work(level::AggregationLevel) = CycleWork(zeros(size(level.matrix,1)),
    zeros(level.coarse_count),zeros(level.coarse_count),work(level.coarse))

function precondition!(z,level::CoarseSolve,r,::Nothing)
    z .= level.factor\r
end

function precondition!(z,level::AggregationLevel,r,w)
    A, D = level.matrix, level.invdiag
    omega = 2/3
    @. z = omega*D*r
    mul!(w.product,A,z)
    @. z += omega*D*(r-w.product)
    mul!(w.product,A,z)
    fill!(w.coarse_rhs,0.)
    @inbounds for j in eachindex(z)
        w.coarse_rhs[level.groups[j]] += r[j]-w.product[j]
    end
    precondition!(w.coarse_solution,level.coarse,w.coarse_rhs,w.coarse_work)
    @inbounds for j in eachindex(z)
        z[j] += w.coarse_solution[level.groups[j]]
    end
    for _ in 1:2
        mul!(w.product,A,z)
        @. z += omega*D*(r-w.product)
    end
end

function potential(A,levels,rhs)
    x = zeros(length(rhs))
    r = copy(rhs)
    z,p,product = (similar(rhs) for _ in 1:3)
    tolerance = 4e-14*max(norm(rhs),1.)
    norm(r)<=tolerance && return x,0
    w = work(levels)
    precondition!(z,levels,r,w)
    copyto!(p,z)
    rz = dot(r,z)
    for iteration in 1:600
        mul!(product,A,p)
        curvature = dot(p,product)
        curvature>0 && rz>0 || error("Nonpositive preconditioned graph curvature")
        step = rz/curvature
        @. x += step*p
        @. r -= step*product
        norm(r)<=tolerance && return x,iteration
        precondition!(z,levels,r,w)
        next = dot(r,z)
        @. p = z+(next/rz)*p
        rz = next
    end
    error("Preconditioned graph projection did not converge")
end

function project(A,levels,left,right,rates,name)
    phi, iterations = potential(A,levels,divergence(left,right,rates,size(A,1)))
    @inbounds for e in eachindex(rates)
        rates[e] += phi[right[e]]-phi[left[e]]
    end
    println(name," completed in ",iterations," iterations")
    flush(stdout)
    rates
end

"""Same constrained least-squares projection as PorousFaceTransport."""
function balance(left,right,axis,raw,coords,directional_sums)
    cells = size(coords,2)
    scale = maximum(abs,raw)
    scale>0 || error("Expected a nonzero carrier field")
    matrix = graph_laplacian(left,right,cells)
    levels = hierarchy(matrix,coords)
    println("Built graph projection for ",cells," cells")
    flush(stdout)
    # Reuse the same graph and coarse factor for all four projections.
    balanced = project(matrix,levels,left,right,raw/scale,"raw")
    basis = [project(matrix,levels,left,right,Float64.(axis.==d),"axis $d") for d in 1:3]
    totals(flow) = [sum(flow[e] for e in eachindex(flow) if axis[e]==d) for d in 1:3]
    coupling = hcat(map(totals,basis)...)
    coefficients = coupling\(directional_sums/scale-totals(balanced))
    for d in 1:3
        @. balanced += basis[d]*coefficients[d]
    end
    balanced .*= scale
    net = divergence(left,right,balanced,cells)
    maximum(abs,net)<=1e-10*maximum(abs,balanced) || error("Projected flow does not close locally")
    balanced
end

function main(args)
    length(args)==5 || error("Use directory liquid_cell_count sum_Qx sum_Qy sum_Qz")
    dir = args[1]
    cells = parse(Int,args[2])
    totals = parse.(Float64,args[3:end])
    BLAS.set_num_threads(1)
    load(name,T) = reinterpret(T,read(joinpath(dir,name)))
    left = load("left.i32",Int32)
    right = load("right.i32",Int32)
    axis = load("axis.u8",UInt8)
    raw = load("unbalanced.f64",Float64)
    coords = reshape(load("coords.i32",Int32),3,cells)
    q = balance(left,right,axis,raw,coords,totals)
    write(joinpath(dir,"balanced.f64"),q)
end

if abspath(PROGRAM_FILE)==@__FILE__
    main(ARGS)
end

end
