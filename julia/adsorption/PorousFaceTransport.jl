module PorousFaceTransport

using LinearAlgebra
using SparseArrays

export balance_faces, face_operator, evolve, exponential_reference

"""Return net outward volume flow at each cell from shared oriented faces."""
function divergence(left, right, flow, count)
    net = zeros(count)
    @inbounds for e in eachindex(flow)
        net[left[e]] += flow[e]
        net[right[e]] -= flow[e]
    end
    net
end

function graph_laplacian(left, right, count)
    weight = ones(length(left))
    matrix = sparse(vcat(left,right,left,right), vcat(left,right,right,left),
                    vcat(weight,weight,-weight,-weight), count, count)
    any(iszero, diag(matrix)) && error("Isolated pore cell")
    matrix[1,1] += 1
    matrix
end

# The pinned Laplacian has the constant null mode removed. Jacobi CG keeps
# the memory cost proportional to the pore graph.
function potential(matrix, invdiag, rhs)
    x = zeros(length(rhs))
    r = copy(rhs)
    z = r .* invdiag
    p = copy(z)
    product = similar(r)
    rz = dot(r,z)
    tolerance = 4e-14*max(norm(rhs), 1.0)
    norm(r) <= tolerance && return x
    for _ in 1:6000
        mul!(product,matrix,p)
        curvature = dot(p,product)
        curvature > 0 || error("Nonpositive graph Laplacian curvature")
        step = rz/curvature
        @. x += step*p
        @. r -= step*product
        norm(r) <= tolerance && return x
        @. z = r*invdiag
        next = dot(r,z)
        @. p = z+(next/rz)*p
        rz = next
    end
    error("Graph projection did not converge")
end

function project(matrix, invdiag, left, right, flow)
    phi = potential(matrix,invdiag,divergence(left,right,flow,size(matrix,1)))
    [flow[e]-phi[left[e]]+phi[right[e]] for e in eachindex(flow)]
end

"""Least-squares face correction with local balance and all directional sums fixed."""
function balance_faces(left, right, axis, raw, cells, directional_sums)
    scale = maximum(abs,raw)
    scale > 0 || throw(ArgumentError("A nonzero saved carrier field is required"))
    matrix = graph_laplacian(left,right,cells)
    invdiag = 1 ./ diag(matrix)
    normalized = raw/scale
    balanced = project(matrix,invdiag,left,right,normalized)
    # Project the directional indicator columns onto the same divergence-free
    # subspace. They span the least-norm adjustment of the net periodic flows.
    basis = hcat([project(matrix,invdiag,left,right,Float64.(axis .== d)) for d in 1:3]...)
    totals(flow) = [sum(flow[e] for e in eachindex(flow) if axis[e]==d) for d in 1:3]
    coupling = hcat([totals(view(basis,:,d)) for d in 1:3]...)
    balanced += basis*(coupling\(directional_sums/scale-totals(balanced)))
    balanced *= scale
    net = divergence(left,right,balanced,cells)
    maximum(abs,net) <= 1e-10*maximum(abs,balanced) || error("Projected flow does not close locally")
    balanced
end

"""Concentration-rate matrix for one centered shared face flux per edge."""
function face_operator(left, right, flow, cells, dx, diffusivity)
    diff_rate = diffusivity/dx^2
    adv_rate = flow/(2dx^3)
    sparse(vcat(left,left,right,right),vcat(left,right,left,right),
           vcat(-diff_rate .- adv_rate,diff_rate .- adv_rate,
                 diff_rate .+ adv_rate,-diff_rate .+ adv_rate),cells,cells)
end

function dissipation(c, left, right, diffusion_rate)
    total = 0.0
    @inbounds @simd for e in eachindex(left)
        total += abs2(c[left[e]]-c[right[e]])
    end
    diffusion_rate*total/length(c)
end

operator_bound(matrix) = maximum(vec(sum(abs,matrix;dims=2)))

"""RK4 time evolution with a matched quadrature of physical variance dissipation."""
function evolve(matrix, initial, times, left, right, diffusion_rate; time_factor=1)
    c = copy(initial)
    k1, k2, k3, k4, stage = (similar(c) for _ in 1:5)
    snapshots = Vector{Float64}[]
    history = NTuple{5,Float64}[]
    dissipated = 0.0
    function record(time)
        lo, hi = extrema(c)
        all(isfinite, (lo,hi,dissipated)) || error("Non-finite tracer state")
        push!(history,(time,sum(c)/length(c),lo,hi,dissipated))
    end
    record(0.)
    previous = 0.0
    dt_max = time_factor*.75/operator_bound(matrix)
    for t in times
        steps = ceil(Int,(t-previous)/dt_max)
        dt = (t-previous)/steps
        for step in 1:steps
            d1 = dissipation(c,left,right,diffusion_rate)
            mul!(k1,matrix,c)
            @. stage = c+dt*k1/2
            d2 = dissipation(stage,left,right,diffusion_rate)
            mul!(k2,matrix,stage)
            @. stage = c+dt*k2/2
            d3 = dissipation(stage,left,right,diffusion_rate)
            mul!(k3,matrix,stage)
            @. stage = c+dt*k3
            d4 = dissipation(stage,left,right,diffusion_rate)
            mul!(k4,matrix,stage)
            @. c += dt*(k1+2k2+2k3+k4)/6
            dissipated += dt*(d1+2d2+2d3+d4)/6
            record(step==steps ? t : previous+step*dt)
        end
        push!(snapshots,copy(c))
        previous = t
    end
    (;snapshots,history)
end

"""Independent 18-term Taylor exponential action, with each argument norm <= 0.75."""
function exponential_reference(matrix, initial, times)
    c = copy(initial)
    term, scratch = similar(c), similar(c)
    out = Vector{Float64}[]
    previous = 0.0
    dt_max = .75/operator_bound(matrix)
    for t in times
        steps = ceil(Int,(t-previous)/dt_max)
        dt = (t-previous)/steps
        for _ in 1:steps
            copyto!(term,c)
            for k in 1:18
                mul!(scratch,matrix,term)
                @. term = (dt/k)*scratch
                c .+= term
            end
        end
        push!(out,copy(c))
        previous = t
    end
    out
end

end
