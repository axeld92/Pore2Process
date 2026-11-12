module BoundedFaceTransport

export evolve_bounded

"""Signed neighbor beyond either side of a face, or zero across a wall."""
function extended_stencil(left, right, axis, count)
    before, after = zeros(Int32,count,3), zeros(Int32,count,3)
    for e in eachindex(left)
        before[right[e],axis[e]] = left[e]
        after[left[e],axis[e]] = right[e]
    end
    [before[left[e],axis[e]] for e in eachindex(left)],
    [after[right[e],axis[e]] for e in eachindex(right)]
end

function stability_cap(left, right, rate, diff_rate, count)
    outgoing = zeros(count)
    for e in eachindex(rate)
        outgoing[left[e]] += diff_rate+max(rate[e],0.)
        outgoing[right[e]] += diff_rate+max(-rate[e],0.)
    end
    .8/maximum(outgoing)
end

function workspace(c, faces)
    (; lower=similar(c), upper=similar(c), positive=similar(c), negative=similar(c),
       up=similar(c), down=similar(c), antidiffusion=zeros(faces))
end

# A forward Euler update with paired conservative flux corrections.
# Bounds are imposed on each correction factor, not on concentrations.
function limited_step!(out, c, grid, work, dt, d)
    copyto!(out,c)
    copyto!(work.lower,c)
    copyto!(work.upper,c)
    fill!(work.positive,0.)
    fill!(work.negative,0.)
    energy = 0.0
    @inbounds for e in eachindex(grid.rate)
        i, j, q = grid.left[e], grid.right[e], grid.rate[e]
        a, b = c[i], c[j]
        work.lower[i] = min(work.lower[i],b)
        work.lower[j] = min(work.lower[j],a)
        work.upper[i] = max(work.upper[i],b)
        work.upper[j] = max(work.upper[j],a)
        diff = a-b
        energy += diff*diff
        donor = ifelse(q>=0,a,b)
        low = dt*(q*donor+d*diff)
        out[i] -= low
        out[j] += low
        high = (a+b)/2
        up, down = grid.before[e], grid.after[e]
        if up!=0 && down!=0
            high += (a+b-c[up]-c[down])/12
        end
        delta = dt*q*(high-donor)
        work.antidiffusion[e] = delta
        gain, loss = max(delta,0.), max(-delta,0.)
        work.positive[i] += loss
        work.positive[j] += gain
        work.negative[i] += gain
        work.negative[j] += loss
    end
    @inbounds for i in eachindex(c)
        work.up[i] = work.positive[i]>0 ? clamp((work.upper[i]-out[i])/work.positive[i],0.,1.) : 1.
        work.down[i] = work.negative[i]>0 ? clamp((out[i]-work.lower[i])/work.negative[i],0.,1.) : 1.
    end
    @inbounds for e in eachindex(grid.rate)
        i, j, delta = grid.left[e], grid.right[e], work.antidiffusion[e]
        factor = delta>=0 ? min(work.down[i],work.up[j]) : min(work.up[i],work.down[j])
        correction = factor*delta
        out[i] -= correction
        out[j] += correction
    end
    d*energy/length(c)
end

"""Conservative flux-corrected tracer transport with SSP RK3 integration."""
function evolve_bounded(left, right, axis, flow, initial, dx, diffusivity, times; time_factor=1)
    c = copy(initial)
    rate = flow/dx^3
    d = diffusivity/dx^2
    before, after = extended_stencil(left,right,axis,length(c))
    grid = (;left,right,rate,before,after)
    dt_max = time_factor*stability_cap(left,right,rate,d,length(c))
    work = workspace(c,length(flow))
    a, b = similar(c), similar(c)
    snapshots = Vector{Float64}[]
    history = NTuple{5,Float64}[]
    dissipated = 0.
    function record(t)
        lo, hi = extrema(c)
        all(isfinite,(lo,hi,dissipated)) || error("Non-finite tracer state")
        push!(history,(t,sum(c)/length(c),lo,hi,dissipated))
    end
    record(0.)
    previous = 0.
    for t in times
        steps = ceil(Int,(t-previous)/dt_max)
        dt = (t-previous)/steps
        for step in 1:steps
            d1 = limited_step!(a,c,grid,work,dt,d)
            d2 = limited_step!(b,a,grid,work,dt,d)
            @. b = c+(b-c)/4
            d3 = limited_step!(a,b,grid,work,dt,d)
            @. c += 2*(a-c)/3
            dissipated += dt*(d1+d2+4d3)/6
            record(step==steps ? t : previous+step*dt)
        end
        push!(snapshots,copy(c))
        previous = t
    end
    (;snapshots,history)
end

end
