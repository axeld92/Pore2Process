module ThreadedBoundedFaceTransport

using ..BoundedFaceTransport: extended_stencil, stability_cap

export prepare, workspace, advance!, evolve_threaded

function chunks(n)
    count = Threads.nthreads()
    [((i-1)*n÷count+1):(i*n÷count) for i in 1:count]
end

function prepare(left, right, axis, flow, count, dx, diffusivity)
    rate = flow/dx^3
    d = diffusivity/dx^2
    cap = stability_cap(left,right,rate,d,count)
    before, after = extended_stencil(left,right,axis,count)
    # Edge IDs follow the original serial face traversal. Positive IDs
    # indicate the left cell, negative IDs the right cell.
    edges = zeros(Int32,6,count)
    degree = zeros(UInt8,count)
    for e in eachindex(left)
        i, j = left[e], right[e]
        degree[i] += 1
        degree[j] += 1
        edges[degree[i],i] = e
        edges[degree[j],j] = -e
    end
    face_chunks = chunks(length(flow))
    cell_chunks = chunks(count)
    (;left,right,rate,d,cap,before,after,edges,degree,face_chunks,cell_chunks)
end

function workspace(c, grid)
    (;low=zeros(length(grid.rate)), anti=zeros(length(grid.rate)), up=similar(c),
      down=similar(c), energy=zeros(length(grid.face_chunks)), a=similar(c), b=similar(c))
end

function face_increments!(c,grid,work,dt)
    Threads.@threads :static for block in eachindex(grid.face_chunks)
        energy = 0.
        @inbounds for e in grid.face_chunks[block]
            a,b = c[grid.left[e]],c[grid.right[e]]
            q = grid.rate[e]
            diff = a-b
            energy += diff*diff
            donor = ifelse(q>=0,a,b)
            work.low[e] = dt*(q*donor+grid.d*diff)
            high = (a+b)/2
            up, down = grid.before[e],grid.after[e]
            if up!=0 && down!=0
                high += (a+b-c[up]-c[down])/12
            end
            work.anti[e] = dt*q*(high-donor)
        end
        work.energy[block] = energy
    end
end

function correction_factors!(out,c,grid,work)
    Threads.@threads :static for block in eachindex(grid.cell_chunks)
        @inbounds for i in grid.cell_chunks[block]
            updated = c[i]
            lower, upper = updated,updated
            positive, negative = 0.,0.
            for slot in 1:grid.degree[i]
                signed = grid.edges[slot,i]
                e = abs(signed)
                neighbor = signed>0 ? grid.right[e] : grid.left[e]
                lower = min(lower,c[neighbor])
                upper = max(upper,c[neighbor])
                if signed>0
                    updated -= work.low[e]
                    delta = -work.anti[e]
                else
                    updated += work.low[e]
                    delta = work.anti[e]
                end
                positive += max(delta,0.)
                negative += max(-delta,0.)
            end
            out[i] = updated
            work.up[i] = positive>0 ? clamp((upper-updated)/positive,0.,1.) : 1.
            work.down[i] = negative>0 ? clamp((updated-lower)/negative,0.,1.) : 1.
        end
    end
end

function apply_corrections!(out,grid,work)
    Threads.@threads :static for block in eachindex(grid.cell_chunks)
        @inbounds for i in grid.cell_chunks[block]
            updated = out[i]
            for slot in 1:grid.degree[i]
                signed = grid.edges[slot,i]
                e = abs(signed)
                delta = work.anti[e]
                if signed>0
                    j = grid.right[e]
                    factor = delta>=0 ? min(work.down[i],work.up[j]) : min(work.up[i],work.down[j])
                    updated -= factor*delta
                else
                    j = grid.left[e]
                    factor = delta>=0 ? min(work.down[j],work.up[i]) : min(work.up[j],work.down[i])
                    updated += factor*delta
                end
            end
            out[i] = updated
        end
    end
end

function limited_step!(out,c,grid,work,dt)
    face_increments!(c,grid,work,dt)
    correction_factors!(out,c,grid,work)
    apply_corrections!(out,grid,work)
    grid.d*sum(work.energy)/length(c)
end

function advance!(c,grid,work,dt)
    a,b = work.a,work.b
    d1 = limited_step!(a,c,grid,work,dt)
    d2 = limited_step!(b,a,grid,work,dt)
    @. b = c+(b-c)/4
    d3 = limited_step!(a,b,grid,work,dt)
    @. c += 2*(a-c)/3
    dt*(d1+d2+4d3)/6
end

function evolve_threaded(grid,initial,times,counts; progress=false)
    c = copy(initial)
    work = workspace(c,grid)
    history = NTuple{5,Float64}[]
    snapshots = Vector{Float64}[]
    dissipated = 0.
    function record(t)
        lo,hi = extrema(c)
        all(isfinite,(lo,hi,dissipated)) || error("Non-finite tracer state")
        push!(history,(t,sum(c)/length(c),lo,hi,dissipated))
    end
    record(0.)
    previous = 0.
    done = 0
    for (t,steps) in zip(times,counts)
        dt = (t-previous)/steps
        0 < dt <= grid.cap || error("Matched time step exceeds the stable cap")
        for step in 1:steps
            dissipated += advance!(c,grid,work,dt)
            record(step==steps ? t : previous+step*dt)
            done += 1
            if progress && (done%64==0 || step==steps)
                println("Tracer steps ",done," / ",sum(counts))
                flush(stdout)
            end
        end
        push!(snapshots,copy(c))
        previous = t
    end
    (;snapshots,history)
end

end
