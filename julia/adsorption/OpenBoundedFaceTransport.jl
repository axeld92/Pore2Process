module OpenBoundedFaceTransport

using ..ThreadedBoundedFaceTransport: face_increments!, apply_corrections!, workspace

export evolve_open

function boundary_cap(grid, boundary)
    outgoing = max.(-boundary,0.)
    for e in eachindex(grid.rate)
        q = grid.rate[e]
        outgoing[grid.left[e]] += grid.d+max(q,0.)
        outgoing[grid.right[e]] += grid.d+max(-q,0.)
    end
    .8/maximum(outgoing)
end

function correction_factors!(out,c,grid,work,boundary,feed,dt)
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
            q = boundary[i]
            if q>0
                updated += dt*q*feed
                lower = min(lower,feed)
                upper = max(upper,feed)
            elseif q<0
                updated += dt*q*c[i]
            end
            out[i] = updated
            work.up[i] = positive>0 ? clamp((upper-updated)/positive,0.,1.) : 1.
            work.down[i] = negative>0 ? clamp((updated-lower)/negative,0.,1.) : 1.
        end
    end
end

function limited_open!(out,c,grid,work,boundary,feed,dt)
    face_increments!(c,grid,work,dt)
    correction_factors!(out,c,grid,work,boundary,feed,dt)
    apply_corrections!(out,grid,work)
end

function exit_rate(c,indices,rates,volume)
    total = 0.
    @inbounds @simd for j in eachindex(indices)
        total += c[indices[j]]*rates[j]
    end
    volume*total
end

"""Advance a finite rectangular pulse with conservative inlet and outlet fluxes."""
function evolve_open(grid,initial,boundary,dx,dt,pulse_steps,steps; progress=false,stream=nothing)
    length(initial)==length(boundary)==length(grid.degree) || throw(DimensionMismatch())
    dt>0 && dt<=boundary_cap(grid,boundary) || error("Time step above open stability cap")
    0<=pulse_steps<=steps || error("Pulse must fit in the observation window")
    c = copy(initial)
    work = workspace(c,grid)
    a,b = work.a,work.b
    volume = dx^3
    inlet = sum(x->max(x,0.),boundary)*volume
    exit = findall(<(0),boundary)
    rates = -boundary[exit]
    inlet>0 && !isempty(exit) || error("Inlet and outlet are required")
    introduced = removed = moment1 = moment2 = 0.
    history = NTuple{9,Float64}[]
    function record(time,dm)
        lo,hi = extrema(c)
        all(isfinite,(lo,hi,removed,moment1,moment2)) || error("Non-finite tracer state")
        row = (time,volume*sum(c),lo,hi,introduced,removed,moment1,moment2,dm)
        push!(history,row)
        if stream!==nothing
            println(stream,join(row,','))
        end
    end
    record(0.,0.)
    for step in 1:steps
        feed = step<=pulse_steps ? 1. : 0.
        t = (step-1)*dt
        f1 = exit_rate(c,exit,rates,volume)
        limited_open!(a,c,grid,work,boundary,feed,dt)
        f2 = exit_rate(a,exit,rates,volume)
        limited_open!(b,a,grid,work,boundary,feed,dt)
        @. b = c+(b-c)/4
        f3 = exit_rate(b,exit,rates,volume)
        limited_open!(a,b,grid,work,boundary,feed,dt)
        @. c += 2*(a-c)/3
        t2,t3 = t+dt,t+dt/2
        dm = dt*(f1+f2+4*f3)/6
        introduced += dt*inlet*feed
        removed += dm
        moment1 += dt*(f1*t+f2*t2+4*f3*t3)/6
        moment2 += dt*(f1*t^2+f2*t2^2+4*f3*t3^2)/6
        record(step*dt,dm)
        if progress && (step%128==0 || step==steps)
            println("Outlet tracer steps ",step," / ",steps,", recovered ",removed/introduced)
            flush(stdout)
            stream===nothing || flush(stream)
        end
    end
    (;history,concentration=c)
end

end
