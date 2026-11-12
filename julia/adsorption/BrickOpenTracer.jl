module BrickOpenTracer

using ..OpenBoundedFaceTransport: exit_rate

export prepare, evolve

chunks(n) = [((i-1)*n÷Threads.nthreads()+1):(i*n÷Threads.nthreads()) for i in 1:Threads.nthreads()]

"""Reconstruct tracer neighbours and face rates within each liquid voxel."""
function prepare(before,after,inflow,outflow,ratio,dx,diffusivity)
    parents = size(before,2)
    width = ratio^3
    count = parents*width
    strides = (ratio^2,ratio,1)
    backward = zeros(Int32,3,count)
    forward = similar(backward)
    rate = Array{Float64}(undef,ratio,3,parents)
    boundary = zeros(count)
    offset = Matrix{Int}(undef,3,width)
    for s in 1:width, a in 1:3
        offset[a,s] = (s-1)÷strides[a]%ratio+1
    end
    face_divisor = ratio^2
    volume = dx^3
    d = diffusivity/dx^2
    parent_chunks = chunks(parents)
    block_max = zeros(length(parent_chunks))
    Threads.@threads :static for block in eachindex(parent_chunks)
        largest = 0.
        @inbounds for p in parent_chunks[block]
            for a in 1:3, o in 1:ratio
                s = o/ratio
                rate[o,a,p] = (((1-s)*inflow[a,p]+s*outflow[a,p])/face_divisor)/volume
            end
            for s in 1:width
                i = (p-1)*width+s
                depletion = 0.
                for a in 1:3
                    o,stride = offset[a,s],strides[a]
                    if o>1
                        backward[a,i] = i-stride
                        qin = rate[o-1,a,p]
                    else
                        other = before[a,p]
                        backward[a,i] = other>0 ? (other-1)*width+s+(ratio-1)*stride : 0
                        qin = (inflow[a,p]/face_divisor)/volume
                        if other==0 && a==1
                            boundary[i] += qin
                        end
                    end
                    if o<ratio
                        forward[a,i] = i+stride
                    else
                        other = after[a,p]
                        forward[a,i] = other>0 ? (other-1)*width+s-(ratio-1)*stride : 0
                        if other==0 && a==1
                            boundary[i] -= rate[o,a,p]
                        end
                    end
                    depletion += max(-qin,0.)+max(rate[o,a,p],0.)
                    depletion += d*((backward[a,i]!=0)+(forward[a,i]!=0))
                end
                largest = max(largest,depletion)
            end
        end
        block_max[block] = largest
    end
    (;ratio,width,count,rate,d,volume,offset,backward,forward,boundary,parent_chunks,
      cell_chunks=chunks(count),max_depletion=maximum(block_max))
end

function face_increments!(c,grid,low,anti,dt)
    Threads.@threads :static for block in eachindex(grid.parent_chunks)
        @inbounds for p in grid.parent_chunks[block], s in 1:grid.width
            i = (p-1)*grid.width+s
            a = c[i]
            for direction in 1:3
                j = grid.forward[direction,i]
                if j==0
                    low[direction,i] = 0.
                    anti[direction,i] = 0.
                    continue
                end
                b = c[j]
                q = grid.rate[grid.offset[direction,s],direction,p]
                donor = ifelse(q>=0,a,b)
                low[direction,i] = dt*(q*donor+grid.d*(a-b))
                high = (a+b)/2
                up,down = grid.backward[direction,i],grid.forward[direction,j]
                if up!=0 && down!=0
                    high += (a+b-c[up]-c[down])/12
                end
                anti[direction,i] = dt*q*(high-donor)
            end
        end
    end
end

@inline function accumulate(updated,lower,upper,positive,negative,c,q,delta)
    updated+q,min(lower,c),max(upper,c),positive+max(delta,0.),negative+max(-delta,0.)
end

function factors!(out,c,grid,low,anti,up,down,feed,dt)
    Threads.@threads :static for block in eachindex(grid.cell_chunks)
        @inbounds for i in grid.cell_chunks[block]
            updated = c[i]
            lower,upper = updated,updated
            positive,negative = 0.,0.
            for a in 1:3
                j,k = grid.backward[a,i],grid.forward[a,i]
                if j!=0 && j<i
                    updated,lower,upper,positive,negative = accumulate(
                        updated,lower,upper,positive,negative,c[j],low[a,j],anti[a,j])
                end
                if k!=0
                    updated,lower,upper,positive,negative = accumulate(
                        updated,lower,upper,positive,negative,c[k],-low[a,i],-anti[a,i])
                end
                if j>i
                    updated,lower,upper,positive,negative = accumulate(
                        updated,lower,upper,positive,negative,c[j],low[a,j],anti[a,j])
                end
            end
            q = grid.boundary[i]
            if q>0
                updated += dt*q*feed
                lower = min(lower,feed)
                upper = max(upper,feed)
            elseif q<0
                updated += dt*q*c[i]
            end
            out[i] = updated
            up[i] = positive>0 ? clamp((upper-updated)/positive,0.,1.) : 1.
            down[i] = negative>0 ? clamp((updated-lower)/negative,0.,1.) : 1.
        end
    end
end

@inline function correction(delta,i,j,up,down)
    factor = delta>=0 ? min(down[i],up[j]) : min(up[i],down[j])
    factor*delta
end

function apply!(out,grid,anti,up,down)
    Threads.@threads :static for block in eachindex(grid.cell_chunks)
        @inbounds for i in grid.cell_chunks[block]
            updated = out[i]
            for a in 1:3
                j,k = grid.backward[a,i],grid.forward[a,i]
                if j!=0 && j<i
                    updated += correction(anti[a,j],j,i,up,down)
                end
                if k!=0
                    updated -= correction(anti[a,i],i,k,up,down)
                end
                if j>i
                    updated += correction(anti[a,j],j,i,up,down)
                end
            end
            out[i] = updated
        end
    end
end

function step!(out,c,grid,work,feed,dt)
    face_increments!(c,grid,work.low,work.anti,dt)
    factors!(out,c,grid,work.low,work.anti,work.up,work.down,feed,dt)
    apply!(out,grid,work.anti,work.up,work.down)
end

function inventory(c,grid,sums,lows,highs)
    Threads.@threads :static for block in eachindex(grid.cell_chunks)
        cells = view(c,grid.cell_chunks[block])
        sums[block] = sum(cells)
        lows[block],highs[block] = extrema(cells)
    end
    grid.volume*sum(sums),minimum(lows),maximum(highs)
end

"""Integrate the open tracer with the same corrections on a structured face layout."""
function evolve(grid,initial,dt,pulse_steps,steps; stream=nothing,progress=false)
    0<dt*grid.max_depletion<1 || error("Time step exceeds the monotone limit")
    c = copy(initial)
    a,b = similar(c),similar(c)
    work = (;low=zeros(3,grid.count),anti=zeros(3,grid.count),up=similar(c),down=similar(c))
    introduced = removed = m1 = m2 = 0.
    inlet = grid.volume*sum(x->max(x,0.),grid.boundary)
    exit = findall(<(0),grid.boundary)
    rates = -grid.boundary[exit]
    sums,lows,highs = (zeros(length(grid.cell_chunks)) for _ in 1:3)
    history = NTuple{9,Float64}[]
    function record(t,dm)
        liquid,lo,hi = inventory(c,grid,sums,lows,highs)
        all(isfinite,(liquid,lo,hi,removed,m1,m2)) || error("Non-finite tracer state")
        row = (t,liquid,lo,hi,introduced,removed,m1,m2,dm)
        push!(history,row)
        stream===nothing || println(stream,join(row,','))
    end
    record(0.,0.)
    for step in 1:steps
        feed = step<=pulse_steps ? 1. : 0.
        t = (step-1)*dt
        f1 = exit_rate(c,exit,rates,grid.volume)
        step!(a,c,grid,work,feed,dt)
        f2 = exit_rate(a,exit,rates,grid.volume)
        step!(b,a,grid,work,feed,dt)
        @. b = c+(b-c)/4
        f3 = exit_rate(b,exit,rates,grid.volume)
        step!(a,b,grid,work,feed,dt)
        @. c += 2*(a-c)/3
        t2,t3 = t+dt,t+dt/2
        dm = dt*(f1+f2+4*f3)/6
        introduced += dt*inlet*feed
        removed += dm
        m1 += dt*(f1*t+f2*t2+4*f3*t3)/6
        m2 += dt*(f1*t^2+f2*t2^2+4*f3*t3^2)/6
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
