"""Generate a reproducible synthetic hard-sphere packing, periodic in the tube axis.

Lengths in this generator are in bead diameters. Grow soft spheres in small
increments and minimize overlaps at each diameter. The final audit uses the
actual minimum-image distances and the cylindrical wall.
"""
module GlassPacking

using DelimitedFiles, LinearAlgebra, Random

function neighbours(x, period, radius)
    reach = 1.3
    nx = floor(Int, period / reach)
    nx >= 3 || error("Packing period must span at least four diameters")
    nt = ceil(Int, 2radius / reach) + 2
    buckets = [Int[] for i in 1:nx, j in 1:nt, k in 1:nt]
    cell(i) = (1 + floor(Int, mod(x[1,i], period) * nx / period),
               2 + floor(Int, (x[2,i]+radius)/reach), 2 + floor(Int, (x[3,i]+radius)/reach))
    for i in axes(x,2)
        push!(buckets[cell(i)...], i)
    end
    pairs = Tuple{Int,Int}[]
    for i in axes(x,2)
        a,b,c = cell(i)
        for da in -1:1, db in -1:1, dc in -1:1, j in buckets[mod1(a+da,nx),b+db,c+dc]
            j > i || continue
            dx = rem(x[1,i]-x[1,j], period, RoundNearest)
            dx^2 + (x[2,i]-x[2,j])^2 + (x[3,i]-x[3,j])^2 < reach^2 && push!(pairs,(i,j))
        end
    end
    pairs
end

function forces!(f, x, pairs, period, radius, diameter)
    fill!(f,0)
    overlap = 0.0
    @inbounds for (i,j) in pairs
        dx = rem(x[1,i]-x[1,j],period,RoundNearest)
        dy, dz = x[2,i]-x[2,j], x[3,i]-x[3,j]
        r2 = dx*dx+dy*dy+dz*dz
        if r2 < diameter^2
            r = sqrt(r2)
            overlap = max(overlap,diameter-r)
            scale = (diameter-r)/r
            f[1,i] += dx*scale; f[1,j] -= dx*scale
            f[2,i] += dy*scale; f[2,j] -= dy*scale
            f[3,i] += dz*scale; f[3,j] -= dz*scale
        end
    end
    @inbounds for i in axes(x,2)
        r = hypot(x[2,i], x[3,i])
        penetration = r + diameter/2 - radius
        if penetration > 0
            overlap = max(overlap,penetration)
            f[2,i] -= x[2,i]*penetration/r
            f[3,i] -= x[3,i]*penetration/r
        end
    end
    overlap
end

function shell_seed(radius, period, n, rng)
    # A confined random initial placement can jam below this solids loading.
    # Seed staggered cylindrical shells, then remove random sites and displace
    # every centre before overlap relaxation. Each seed is a new packing.
    rows = 2round(Int,period / (2*.875))
    sites = NTuple{3,Float64}[]
    for r in (radius-.5):(-.87):0
        around = max(1,floor(Int,2pi*r))
        phase, start = rand(rng), rand(rng)
        for i in 0:rows-1, j in 0:around-1
            angle = 2pi*(phase+(j+.5*(i%2))/around)
            push!(sites,((i+start)*period/rows,r*cos(angle),r*sin(angle)))
        end
    end
    length(sites) >= n || error("Too few shell sites at this period")
    selected = sites[randperm(rng,length(sites))[1:n]]
    x = reduce(hcat,collect.(selected))
    x .+= .16 .* (2rand(rng,3,n) .- 1)
    x
end

function confined_seed(radius, period, n, rng)
    # Template for the 6.75-diameter tube. At a seven-diameter period these
    # shell inventories give the target solids fraction. Wind each shell to
    # maximize its own minimum bead separation, then relax all positions.
    # Wind each cylindrical layer around the periodic box. Choose within-layer
    # spacing solely by sphere separation, before any flow or tracer calculation.
    layers = (radius-.5, radius-1.405, radius-2.33, 0.)
    base_counts = [143,93,49,7]
    counts = round.(Int,base_counts .* (period/7))
    counts[2] += n-sum(counts)
    sites = zeros(3,n)
    previous = 0
    for (r,count) in zip(layers,counts)
        distances(k) = minimum(hypot(m*period/count, 2r*sinpi(m*k/count))
                               for m in 1:ceil(Int,1.1count/period))
        winding = argmax(distances,0:div(count,2)) * rand(rng,(-1,1))
        phase, start = rand(rng),rand(rng)
        for j in 0:count-1
            angle = 2pi*(phase+winding*j/count)
            sites[:,previous+j+1] .= ((j+start)*period/count,r*cos(angle),r*sin(angle))
        end
        previous += count
    end
    sites .+= .05 .* (2rand(rng,3,n) .- 1)
    sites
end

function grow(radius, period, porosity, seed)
    n = round(Int, 6radius^2*period*(1-porosity))
    rng = Xoshiro(seed)
    x = radius < 4 ? confined_seed(radius,period,n,rng) : shell_seed(radius,period,n,rng)
    v, f, anchor = zero(x), zero(x), copy(x)
    pairs = neighbours(x,period,radius)
    for d in (.65,.80,.90,.95,.975,.99,1.)
        fill!(v,0)
        dt, mix, positive = .08, .1, 0
        best_overlap = Inf
        best = copy(x)
        for step in 1:90000
            if maximum(i->sum(k->(x[k,i]-anchor[k,i])^2,1:3),1:n) > .14^2
                pairs = neighbours(x,period,radius)
                copyto!(anchor,x)
            end
            overlap = forces!(f,x,pairs,period,radius,d)
            overlap < 1e-7 && break
            step < 90000 || error("Packing stalled at diameter $d with overlap $overlap")
            if step % 1800 == 0 && overlap > 1e-5
                if overlap < best_overlap
                    best_overlap = overlap
                    copyto!(best,x)
                else
                    copyto!(x,best)
                end
                x .+= .12 .* randn(rng,size(x))
                fill!(v,0)
                dt, mix, positive = .08, .1, 0
                continue
            end
            @. v += dt*f
            if dot(v,f) > 0
                positive += 1
                ratio = norm(v)/norm(f)
                @. v = (1-mix)*v + mix*ratio*f
                if positive > 5
                    dt = min(1.1dt,.30)
                    mix *= .99
                end
            else
                positive = 0
                fill!(v,0)
                dt *= .5
                mix = .1
            end
            @. x += dt*v
        end
    end
    x[1,:] .= mod.(x[1,:],period)
    overlap = forces!(f,x,neighbours(x,period,radius),period,radius,1.)
    overlap < 1e-7 || error("Final wall or bead overlaps exceed tolerance")
    println("spheres=$n porosity=$(1-n/(6radius^2*period)) max_overlap_dp=$overlap")
    x
end

if abspath(PROGRAM_FILE) == @__FILE__
    radius, period, porosity = parse.(Float64,ARGS[1:3])
    BLAS.set_num_threads(1)
    x = grow(radius, period, porosity, parse(Int,ARGS[4]))
    writedlm(ARGS[5],permutedims(x))
end
end
