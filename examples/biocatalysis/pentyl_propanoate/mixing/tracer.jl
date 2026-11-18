using Random

include(joinpath(@__DIR__, "../../../adsorption/path_exchange/track.jl"))
using .Track: Flow, advect, diffusion_axis

@inline function axial_step(x0, p0, x1, p1, indices, width)
    change = (indices[p1] - indices[p0]) + (x1[1] - x0[1])
    change - width * round(change / width)
end

"""Track an inert flux-weighted marker on the saved, conservative periodic field."""
function track_marker(flow, x, p, indices, width, length_cells, rng, sigma, dt, steps)
    distance = 0.
    returned = false
    for step in 1:steps
        for part in (1, 2, 3)
            time = (step - 1) * dt + (part > 1) * dt / 2
            if part == 2
                xnew, pnew = x, p
                for axis in 1:3
                    xnew, pnew, _ = diffusion_axis(
                        flow, xnew, pnew, Inf, axis, sigma * randn(rng), time)
                end
            else
                xnew, pnew, _, exit = advect(flow, x, p, Inf, 0., dt / 2)
                isnan(exit) || error("Axial wraps must remain connected in a periodic graph")
            end
            increment = axial_step(x, p, xnew, pnew, indices, width)
            if distance + increment >= length_cells
                arrive = part == 2 ? time : time + (dt / 2) * (length_cells - distance) / increment
                return arrive, Float64(returned)
            end
            distance += increment
            returned |= distance < -1e-8
            x, p = xnew, pnew
        end
    end
    NaN, Float64(returned)
end

function sample(flow, indices, width, cut, dx, length, D, dt, stop, count, seed)
    size(flow.low, 2) == Base.length(indices) || error("Incomplete periodic voxel index")
    all(i -> 0 <= i < width, indices) || error("Invalid axial index")
    every_speed = maximum(abs, flow.low) + maximum(abs, flow.high)
    # Track's advection integrates exactly within each voxel. Its net motion
    # between samples must stay inside half a period to unwrap uniquely.
    every_speed * dt <= width / 4 || error("Advection step is too large to unwrap")
    sigma = sqrt(2 * D * dt) / dx
    12sigma < 1 || error("Brownian step is too large to resolve reflection")
    cutoff = length / dx
    steps = round(Int, stop / dt)
    inlet = findall(i -> indices[i] == cut && flow.low[1, i] > 0., eachindex(indices))
    isempty(inlet) && error("No forward liquid flux at the release plane")
    cdf = cumsum(flow.low[1, inlet])
    result = Matrix{Float64}(undef, 2, count)
    Threads.@threads for i in 1:count
        rng = Xoshiro(seed + i)
        p = inlet[searchsortedfirst(cdf, cdf[end] * rand(rng))]
        x = (0., rand(rng), rand(rng))
        result[:, i] .= track_marker(flow, x, p, indices, width, cutoff, rng, sigma, dt, steps)
    end
    result
end

function main(args)
    length(args) == 11 || error("Use compact output width cut dx L D dt stop markers seed")
    folder, path = args[1:2]
    width, cut = parse.(Int, args[3:4])
    dx, L, D, dt, stop = parse.(Float64, args[5:9])
    count, seed = parse.(Int, args[10:11])
    0 <= cut < width && width > 0 && count > 0 || error("Invalid release cut or marker count")
    all(isfinite, (dx, L, D, dt, stop)) && min(dx, L, D, dt) > 0 && stop >= dt ||
        error("Invalid tracer geometry or times")
    Base.ENDIAN_BOM == 0x04030201 || error("Expected little-endian storage")
    load(name, T) = reinterpret(T, read(joinpath(folder, name)))
    before = reshape(load("before.i32", Int32), 3, :)
    after = reshape(load("after.i32", Int32), 3, :)
    incoming = reshape(load("incoming.f64", Float64), 3, :) ./ dx^3
    outgoing = reshape(load("outgoing.f64", Float64), 3, :) ./ dx^3
    size(before) == size(after) == size(incoming) == size(outgoing) ||
        error("Incomplete periodic velocity graph")
    indices = load("x_index.i32", Int32)
    result = sample(Flow(before, after, incoming, outgoing), indices, width, cut, dx, L, D,
                    dt, stop, count, seed)
    write(path, result)
    println("Recovered ", sum(isfinite, view(result, 1, :)), " of ", count, " blue dextran markers")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
