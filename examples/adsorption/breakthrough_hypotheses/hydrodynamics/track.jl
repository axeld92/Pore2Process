# Inert Julia parcels on the balanced OpenLB faces. Record the first absorbing
# carbon encounter on the same trajectories, with the tube wall reflecting.
include(joinpath(@__DIR__, "../../carbon_contact/track_contacts.jl"))
using .TrackContacts: Contacts, FaceEvent, Flow, diffusion_contacts
using .TrackContacts.Track: crossing, destination, drift
using Random

const WALL = UInt8(1)
const CORE = UInt8(2)
const CARBON = UInt8(2)
const TUBE = UInt8(3)
const CUT = UInt8(4)

struct Grid
    flow::Flow
    material::Matrix{UInt8}
    face_wall::Matrix{UInt8}
    coordinates::Matrix{Int32}
    centre_cells::Float64
    cut_radius_cells::Float64
end

@inline function zone(grid, x, voxel)
    @inbounds squared = (grid.coordinates[2,voxel] + x[2] - grid.centre_cells)^2 +
                       (grid.coordinates[3,voxel] + x[3] - grid.centre_cells)^2
    squared >= grid.cut_radius_cells^2 ? WALL : CORE
end

mutable struct Exposure
    grid::Grid
    time::Float64
    first_carbon::Float64
    carbon_zone::UInt8
    first_tube::Float64
    wall_carbon_local_m::Float64
    core_carbon_local_m::Float64
end

@inline function Base.push!(e::Exposure, event::FaceEvent)
    face, voxel = mod1(Int(event.face_id), 6), cld(Int(event.face_id), 6)
    @inbounds material = e.grid.material[face, voxel]
    if material == CARBON
        in_wall = e.grid.face_wall[face,voxel] == 1
        if isnan(e.first_carbon)
            e.first_carbon = e.time
            e.carbon_zone = in_wall ? WALL : CORE
        end
        if in_wall
            e.wall_carbon_local_m += event.reflection_m
        else
            e.core_carbon_local_m += event.reflection_m
        end
    elseif material == TUBE
        isnan(e.first_tube) && (e.first_tube = e.time)
    elseif material != CUT
        error("Reflection is not on a carbon, tube or inlet/outlet face")
    end
    e
end

# Use the established crossings and analytic drift within one pore voxel.
# The physical control section can have a small measured reflux at either end.
@inline function advance(grid, x, voxel, t, dt)
    stop = t + dt
    flow = grid.flow
    while t < stop
        @inbounds begin
            times = ntuple(a -> crossing(x[a], flow.low[a,voxel], flow.high[a,voxel], flow.slope[a,voxel]), 3)
            axis = times[1][1] <= times[2][1] ? 1 : 2
            axis = times[3][1] < times[axis][1] ? 3 : axis
            hit, side = times[axis]
            duration = min(hit, stop-t)
            x = ntuple(a -> drift(x[a], flow.low[a,voxel], flow.high[a,voxel], flow.slope[a,voxel], duration), 3)
            if hit > stop-t
                return x, voxel, 0, stop
            end
            t += hit
            next = destination(flow, voxel, axis, side)
            if next == 0
                axis == 1 && grid.material[side > 0 ? 2 : 1,voxel] == CUT || error("Advective exit into solid")
                return x, voxel, side, t
            end
            voxel = Int(next)
            x = Base.setindex(x, side > 0 ? 0. : 1., axis)
        end
    end
    x, voxel, 0, stop
end

@inline expose(e, initial_zone, final_zone, exit, time, wall_s, surviving_wall_s) =
    (initial_zone, final_zone, exit, time, e.first_carbon, e.carbon_zone, e.first_tube,
     e.wall_carbon_local_m, e.core_carbon_local_m, wall_s, surviving_wall_s)

function trace(grid, x, voxel, rng, dx, sigma, dt, steps)
    initial_zone = zone(grid, x, voxel)
    e = Exposure(grid, 0., NaN, 0, NaN, 0., 0.)
    wall_s = 0.
    surviving_wall_s = 0.
    half = dt / 2
    for step in 1:steps
        t = (step-1)*dt
        x, voxel, exit, endtime = advance(grid, x, voxel, t, half)
        in_wall = zone(grid, x, voxel) == WALL
        entry_wall_time = (endtime-t) * in_wall
        wall_s += entry_wall_time
        isnan(e.first_carbon) && (surviving_wall_s += entry_wall_time)
        exit != 0 && return expose(e, initial_zone, zone(grid, x, voxel), exit, endtime, wall_s, surviving_wall_s)
        e.time = t + half
        stats = Contacts()
        for axis in 1:3
            x, voxel, _, stats = diffusion_contacts(grid.flow, x, voxel, NaN, axis, sigma*randn(rng),
                                                     e.time, dx, stats, e, 1)
        end
        x, voxel, exit, endtime = advance(grid, x, voxel, e.time, half)
        in_wall_after = zone(grid, x, voxel) == WALL
        later_wall_time = (endtime-e.time) * in_wall_after
        wall_s += later_wall_time
        isnan(e.first_carbon) && (surviving_wall_s += later_wall_time)
        exit != 0 && return expose(e, initial_zone, zone(grid, x, voxel), exit, endtime, wall_s, surviving_wall_s)
    end
    expose(e, initial_zone, zone(grid, x, voxel), 0, steps*dt, wall_s, surviving_wall_s)
end

function main(args)
    length(args) == 10 || error("Use directory output dx D dt stop n seed transverse_width wall_radius_m")
    directory, output = args[1:2]
    dx, D, dt, stop = parse.(Float64, args[3:6])
    count, seed, width = parse.(Int, args[7:9])
    cut = parse(Float64, args[10])
    dx > 0 && D > 0 && dt > 0 && stop >= dt && count > 0 && cut > 0 && width > 0 ||
        error("Invalid carrier grid or tracer parameters")
    load(name, T, n) = reshape(reinterpret(T, read(joinpath(directory, name))), n, :)
    before = load("before.i32", Int32, 3)
    after = load("after.i32", Int32, 3)
    low = load("incoming.f64", Float64, 3) ./ dx^3
    high = load("outgoing.f64", Float64, 3) ./ dx^3
    coordinates = load("coords.i32", Int32, 3)
    material = load("material.u8", UInt8, 6)
    wall = load("carbon_face_in_shell.u8", UInt8, 6)
    n = size(before, 2)
    size(after) == size(low) == size(high) == size(coordinates) == size(before) &&
        size(material) == size(wall) == (6, n) || error("Mismatched compact face fields")
    # The established Brownian folding routine receives the whole period here.
    # Its contact observer sees physical voxel IDs, and no internal axial seam.
    flow = Flow(before, after, low, high, high-low, n)
    grid = Grid(flow, material, wall, coordinates, width/2, cut/dx)
    sources = findall(p -> material[1,p] == CUT && low[1,p] > 0, 1:n)
    isempty(sources) && error("Invalid selected inlet plane")
    weights = cumsum(low[1,sources])
    steps = round(Int, stop/dt)
    sigma = sqrt(2D*dt) / dx
    result = Matrix{Float64}(undef, 11, count)
    Threads.@threads for marker in 1:count
        rng = Xoshiro(seed + marker)
        voxel = sources[searchsortedfirst(weights, rand(rng)*last(weights))]
        result[:,marker] .= trace(grid, (0., rand(rng), rand(rng)), voxel, rng, dx, sigma, dt, steps)
    end
    write(output, result)
    println("Traced ", count, " markers with ", sum(result[3,:] .== 1), " outlet and ",
            sum(result[3,:] .== -1), " inlet returns")
end

main(ARGS)
