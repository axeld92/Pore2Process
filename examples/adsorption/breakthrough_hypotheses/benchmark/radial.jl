using LinearAlgebra

const OUTPUT = "results/adsorption/breakthrough_hypotheses/benchmark"
const CUTS = [.25, .5, .75, 1., 1.5, 2., 2.5]
const BINS_PER_DIAMETER = 160

function read_profile(path, columns)
    pairs = Tuple{Float64, Float64}[]
    for line in Iterators.drop(eachline(path), 2)
        fields = split(line, ',')
        isempty(strip(fields[columns[1]])) && continue
        push!(pairs, (parse(Float64, fields[columns[1]]), parse(Float64, fields[columns[2]])))
    end
    sort!(pairs)
    distances, values = first.(pairs), last.(pairs)
    @assert all(diff(distances) .> 0)
    return distances, values
end

function linear(x, y, at)
    i = clamp(searchsortedlast(x, at), 1, length(x) - 1)
    weight = clamp((at - x[i]) / (x[i + 1] - x[i]), 0, 1)
    return (1 - weight) * y[i] + weight * y[i + 1]
end

wall_area(s, radius) = 1 - (1 - s / radius)^2

function radial_grid(radius, bins_per_diameter)
    n = round(Int, radius * bins_per_diameter)
    faces = range(0, radius; length=n + 1)
    r = (faces[1:end-1] .+ faces[2:end]) ./ 2
    dr = step(faces)
    ring_area = diff(faces .^ 2)
    return (; faces, s=radius .- r, dr, ring_area, w=ring_area ./ radius^2)
end

function cumulative(grid, u, distances)
    # Integrate the piecewise constant superficial velocity over the ring area.
    areas = [0.; cumsum(reverse(grid.w))]
    flow = [0.; cumsum(reverse(grid.w .* u))]
    return [linear(areas, flow ./ flow[end], wall_area(s, grid.faces[end])) for s in distances]
end

function ergun_coefficients(porosity, reynolds)
    solid = 1 .- porosity
    return (150 .* solid .^ 2 ./ porosity .^ 3, 1.75 * reynolds .* solid ./ porosity .^ 3)
end

function flow_profile(grid, porosity, reynolds, viscosity_ratio)
    # Dimensionless flow with U=u/U_in and G=(-∂p/∂x)*d_p^2/(μ*U_in).
    # The diffusion term is -(μ_eff/μ) (r U')'/r, with symmetry at r=0 and
    # U=0 at the wall. Viscosity ratio zero gives independent Ergun annuli.
    @assert all(0 .< porosity .<= 1)
    alpha, beta = ergun_coefficients(porosity, reynolds)
    faces, dr, area = grid.faces, grid.dr, grid.ring_area
    conductance = 2 * viscosity_ratio .* faces ./ dr
    diagonal = (conductance[1:end-1] .+ conductance[2:end]) ./ area
    diagonal[end] += conductance[end] / area[end]  # Wall is one half bin away.
    lower = -conductance[2:end-1] ./ area[2:end]
    upper = -conductance[2:end-1] ./ area[1:end-1]
    transport = Tridiagonal(lower, diagonal .+ alpha, upper)
    u = zeros(length(alpha))
    g = 0.
    # Solve for U and the common pressure gradient together, holding Q fixed.
    for _ in 1:30
        resistance = transport * u .+ beta .* u .^ 2
        residual = g .- resistance
        flow_error = 1 - dot(grid.w, u)
        if maximum(abs, residual) <= 1e-9 * max(1, g) && abs(flow_error) <= 1e-12
            @assert all(u .>= 0) && g > 0
            return (; u, g, relative_force_residual=maximum(abs, residual) / g, flow_error)
        end
        jacobian = Tridiagonal(lower, transport.d .+ 2 .* beta .* u, upper)
        corrections = jacobian \ hcat(residual, ones(length(u)))
        du, response = eachcol(corrections)
        dg = (flow_error - dot(grid.w, du)) / dot(grid.w, response)
        u .+= du .+ dg .* response
        g += dg
    end
    error("Radial flow did not converge at Re=$reynolds")
end

function two_annuli(grid, porosity, cutoff)
    wall = grid.s .<= cutoff
    wall_mean = dot(grid.w[wall], porosity[wall]) / sum(grid.w[wall])
    core_mean = dot(grid.w[.!wall], porosity[.!wall]) / sum(grid.w[.!wall])
    return ifelse.(wall, wall_mean, core_mean)
end

viscosity_ratio(reynolds) = 2 * exp(.0035 * reynolds)

function rmse(value, reference)
    return sqrt(sum(abs2, value .- reference) / length(reference))
end

function verify_limit()
    grid = radial_grid(3.5, BINS_PER_DIAMETER)
    pipe = flow_profile(grid, ones(length(grid.s)), 0., 1.)
    pipe_error = maximum(abs, cumulative(grid, pipe.u, CUTS) .- wall_area.(CUTS, 3.5) .^ 2)
    @assert pipe_error < 1e-5
    ergun = flow_profile(grid, fill(.427, length(grid.s)), 10., 0.)
    @assert maximum(abs, ergun.u .- 1) < 1e-10
    return pipe_error
end

function score_mri(directory, rows, checks, reynolds, local_panel, cumulative_panel)
    local_profile = joinpath(directory, "Figure_2_$local_panel.csv")
    cumulative_profile = joinpath(directory, "Figure_2_$cumulative_panel.csv")
    porosity_distance, porosity = read_profile(local_profile, (6, 7))
    velocity_distance, velocity = read_profile(local_profile, (4, 5))
    mri_distance, mri_core = read_profile(cumulative_profile, (3, 4))
    cfd_distance, cfd_core = read_profile(cumulative_profile, (1, 2))
    push!(cfd_distance, 3.5)
    push!(cfd_core, 0.)
    mri = s -> 1 - linear(mri_distance, mri_core, s) / mri_core[1]
    cfd = s -> 1 - linear(cfd_distance, cfd_core, s) / cfd_core[1]
    grid = radial_grid(3.5, BINS_PER_DIAMETER)
    epsilon = [linear(porosity_distance, porosity, s) for s in grid.s]
    reduced = flow_profile(grid, epsilon, reynolds, viscosity_ratio(reynolds))
    zoned = two_annuli(grid, epsilon, .5)
    annuli = flow_profile(grid, zoned, reynolds, 0.)
    sheared_annuli = flow_profile(grid, zoned, reynolds, viscosity_ratio(reynolds))
    # Published local velocity samples are centers of .0651-dp bins. Integrate
    # as piecewise constant with exact cylindrical weights to audit Q(r).
    boundaries = [0; (velocity_distance[1:end-1] .+ velocity_distance[2:end]) ./ 2; 3.5]
    @assert maximum(velocity_distance) <= 3.5 && all(diff(boundaries) .> 0)
    integrated_areas = wall_area.(boundaries, 3.5)
    integrated_flow = [0; cumsum(diff(integrated_areas) .* velocity)]
    integrate_velocity = s -> linear(integrated_areas, integrated_flow ./ integrated_flow[end], wall_area(s, 3.5))
    # Compare on the published MRI sample locations inside the physical tube.
    positions = mri_distance[mri_distance .<= 3.5]
    reference = mri.(positions)
    fine = radial_grid(3.5, 2 * BINS_PER_DIAMETER)
    fine_epsilon = [linear(porosity_distance, porosity, s) for s in fine.s]
    fine_flow = flow_profile(fine, fine_epsilon, reynolds, viscosity_ratio(reynolds))
    coarse_values = cumulative(grid, reduced.u, CUTS)
    fine_values = cumulative(fine, fine_flow.u, CUTS)
    change = maximum(abs, coarse_values .- fine_values)
    @assert change < 1e-3
    ergun_values = cumulative(grid, annuli.u, CUTS)
    sheared_annuli_values = cumulative(grid, sheared_annuli.u, CUTS)
    for (index, s) in enumerate(CUTS)
        println(rows, join((reynolds, s, wall_area(s, 3.5), mri(s), cfd(s),
                            integrate_velocity(s), ergun_values[index], sheared_annuli_values[index],
                            fine_values[index]), ','))
    end
    println(checks, join((reynolds, dot(grid.w, epsilon), integrated_flow[end],
                          rmse(cfd.(positions), reference), rmse(integrate_velocity.(positions), reference),
                          rmse(cumulative(grid, annuli.u, positions), reference),
                          rmse(cumulative(grid, sheared_annuli.u, positions), reference),
                          rmse(cumulative(fine, fine_flow.u, positions), reference),
                          change, fine_flow.relative_force_residual, fine_flow.flow_error), ','))
end

function alhamed_porosity(grid, mean_porosity)
    # Kwapinski et al. Eq. (1), based on Giese et al. liquid measurements.
    smooth = 1 .+ 1.36 .* exp.(-5 .* grid.s)
    smooth .*= mean_porosity / dot(grid.w, smooth)
    # de Klerk, as reproduced in Das et al. Eq. (15). The core offset is
    # chosen so that the original bulk liquid and solid volumes are retained.
    s = grid.s
    outer = s .<= .637
    oscillation = .29 .* exp.(-.6 .* s) .* cos.(2.3 * pi .* (s .- .16)) .+ .15 .* exp.(-.9 .* s)
    shape = ifelse.(outer, 2.14 .* s .^ 2 .- 2.53 .* s .+ 1, oscillation)
    core_offset = (mean_porosity - dot(grid.w, shape)) / sum(grid.w[.!outer])
    shape[.!outer] .+= core_offset
    return (("smooth_Giese_Kwapinski", smooth), ("oscillatory_de_Klerk", shape))
end

function alhamed_scenarios(directory)
    particle_diameter = .0008
    tube_radius = .0055
    viscosity = .001 / 1.2  # Same water-viscosity assumption as the stock-TRT scenarios.
    density = 1000.
    mean_porosity = .427
    radius = tube_radius / particle_diameter
    grid = radial_grid(radius, 2 * BINS_PER_DIAMETER)
    open(joinpath(directory, "near_wall_scenarios.csv"), "w") do io
        println(io, "porosity_shape,flow_mL_min,reynolds,wall_width_dp,wall_area_fraction,wall_voidage,wall_solid_fraction,brinkman_wall_flow_fraction,two_annuli_ergun_wall_flow_fraction")
        for (name, epsilon) in alhamed_porosity(grid, mean_porosity)
            @assert abs(dot(grid.w, epsilon) - mean_porosity) <= 1e-12
            for rate in (23.3, 62.5)
                superficial_velocity = rate * 1e-6 / (60 * pi * tube_radius^2)
                reynolds = density * superficial_velocity * particle_diameter / viscosity
                reduced = flow_profile(grid, epsilon, reynolds, viscosity_ratio(reynolds))
                annuli = flow_profile(grid, two_annuli(grid, epsilon, .5), reynolds, 0.)
                for s in (.25, .5, .75, 1.)
                    wall = grid.s .<= s
                    wall_porosity = dot(grid.w[wall], epsilon[wall])
                    fraction = wall_area(s, radius)
                    solid_fraction = (fraction - wall_porosity) / (1 - mean_porosity)
                    wall_flow = only(cumulative(grid, reduced.u, [s]))
                    annuli_flow = only(cumulative(grid, annuli.u, [s]))
                    println(io, join((name, rate, reynolds, s, fraction, wall_porosity / fraction,
                                      solid_fraction, wall_flow, annuli_flow), ','))
                end
            end
        end
    end
end

function main(directory)
    @assert isdir(directory) "Run benchmark/fetch.py to get the source Figure 2 CSVs"
    pipe_error = verify_limit()
    open(joinpath(directory, "radial_benchmark.csv"), "w") do rows
        println(rows, "reynolds,wall_width_dp,wall_area_fraction,mri_wall_flow_fraction,pr_cfd_wall_flow_fraction,integrated_mri_wall_flow_fraction,two_annuli_ergun_wall_flow_fraction,two_annuli_brinkman_wall_flow_fraction,brinkman_wall_flow_fraction")
        open(joinpath(directory, "benchmark_checks.csv"), "w") do checks
            println(checks, "reynolds,area_mean_mri_voidage,integrated_mri_u_over_u_in,pr_cfd_rmse,integrated_mri_rmse,two_annuli_ergun_rmse,two_annuli_brinkman_rmse,brinkman_rmse,brinkman_grid_change,brinkman_relative_force_residual,brinkman_flow_error")
            for (reynolds, a, b) in ((10, 'a', 'b'), (30, 'c', 'd'), (50, 'e', 'f'))
                score_mri(directory, rows, checks, reynolds, a, b)
            end
        end
    end
    alhamed_scenarios(directory)
    println("Pipe-limit maximum cumulative error: ", pipe_error)
    println("Wrote radial_benchmark.csv, benchmark_checks.csv and near_wall_scenarios.csv in ", directory)
end

main(isempty(ARGS) ? OUTPUT : only(ARGS))
