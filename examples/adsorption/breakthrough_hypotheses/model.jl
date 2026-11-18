module BreakthroughHypotheses

using LinearAlgebra
using ..MemoryColumn: advect!, prepare_intervals!, schedule
using ..FilmParticle
using ..AccessCohorts

function fractions(values)
    f = Vector{Float64}(values)
    all(>(0), f) && isapprox(sum(f), 1; atol=1e-12, rtol=0) ||
        throw(ArgumentError("Expected positive fractions that sum to one"))
    f
end

function mix_regions!(incoming, flow, mixing)
    wall, core = flow
    for i in axes(incoming, 1)
        exchange = mixing * (incoming[i, 2] - incoming[i, 1])
        incoming[i, 1] += core * exchange
        incoming[i, 2] -= wall * exchange
    end
end

function solve(rules, flow, contact, zone_flow, zone_void, zone_carbon, grouping;
        blocks, terminal_sections, liquid_volume_m3, dry_mass_kg, flow_m3_s,
        inlet_concentration_kg_m3, direct_fraction, slow_ratio, mixing=0.,
        startup_s=0., startup_dt_s=.01, capacity_kg_kg=.1115, affinity_m3_kg=139.,
        linear_partition_m3_kg=nothing, kwargs...)
    p = FilmParticle.parameters(; capacity_kg_kg, affinity_m3_kg, linear_partition_m3_kg, kwargs...)
    volumes, waters, carbons = fractions(zone_void), fractions(zone_flow), fractions(zone_carbon)
    zones = length(waters)
    n, kinds = size(contact)
    Q, cin = flow_m3_s, inlet_concentration_kg_m3
    p.qinit == 0 && all(v -> isfinite(v) && v > 0, (Q, cin, liquid_volume_m3, dry_mass_kg)) ||
        throw(ArgumentError("Expected positive inventories, feed and initially clean carbon"))
    n == length(flow) && kinds == length(rules) && 0 <= terminal_sections < kinds && blocks >= 1 ||
        throw(ArgumentError("Invalid local transport or column geometry"))
    zones in (1, 2) && length(volumes) == length(carbons) == zones && 0 <= mixing <= 1 ||
        throw(ArgumentError("Invalid region split or exchange"))
    labels = Vector{Int}(grouping)
    length(labels) == n*zones && sort(unique(labels)) in ([1], [1, 2]) ||
        throw(ArgumentError("Invalid adsorption-history groups"))
    all(>=(0), contact) && all(isapprox.(vec(sum(contact; dims=1)), 1.; atol=1e-12, rtol=0)) ||
        throw(ArgumentError("Each axial section must use its full particle area once"))
    all(>(0), flow) && isapprox(sum(flow), 1; atol=1e-12, rtol=0) ||
        throw(ArgumentError("Expected positive path flow shares that sum to one"))
    groups = [findall(==(g), labels) for g in 1:maximum(labels)]
    kinds_in_column = [repeat(1:kinds, blocks); 1:terminal_sections]
    sections = length(kinds_in_column)
    counts = [count(==(kind), kinds_in_column) for kind in 1:kinds]
    weights = vcat([f .* contact for f in carbons]...)
    flows = vcat([Q*f .* flow for f in waters]...)
    mass = dry_mass_kg / sections
    area = 3mass / (p.rho*p.R)
    tau = sum(dot(rules[k].share, rules[k].seconds) for k in kinds_in_column)
    scales = liquid_volume_m3/Q/tau .* volumes ./ waters
    max_delay = maximum(scales) * maximum(maximum(r.seconds) for r in rules)
    time = schedule(p.tend, p.dt, startup_s, startup_dt_s)
    window = min(length(time), ceil(Int, max_delay/minimum(diff(time)))+3)
    slots = [mod1(k, window) for k in eachindex(time)]

    beds = [let
                 allocation = sum(view(weights, indices, k))
                 allocation == 0 ? nothing : AccessCohorts.setup(p, weights[indices, k:k]/allocation, counts[k],
                     mass*allocation, area*allocation, flows[indices], cin, direct_fraction, slow_ratio)
            end for indices in groups, k in 1:kinds]
    c, incoming = zeros(n, zones), zeros(n, zones)
    contacts = [view(vec(c), indices) for indices in groups]
    history = zeros(n, window, sections, zones)
    histories = [view(history, :, :, :, z) for z in 1:zones]
    delivered = [zeros(length(rules[k].share)) for k in kinds_in_column, _ in 1:zones]
    left = [zeros(Int, length(r.share)) for r in rules, _ in 1:zones]
    right, alpha = deepcopy(left), [zeros(length(r.share)) for r in rules, _ in 1:zones]
    outlet, passed, liquid, solid, residual = (zeros(length(time)) for _ in 1:5)
    cohorts = zeros(2, length(time))
    maxflux, uptake, release = 0., 0., 0.
    minc, maxc, maxq = 0., 0., 0.
    old_h = -1.
    op = FilmParticle.operator(p.nr, p.R, p.Ds, p.kf, time[2])
    t0 = time_ns()
    for step in 1:length(time)-1
        cur, prev = slots[step+1], slots[step]
        h = time[step+1]-time[step]
        if !isapprox(h, old_h; rtol=0., atol=8eps(p.tend))
            op = FilmParticle.operator(p.nr, p.R, p.Ds, p.kf, h)
            old_h = h
        end
        for bed in beds
            isnothing(bed) || AccessCohorts.diffuse!(bed, op)
        end
        for z in 1:zones, k in 1:kinds
            prepare_intervals!(left[k,z], right[k,z], alpha[k,z], rules[k].seconds,
                               time, step+1, slots, scales[z])
        end
        fill!(incoming, cin)
        for (section,k) in enumerate(kinds_in_column)
            zones == 2 && mix_regions!(incoming, waters, mixing)
            for z in 1:zones
                local_history = histories[z]
                if k == 1
                    boundary = dot(flow, view(incoming,:,z))
                    local_history[1,cur,section] = muladd(h, boundary, local_history[1,prev,section])
                else
                    @inbounds for i in 1:n
                        local_history[i,cur,section] = muladd(h, incoming[i,z], local_history[i,prev,section])
                    end
                end
                liquid[step+1] += Q*waters[z]*advect!(view(c,:,z), local_history, delivered[section,z],
                    rules[k], section, cur, left[k,z], right[k,z], alpha[k,z], h, flow)
            end
            for g in eachindex(groups)
                bed = beds[g,k]
                isnothing(bed) && continue
                uptake_result = AccessCohorts.contact!(bed, contacts[g], cld(section, kinds), 1, h, op,
                    capacity_kg_kg, affinity_m3_kg, linear_partition_m3_kg)
                for i in 1:2
                    cohorts[i,step+1] += uptake_result.adsorbed[i]
                end
                maxq = max(maxq, uptake_result.maxq)
                maxflux = max(maxflux, uptake_result.maxflux)
                uptake += uptake_result.uptake
                release += uptake_result.release
            end
            minc = min(minc, minimum(c))
            maxc = max(maxc, maximum(c))
            copyto!(incoming, c)
        end
        outlet[step+1] = sum(waters[z]*dot(flow,view(c,:,z)) for z in 1:zones)
        passed[step+1] = passed[step]+h*Q*outlet[step+1]
        solid[step+1] = sum(view(cohorts,:,step+1))
        residual[step+1] = liquid[step+1]+solid[step+1]+passed[step+1]-Q*cin*time[step+1]
    end
    (; time_s=time, outlet_concentration_kg_m3=outlet, cumulative_outlet_kg=passed,
       liquid_inventory_kg=liquid, adsorbed_inventory_kg=solid, mass_residual_kg=residual,
       normalized_mass_residual=maximum(abs,residual)/(Q*cin*time[end]),
       normalized_interface_flux_mismatch=maxflux/max(p.kf*cin,1e-30),
       concentration_range_kg_m3=[minc,maxc], loading_range_kg_kg=[0.,maxq],
       liquid_volume_m3, physical_adsorbent_mass_kg=dry_mass_kg, external_particle_area_m2=sections*area,
       total_film_conductance_m3_s=sections*area*p.kf,
       cumulative_film_uptake_kg=uptake, cumulative_film_release_kg=release,
       delay_scale=scales[1], maximum_delay_s=max_delay, sections,
       adsorbed_inventory_by_cohort_kg=cohorts,
       carbon_by_cohort_kg=dry_mass_kg.*[direct_fraction,1-direct_fraction],
       equivalent_film_by_cohort_m_s=p.kf.*[1.,slow_ratio],
       effective_conductance_m3_s=sections*area*p.kf*(direct_fraction+(1-direct_fraction)*slow_ratio),
       equilibrium_adsorbed_ceiling_kg=dry_mass_kg*p.qe(cin)*(slow_ratio>0 ? 1. : direct_fraction),
       wall_s=(time_ns()-t0)/1e9)
end

end
