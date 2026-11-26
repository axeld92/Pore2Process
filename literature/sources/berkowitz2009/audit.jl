module BerkowitzAudit

using SHA
using Statistics
using TOML

const HERE = @__DIR__

"""Least-squares conversion of vector coordinates using every labelled major tick."""
function axis_map(points, labels)
    x = Float64.(points)
    y = Float64.(labels)
    slope = sum((x .- mean(x)) .* (y .- mean(y))) / sum(abs2, x .- mean(x))
    intercept = mean(y) - slope * mean(x)
    return (; slope, intercept, max_tick_error=maximum(abs.(slope .* x .+ intercept .- y)))
end
physical(axis, point) = axis.intercept + axis.slope * point

function read_series(series, axes, tolerance_pt)
    points = sort(series["points_pt"]; by=first)
    groups = Vector{Vector{Int}}()
    for (index, point) in pairs(points)
        if isempty(groups) || point[1] - points[first(last(groups))][1] > tolerance_pt
            push!(groups, Int[])
        end
        push!(last(groups), index)
    end
    times = [mean(physical(axes.time, points[i][1]) for i in group) for group in groups]
    values = [[physical(axes.ec, points[i][2]) for i in group] for group in groups]
    low, high = minimum.(values), maximum.(values)
    return (; id=series["id"], flow=series["flow_ml_min"], direction=series["direction"],
            times, low, high, middle=(low + high) / 2, markers=length(points),
            multiple_count=count(g -> length(g) > 1, groups))
end

function lerp_cross(x1, x2, y1, y2, level)
    return x1 + (level - y1) * (x2 - x1) / (y2 - y1)
end

"""Width of the connected above-threshold interval containing the main peak."""
function crossing_width(times, signal, level)
    p = argmax(signal)
    rise = findprev(<(level), signal, p)
    fall = findnext(<(level), signal, p)
    if isnothing(rise) || isnothing(fall)
        throw(ArgumentError("The plotted time window must bracket the half-height crossings"))
    end
    left = lerp_cross(times[rise], times[rise+1], signal[rise], signal[rise+1], level)
    right = lerp_cross(times[fall-1], times[fall], signal[fall-1], signal[fall], level)
    return (; left, right, width=right-left)
end

half_height(times, signal, background) = crossing_width(times, signal, (maximum(signal)+background)/2)

"""First moment over the plotted window after clipping the sample values to background."""
function visible_mean(times, signal, background)
    y = max.(signal .- background, 0.0)
    dt = diff(times)
    area = dt .* (y[1:end-1] .+ y[2:end]) ./ 2
    moment = area .* (times[1:end-1] .+ times[2:end]) ./ 2 .+ dt.^2 .* diff(y) ./ 12
    return sum(moment) / sum(area)
end

function interpolate(times, values, t)
    first(times) <= t <= last(times) || throw(DomainError(t, "Outside the plotted window"))
    i = min(searchsortedlast(times, t), length(times) - 1)
    return values[i] + (t-times[i]) * (values[i+1]-values[i]) / (times[i+1]-times[i])
end

function screen_series(trace, config, common_window)
    t = trace.times
    y = trace.middle
    tie = config["peak_tie_uS"]
    # Every retained time whose upper reading could beat all lower readings.
    possible = t[trace.high .>= maximum(trace.low) - tie]
    nominal = t[y .>= maximum(y) - tie]
    background = config["nominal_background_uS"]
    half = half_height(t, y, background)
    backgrounds = config["background_scenarios_uS"]
    # Every admissible interpolated curve lies between low and high. The most
    # demanding threshold on the lower curve and least demanding threshold on
    # the upper curve enclose the width for every admissible main peak.
    inner = crossing_width(t, trace.low, (maximum(trace.high)+maximum(backgrounds))/2)
    outer = crossing_width(t, trace.high, (maximum(trace.low)+minimum(backgrounds))/2)
    inner.left < minimum(possible) <= maximum(possible) < inner.right ||
        error("The peak candidates span more than one half-height lobe")
    means = [visible_mean(t, series, b) for series in (trace.low, y, trace.high)
             for b in config["background_scenarios_uS"]]
    start, stop = common_window
    common = vcat(start, t[(start .< t) .& (t .< stop)], stop)
    common_y = [interpolate(t, y, time) for time in common]
    return Dict("id"=>trace.id, "direction"=>trace.direction, "flow_ml_min"=>trace.flow,
        "markers"=>trace.markers, "sample_times"=>length(t), "times_with_multiple_markers"=>trace.multiple_count,
        "max_same_time_ec_range_uS"=>maximum(trace.high - trace.low),
        "visible_window_min"=>[first(t), last(t)], "first_last_uS"=>[first(y), last(y)],
        "peak_midrange_uS"=>maximum(y), "peak_midrange_time_min"=>collect(extrema(nominal)),
        "possible_peak_time_min"=>collect(extrema(possible)),
        "half_rise_min"=>half.left, "half_fall_min"=>half.right, "fwhm_min"=>half.width,
        "fwhm_reading_bounds_min"=>[inner.width, outer.width],
        "visible_mean_min"=>visible_mean(t, y, background),
        "visible_mean_reading_scenarios_min"=>collect(extrema(means)),
        "common_window_min"=>[start, stop], "common_window_visible_mean_min"=>visible_mean(common, common_y, background))
end

"""Discrete linear convolution of masses. Index 1 corresponds to zero transit time."""
function convolve(a, b)
    result = zeros(length(a) + length(b) - 1)
    for (i, mass) in pairs(a)
        @inbounds @simd for j in eachindex(b)
            result[i+j-1] += mass * b[j]
        end
    end
    return result
end

function kernel(grid, layer_fraction, dispersion)
    # Dimensionless drift speed is 1 and dispersion is alpha/L, where L is the full length.
    dz = step(grid)
    mass = zeros(length(grid))
    for i in 2:length(grid)
        t = grid[i]
        mass[i] = dz * layer_fraction / sqrt(4π * dispersion * t^3) *
                  exp(-(layer_fraction - t)^2 / (4 * dispersion * t))
    end
    return mass
end

function peak_time(grid, values)
    i = argmax(values)
    left, peak, right = values[i-1:i+1]
    offset = (left-right) / (2 * (left - 2peak + right))
    return grid[i] + offset * step(grid)
end

"""Series reference with a rectangular inlet pulse and two independent layer transit times."""
function baseline(column, dispersivity, resolution, horizon, pulse_fraction)
    dz = pulse_fraction / resolution
    grid = range(0; step=dz, length=ceil(Int, horizon / dz)+1)
    length_mm = column["length_mm"]
    fine = kernel(grid, column["fine_length_mm"]/length_mm, dispersivity["fine_dispersivity_mm"]/length_mm)
    coarse = kernel(grid, column["coarse_length_mm"]/length_mm, dispersivity["coarse_dispersivity_mm"]/length_mm)
    inlet = fill(1.0/resolution, resolution+1)
    inlet[[1,end]] ./= 2
    fc = convolve(convolve(inlet, fine), coarse)
    cf = convolve(convolve(inlet, coarse), fine)
    time = range(0; step=dz, length=length(fc))
    total = sum(fc)
    center = sum(time .* fc) / total
    variance = sum((time .- center).^2 .* fc) / total
    # A collection interval equal to the pulse width, labelled at the interval midpoint.
    averaged = convolve(fc, inlet)
    bin_time = range(-pulse_fraction/2; step=dz, length=length(averaged))
    analytical_variance = 2 * (column["fine_length_mm"]*dispersivity["fine_dispersivity_mm"] +
        column["coarse_length_mm"]*dispersivity["coarse_dispersivity_mm"]) / length_mm^2 + pulse_fraction^2/12
    record = Dict("pulse_steps"=>resolution, "mass"=>total, "mean_pv"=>center, "variance_pv2"=>variance,
        "analytical_mean_pv"=>1+pulse_fraction/2, "analytical_variance_pv2"=>analytical_variance,
        "max_directional_density_difference_per_pv"=>maximum(abs.(fc - cf))/dz,
        "peak_pv"=>peak_time(time, fc), "fwhm_pv"=>half_height(time, fc, 0.0).width,
        "collection_center_peak_pv"=>peak_time(bin_time, averaged),
        "collection_center_fwhm_pv"=>half_height(bin_time, averaged, 0.0).width)
    return record
end

function evaluate()
    source = TOML.parsefile(joinpath(HERE, "evidence.toml"))
    figure = TOML.parsefile(joinpath(HERE, "fig2_markers.toml"))
    config, column, pulse = source["screen"], source["column"], source["pulse"]
    single = source["single_packing"]
    axes = (time=axis_map(figure["time_tick_x_pt"], figure["time_tick_min"]),
            ec=axis_map(figure["ec_tick_y_pt"], figure["ec_tick_uS"]))
    traces = [read_series(s, axes, config["same_time_group_tolerance_pt"]) for s in figure["series"]]
    area_mm2 = π * (column["inner_diameter_mm"]/2)^2
    volume_ml = area_mm2 * column["length_mm"] * column["porosity"] / 1000
    pulse_fraction = pulse["volume_ml"]/volume_ml
    grids = [baseline(column, single, n, config["single_segment_horizon_pore_volumes"], pulse_fraction)
             for n in config["pulse_steps"]]
    finest = last(grids)
    curves = Dict[]
    rates = Dict[]
    for (i, flow) in enumerate(pulse["flow_ml_min"])
        per_rate = filter(t -> t.flow == flow, traces)
        common = (maximum(first(t.times) for t in per_rate), minimum(last(t.times) for t in per_rate))
        rows = [screen_series(t, config, common) for t in per_rate]
        append!(curves, rows)
        transit = volume_ml / flow
        speed = flow * 1000 / (area_mm2 * column["porosity"])
        row = Dict{String,Any}("flow_ml_min"=>flow, "pore_speed_mm_min"=>speed, "bed_mean_min"=>transit,
            "pulse_centroid_min"=>transit+pulse["pulse_min"][i]/2,
            "baseline_peak_min"=>finest["peak_pv"]*transit, "baseline_fwhm_min"=>finest["fwhm_pv"]*transit,
            "center_labelled_collection_peak_min"=>finest["collection_center_peak_pv"]*transit,
            "center_labelled_collection_fwhm_min"=>finest["collection_center_fwhm_pv"]*transit,
            "fine_D_ax_mm2_min"=>speed*single["fine_dispersivity_mm"],
            "coarse_D_ax_mm2_min"=>speed*single["coarse_dispersivity_mm"],
            "dispersivity_basis"=>(flow == single["flow_ml_min"] ? "reported single-packing fit" : "fixed-alpha transfer scenario"))
        for direction in ("FC", "CF")
            runs = filter(t -> t["direction"] == direction, rows)
            row[direction*"_peak_min"] = [minimum(t["possible_peak_time_min"][1] for t in runs),
                                          maximum(t["possible_peak_time_min"][2] for t in runs)]
            row[direction*"_fwhm_min"] = collect(extrema(t["fwhm_min"] for t in runs))
            row[direction*"_fwhm_bounds_min"] = [minimum(t["fwhm_reading_bounds_min"][1] for t in runs),
                                                 maximum(t["fwhm_reading_bounds_min"][2] for t in runs)]
            row[direction*"_fwhm_pumped_ml"] = flow .* row[direction*"_fwhm_min"]
            row[direction*"_visible_mean_min"] = collect(extrema(t["visible_mean_min"] for t in runs))
        end
        row["CF_minus_FC_peak_min"] = [row["CF_peak_min"][1]-row["FC_peak_min"][2],
                                       row["CF_peak_min"][2]-row["FC_peak_min"][1]]
        row["CF_minus_FC_fwhm_pumped_ml"] = [row["CF_fwhm_pumped_ml"][1]-row["FC_fwhm_pumped_ml"][2],
                                              row["CF_fwhm_pumped_ml"][2]-row["FC_fwhm_pumped_ml"][1]]
        row["CF_minus_FC_fwhm_bounds_min"] = [row["CF_fwhm_bounds_min"][1]-row["FC_fwhm_bounds_min"][2],
                                              row["CF_fwhm_bounds_min"][2]-row["FC_fwhm_bounds_min"][1]]
        push!(rates, row)
    end
    inputs = ["evidence.toml", "fig2_markers.toml", "extract_fig2.py", "audit.jl"]
    return Dict("source_doi"=>source["doi"], "pore_volume_ml"=>volume_ml,
        "injected_NaCl_mg"=>pulse["feed_NaCl_mg_L"]*pulse["volume_ml"]/1000,
        "column_to_bead_diameter_ratios"=>[column["inner_diameter_mm"]/column[key]
            for key in ("fine_bead_diameter_mm","coarse_bead_diameter_mm")],
        "time_axis_tick_error_min"=>axes.time.max_tick_error, "ec_axis_tick_error_uS"=>axes.ec.max_tick_error,
        "code_and_input_sha256"=>Dict(name=>bytes2hex(sha256(read(joinpath(HERE, name)))) for name in inputs),
        "curves"=>curves, "rates"=>rates, "baseline_convergence"=>grids)
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 1 || error("Pass the output TOML path")
    result = BerkowitzAudit.evaluate()
    open(ARGS[1], "w") do out
        BerkowitzAudit.TOML.print(out, result; sorted=true)
    end
end
