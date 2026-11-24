module InventoryAudit

using Statistics
using TOML

const REVIEW = @__DIR__
const SOURCES = normpath(joinpath(REVIEW, "..", "..", "sources"))
const protocol = TOML.parsefile(joinpath(REVIEW, "protocol.toml"))
source(id) = TOML.parsefile(joinpath(SOURCES, id, "evidence.toml"))

"""Linear reading of an axis with upper-left image pixel origin."""
axis_value(px, p_lo, p_hi, lo, hi) =
    lo + (px - p_lo) * (hi - lo) / (p_hi - p_lo)

function interp(x, y, t)
    length(x) == length(y) >= 2 && all(diff(x) .> 0) ||
        throw(ArgumentError("Interpolation requires a strictly ordered axis"))
    first(x) <= t <= last(x) || throw(DomainError(t, "Outside reference coverage"))
    j = clamp(searchsortedlast(x, t), 1, length(x)-1)
    y[j] + (t-x[j])/(x[j+1]-x[j])*(y[j+1]-y[j])
end

rmse(r) = sqrt(mean(abs2, r))
range_width(x) = maximum(x) - minimum(x)

function layer_fill(mass_g, volume_ml, density_g_ml, factor)
    all(isfinite, (mass_g, volume_ml, density_g_ml, factor)) ||
        throw(ArgumentError("Non-finite layer inputs"))
    mass_g >= 0 && volume_ml > 0 && density_g_ml > 0 && factor > 0 ||
        throw(ArgumentError("Invalid layer inputs"))
    dry = mass_g / (density_g_ml * volume_ml)
    dry, factor * dry
end

function damstrup_audit(data, p, target)
    f = data["figure2"]
    curves = Dict()
    for s in f["series"]
        t = axis_value.(s["x_px"], f["x0_px"], f["x40_px"], 0., 40.)
        x = axis_value.(s["y_px"], f["y0_px"], f["y60_px"], 0., 60.)
        curves[s["mass_g"]] = (t=t, x=x)
    end
    ref = curves[p["reference_charge_g"]]
    tol_t = f["reading_time_min"]
    tol_x = f["reading_mag_wt_percentage_points"]
    density = data["column"]["dry_bulk_density_g_ml"]
    expan = data["table2"]["expansion_percent_of_dry"] / 100
    expan_sd = data["table2"]["expansion_sd_pp"] / 100
    fill(m, v) = begin
        d, w = layer_fill(m, v, density, expan)
        Dict("dry_fraction" => d, "free_wet_layer_fraction" => w,
             "free_wet_fraction_at_expansion_plus_minus_reported_sd" =>
                 [d*(expan-expan_sd), d*(expan+expan_sd)],
             "expansion_factor_to_fill" => 1/d)
    end
    fills = [merge(Dict("mass_g" => s["mass_g"],
                        "printed_dry_fill_percent"=>s["printed_dry_fill_percent"]),
                   fill(s["mass_g"], data["column"]["volume_ml"])) for s in f["series"]]
    holdouts = []
    for m in p["holdout_charges_g"]
        points = []
        for (t, obs) in zip(curves[m].t, curves[m].x)
            # Allow near-coincident boundary markers within two reading
            # allowances. Do not extrapolate to the unmatched 5 min points.
            first(ref.t)-2tol_t <= t <= last(ref.t)+2tol_t || continue
            tc = clamp(t, first(ref.t), last(ref.t))
            pred = interp(ref.t, ref.x, tc)
            tl = clamp(t-2tol_t, first(ref.t), last(ref.t))
            th = clamp(t+2tol_t, first(ref.t), last(ref.t))
            vals = [interp(ref.t, ref.x, a) for a in
                    vcat(tl, filter(z -> tl < z < th, ref.t), th)]
            lo, hi = obs - maximum(vals) - 2tol_x, obs - minimum(vals) + 2tol_x
            push!(points, Dict("time_min"=>t, "reference_time_min"=>tc,
                              "observed_mag_wt_percent"=>obs,
                              "reference_mag_wt_percent"=>pred,
                              "difference_pp"=>obs-pred,
                              "reading_difference_range_pp"=>[lo, hi],
                              "differs_beyond_reading_only"=>lo>0 || hi<0))
        end
        r = [v["difference_pp"] for v in points]
        push!(holdouts, Dict("mass_g"=>m, "n"=>length(points),
                            "rmse_pp"=>rmse(r), "max_abs_pp"=>maximum(abs, r),
                            "outside_reading_envelope_count"=>
                                count(q->q["differs_beyond_reading_only"],points),
                            "points"=>points))
    end
    matching = filter(r->r[4]==p["scaleup_reference_loading_kg_m3"] &&
                         r[5]==p["scaleup_nominal_time_min"], data["table1"]["rows"])
    base = only(filter(r->r[1]==p["scaleup_reference_diameter_mm"], matching))
    mass(r) = r[3]*r[4]/1000 # mL * kg m^-3 -> g
    scaleup = [Dict("diameter_mm"=>r[1], "length_mm"=>r[2], "mass_g"=>mass(r),
                    "Q_per_mass_ml_min_g"=>r[6]/mass(r),
                    "mag_wt_percent"=>r[7], "mag_sd_pp"=>r[8],
                    "relative_mass"=>mass(r)/mass(base),
                    "relative_Q_times_MAG_proxy"=>r[6]*r[7]/(base[6]*base[7]),
                    "difference_from_reference_pp"=>r[7]-base[7]) for r in matching]
    target_volume = pi*(target["diameter_mm"]/2)^2*target["length_mm"]/1000
    return Dict("wetting"=>fills, "figure2_holdouts"=>holdouts,
                "table1_scaleup"=>scaleup,
                "target_analogue_scenario"=>
                    merge(Dict("volume_ml"=>target_volume, "mass_g"=>target["mass_g"]),
                          fill(target["mass_g"],target_volume)))
end

activity(k, f, t) = f + (1-f)*exp(-k*t)

"""Profile a physical retained fraction analytically at each positive k."""
function profile_activity(k, t, a, allow_retained)
    e = exp.(-k .* t)
    b = 1 .- e
    denom = sum(abs2, b)
    f = allow_retained && denom>0 ? clamp(sum(b .* (a .- e))/denom, 0., 1.) : 0.
    (k=k, retained=f, sse=sum(abs2, activity.(k, f, t) .- a))
end

function fit_activity(t, a, allow_retained; points=600, kmin=1e-6, kmax=10.)
    length(t)==length(a)>=2 && all(isfinite,t) && all(t.>=0) &&
        all(isfinite,a) && all(a.>=0) && 0<kmin<kmax && points>=3 ||
        throw(ArgumentError("Invalid activity fit"))
    logs=range(log(kmin),log(kmax);length=points)
    objective(z)=profile_activity(exp(z),t,a,allow_retained)
    scores=objective.(logs)
    i=argmin([v.sse for v in scores])
    # Include the boundary candidates: no loss and the two grid edges.
    candidates=[scores[i],profile_activity(0.,t,a,allow_retained)]
    lo, hi = logs[max(1,i-1)], logs[min(points,i+1)]
    # Refine within the bracket in log(k) with a deterministic golden search.
    golden=(sqrt(5.)-1)/2
    for _ in 1:90
        l, r = hi-golden*(hi-lo), lo+golden*(hi-lo)
        if objective(l).sse <= objective(r).sse
            hi=r
        else
            lo=l
        end
    end
    push!(candidates,objective((lo+hi)/2))
    candidates[argmin([v.sse for v in candidates])]
end

"""Least-squares natural-log slope with a(0)=1, using only the declared early data."""
early_log_coefficient(t, a) = -sum(t .* log.(a)) / sum(abs2, t)

function rejasse_audit(data, p)
    f = data["figure3"]
    t = Float64.(f["times_min"])
    a0 = f["fresh_activity_umol_min_mg"]
    early = t .<= p["early_log_training_end_min"]
    train = t .<= p["activity_shape_training_end_min"]
    test = only(findall(==(p["late_holdout_min"]),t))
    @assert !train[test]
    pts=p["log_rate_grid_points"]
    kmin=p["log_rate_grid_min_per_min"]
    kmax=p["log_rate_grid_max_per_min"]
    fit(t,a,ret; points=pts) = fit_activity(t,a,ret;points,kmin,kmax)
    series = []
    for s in f["conventional"]
        u = vcat(a0, axis_value.(s["y_positive_time_px"], f["y0_px"], f["y6_px"], 0., 6.))
        a = u ./ a0
        early_k = early_log_coefficient(t[early],a[early])
        # Table 1 is the slope of log10(percentage) vs time, not ln(a).
        table = only(filter(row->row[1]==s["substrate"] && row[2]=="0" &&
                                       row[3]=="conventional",data["table1"]["rows"]))
        author_ln_per_min = log(10)*table[4]/60
        fits = []
        for ret in (false,true)
            r = fit(t[train],a[train],ret)
            fine = fit(t[train],a[train],ret;points=2pts)
            late = 100activity(r.k,r.retained,t[test])
            margin=f["reading_activity_umol_min_mg"]/a0
            # Only perturb the training marker readings, preserving the
            # printed fresh activity. This is a reading scenario, not a CI.
            train_a=a[train]
            train_t=t[train]
            n=length(train_a)-1
            read_preds=Float64[]
            for mask in 0:(2^n-1)
                pert=copy(train_a)
                for j in 1:n
                    pert[j+1]=max(0.,pert[j+1]+(((mask>>(j-1))&1)==0 ? -margin : margin))
                end
                z=fit(train_t,pert,ret)
                push!(read_preds,100activity(z.k,z.retained,t[test]))
            end
            push!(fits,Dict("shape"=>ret ? "fast_plus_retained" : "single_exponential",
                            "k_natural_per_min"=>r.k,"retained_fraction"=>r.retained,
                            "train_rmse_pp"=>100sqrt(r.sse/count(train)),
                            "heldout_30min_percent"=>late,
                            "heldout_error_vs_digitized_pp"=>late-100a[test],
                            "heldout_error_vs_printed_pp"=>late-s["reported_30min_remaining_percent"],
                            "reading_training_scenario_min_max_pp"=>[minimum(read_preds),maximum(read_preds)],
                            "grid_refinement_holdout_change_pp"=>
                                abs(late-100activity(fine.k,fine.retained,t[test]))))
        end
        push!(series,Dict("substrate"=>s["substrate"],"time_min"=>t,
                          "observed_activity_percent"=>100 .* a,
                          "author_30min_remaining_percent"=>s["reported_30min_remaining_percent"],
                          "early_log_k_natural_per_min"=>early_k,
                          "author_k_printed_per_h"=>table[4],
                          "author_k_natural_per_min"=>author_ln_per_min,
                          "author_half_life_h"=>table[5],
                          "half_life_using_ln10_in_h"=>log(2)/author_ln_per_min/60,
                          "early_log_prediction_10_30_percent"=>100 .* exp.(-early_k .* [10.,30.]),
                          "shape_fits"=>fits))
    end
    g=data["figure5"]
    rates=axis_value.(g["y_conventional_px"],g["y0_px"],g["y100_px"],0.,100.)
    Dict("storage_series"=>series,
         "reacting_reuse_initial_rate_umol_min"=>rates,
         "reacting_reuse_relative_initial_rate_percent"=>100 .* rates ./ rates[1])
end

function lerin_audit(data, p)
    out=[]
    for key in ("figure1","figure2","figure3")
        f=data[key]
        convert(y)=axis_value.(y,f["y_low_px"],f["y_high_px"],f["low_percent"],f["high_percent"])
        c=convert(f["conversion_y_px"])
        a=convert(f["activity_y_px"])
        changes=diff(a)
        r=p["comparison_reference_cycle"]
        final=last(f["cycle"])
        i=only(findall(==(r),f["cycle"]))
        # Compare changes directly, not a continuous-age deactivation fit.
        push!(out,Dict("figure"=>key,"product"=>f["product"],"cycle"=>f["cycle"],
                       "conversion_percent"=>c,"activity_percent"=>a,
                       "activity_increase_destinations_beyond_reading"=>
                           [f["cycle"][j+1] for j in eachindex(changes) if changes[j]>2*f["reading_pp"]],
                       "reference_comparison_cycle"=>r,"final_plotted_cycle"=>final,
                       "conversion_change_from_reference_pp"=>c[end]-c[i],
                       "activity_change_from_reference_pp"=>a[end]-a[i],
                       "final_to_reference_conversion_ratio"=>c[end]/c[i],
                       "final_to_reference_activity_ratio"=>a[end]/a[i]))
    end
    out
end

function strnisa_audit(data, p)
    density=data["catalyst"]["density_assumed_mg_ul"]
    q0=data["batch_assay"]["activity_umol_min_mg_N435"]
    table=Dict()
    for row in data["table2"]["rows"]
        m=density*(row[6]-row[7])
        table[row[1]]=(mass=m,void=row[7],gamma=q0*m/row[7],
                       printed_gamma=row[9],bpn=row[10])
    end
    f=data["figure5"]
    curves=Dict(s["reactor"] => (t=Float64.(s["tau_min"]),
                   x=axis_value.(s["y_px"],s["y0_px"],s["y100_px"],0.,100.))
                for s in f["series"])
    ref_id=p["reference_reactor"]
    ref=curves[ref_id]
    # Zero nominal contact time implies zero conversion; this is an
    # analytical anchor, not another experimental calibration point.
    rx=vcat(0.,ref.t)
    ry=vcat(0.,ref.x)
    reference_density=table[ref_id].mass/table[ref_id].void
    rows=[]
    for id in p["held_out_reactors"]
        c=curves[id]
        r=table[id]
        points=[]
        for (time,obs) in zip(c.t,c.x)
            adjusted=time*(r.mass/r.void)/reference_density
            first(rx)<=adjusted<=last(rx) || throw(DomainError(adjusted,"Exposure outside A reference"))
            unadjusted_prediction=interp(rx,ry,time)
            adjusted_prediction=interp(rx,ry,adjusted)
            push!(points,Dict("tau_min"=>time,"observed_percent"=>obs,
                              "A_at_same_nominal_tau_percent"=>unadjusted_prediction,
                              "A_at_same_inventory_per_flow_percent"=>adjusted_prediction,
                              "inventory_adjusted_tau_min"=>adjusted,
                              "nominal_difference_pp"=>obs-unadjusted_prediction,
                              "inventory_adjusted_difference_pp"=>obs-adjusted_prediction))
        end
        error=[v["inventory_adjusted_difference_pp"] for v in points]
        early=filter(z->z["tau_min"]<=p["early_time_max_min"],points)
        late=filter(z->z["tau_min"]>=p["late_time_min"],points)
        push!(rows,Dict("reactor"=>id,"n"=>length(points),
                        "nominal_rmse_pp"=>rmse([z["nominal_difference_pp"] for z in points]),
                        "inventory_adjusted_rmse_pp"=>rmse(error),
                        "inventory_adjusted_max_abs_pp"=>maximum(abs,error),
                        "early_inventory_adjusted_max_abs_pp"=>maximum(abs,[z["inventory_adjusted_difference_pp"] for z in early]),
                        "late_inventory_adjusted_max_abs_pp"=>maximum(abs,[z["inventory_adjusted_difference_pp"] for z in late]),
                        "points"=>points))
    end
    # Same chemistry/temperature, but mean packed-bed rates at finite
    # conversion are NOT the batch initial rate or a residual-activity test.
    cin_umol_ul=data["table2"]["feed_each_mM"]/1000
    protein_fraction=data["catalyst"]["CALB_mass_fraction_assumed"]
    rates=[]
    for id in vcat([ref_id],p["held_out_reactors"])
        r=table[id]
        c=curves[id]
        results=[]
        for t in p["nominal_comparison_times_min"]
            q=r.void/t
            x=interp(c.t,c.x,t)
            prod=cin_umol_ul*q*x/100 # umol/min
            result=Dict("nominal_tau_min"=>t,"Q_from_void_volume_ul_min"=>q,
                        "conversion_percent"=>x,
                        "mean_butanol_consumption_umol_min_mg_N435"=>prod/r.mass,
                        "stoichiometric_BB_if_all_consumption_forms_product_mmol_h"=>prod*0.06)
            if t==last(c.t)
                # Table 2 BPN is independently based on measured BUTYL
                # BUTYRATE at tau=1.4 min, so retain this production path
                # rather than replacing it by 1:1 butanol consumption.
                bb=r.bpn*(protein_fraction*r.mass/1000)*60/t
                result["reported_BPN_mmol_g_CALB"]=r.bpn
                result["BB_production_from_BPN_mmol_h"]=bb
                result["BB_to_stoichiometric_BB_ratio"]=bb/(prod*0.06)
            end
            push!(results,result)
        end
        push!(rates,Dict("reactor"=>id,"inventory_reconstructed_mg"=>r.mass,
                         "calculated_gamma_U_ul_void"=>r.gamma,
                         "printed_gamma_U_ul_void"=>r.printed_gamma,
                         "cases"=>results))
    end
    rtd=[Dict("reactor"=>row[1],"nominal_volume_over_flow_min"=>
              table[row[1]].void/data["table1"]["flow_ul_min"],
              "reported_mean_min"=>row[7],"reported_variance_min2"=>row[8],
              "reported_skewness"=>row[9]) for row in data["table1"]["rows"]]
    Dict("same_contact_time_holdouts"=>rows,"inventory_and_productivity"=>rates,
         "reported_tracer_moments"=>rtd,
         "nominal_A_channel_clearance_um"=>1000*data["table2"]["rows"][1][4]-
                                                 data["catalyst"]["mean_sieved_diameter_um"])
end

function run_audit()
    Dict("schema_version"=>1, "target_packed_bed_conversion_points_used"=>0,
         "damstrup"=>damstrup_audit(source("damstrup2007"),protocol["damstrup"],
                                    protocol["target_fill_scenario"]),
         "rejasse"=>rejasse_audit(source("rejasse2004"),protocol["rejasse"]),
         "lerin"=>lerin_audit(source("lerin2011"),protocol["lerin"]),
         "strnisa"=>strnisa_audit(source("strnisa2018"),protocol["strnisa"]))
end

end # module

if abspath(PROGRAM_FILE)==@__FILE__
    using TOML
    result=InventoryAudit.run_audit()
    if isempty(ARGS)
        TOML.print(stdout,result;sorted=true)
    elseif length(ARGS)==1
        open(only(ARGS), "w") do io
            TOML.print(io,result;sorted=true)
        end
    else
        error("Usage: julia audit.jl [output.toml]")
    end
end
