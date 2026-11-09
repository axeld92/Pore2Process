"""A source operating point, SI units. Flow interpretation is a scenario factor.

The source energy denominator is its measured critical flux. Scenario energy
uses the independently prescribed target flux. Area scale-up is parallel
replication at unchanged axial loss, not an identified plant design.
"""
function crossflow_benchmark_point(velocity_m_s, critical_flux_m_s, axial_loss_Pa,
        source_flow_m3_s, source_area_m2, process_area_m2, target_flux_m_s,
        flow_factor, capacity_fraction)
    u=positive("source velocity_m_s",velocity_m_s)
    jc=positive("source critical_flux_m_s",critical_flux_m_s)
    dp=positive("source axial_loss_Pa",axial_loss_Pa)
    qs=positive("source flow_m3_s",source_flow_m3_s)
    sa=positive("source area_m2",source_area_m2)
    a=positive("process area_m2",process_area_m2)
    j=positive("target flux_m_s",target_flux_m_s)
    factor=positive("declared flow factor",flow_factor)
    margin=positive("critical-flux capacity fraction",capacity_fraction)
    margin<=1 || throw(ArgumentError("Critical-flux capacity fraction must be <= 1."))
    qm=positive("scaled module flow_m3_s",qs*factor*(a/sa))
    qp=positive("target permeate flow_m3_s",j*a)
    source_energy=positive("source axial energy_J_m3",dp*qs/(jc*sa))
    energy=positive("scenario axial energy_J_m3",dp*qm/qp)
    (; velocity_m_s=u,critical_flux_m_s=jc,capacity_fraction=margin,
       admissible_flux_m_s=margin*jc,nominal_admissible=j<=margin*jc,
       axial_pressure_loss_Pa=dp,source_flow_m3_s=qs,flow_factor=factor,
       module_reference_flow_m3_s=qm,permeate_reference_flow_m3_s=qp,
       module_to_permeate_ratio=positive("module/permeate ratio",qm/qp),
       channel_resistance_Pa_s_m3=positive("pointwise channel resistance",dp/qm),
       source_specific_axial_energy_J_m3=source_energy,
       scenario_specific_axial_energy_J_m3=energy)
end

"""Monotonic bracketing with independent velocity and flux reading intervals.

Only definitely-left and definitely-right observations bound the query.
Between measurements this assumes a nondecreasing critical flux, not a
linear interpolation or a validated transition law. No extrapolation.
Missing one-sided information is `nothing`, never infinity/NaN.
"""
function critical_flux_bracket(velocity_low,velocity_high,flux_low,flux_high,
        velocity_m_s,target_flux_m_s,capacity_fraction=1.)
    n=length(velocity_low)
    n>=2 && length(velocity_high)==length(flux_low)==length(flux_high)==n ||
        throw(ArgumentError("Require at least two equal-length source interval arrays."))
    lo=[positive("velocity lower bound",v) for v in velocity_low]
    hi=[positive("velocity upper bound",v) for v in velocity_high]
    jl=[positive("flux lower bound",v) for v in flux_low]
    jh=[positive("flux upper bound",v) for v in flux_high]
    all(lo .<= hi) && all(jl .<= jh) || throw(ArgumentError("Source interval bounds are reversed."))
    all(hi[1:end-1] .< lo[2:end]) && issorted(jl) && issorted(jh) ||
        throw(ArgumentError("Require ordered disjoint velocity intervals and nondecreasing flux bounds."))
    u=positive("query velocity_m_s",velocity_m_s)
    j=positive("target flux_m_s",target_flux_m_s)
    margin=positive("capacity fraction",capacity_fraction)
    margin<=1 || throw(ArgumentError("Capacity fraction must be <= 1."))
    lo[1]<=u<=hi[end] || throw(ArgumentError("Velocity is outside the measured reading domain; extrapolation is unsupported."))
    left=findlast(v->v<=u,hi)
    right=findfirst(v->v>=u,lo)
    lower=isnothing(left) ? nothing : margin*jl[left]
    upper=isnothing(right) ? nothing : margin*jh[right]
    status=!isnothing(lower) && j<=lower ? "admissible" :
        !isnothing(upper) && j>upper ? "inadmissible" : "indeterminate"
    (; lower_flux_m_s=lower,upper_flux_m_s=upper,classification=status,
       lower_source_index=left,upper_source_index=right)
end

"""Passive-solute bookkeeping for the existing lumped water recycle topology.

R = 1-C_permeate/C_module_inlet is an assumed single-pass rejection.
The prescribed module-inlet concentration matches a reference feed series;
fresh concentration is adjusted with circulation to maintain that condition.
This does not change DWSIM's pure-water properties or simulate particle physics.
"""
function crossflow_solute_balance(module_concentration_kg_m3,permeate_flow_m3_s,
        recovery,module_to_permeate_ratio,rejection)
    cm=finite("module concentration_kg_m3",module_concentration_kg_m3)
    cm>=0 || throw(ArgumentError("Concentration must be nonnegative."))
    qp=positive("permeate flow_m3_s",permeate_flow_m3_s)
    y=positive("recovery",recovery)
    y<1 || throw(ArgumentError("Recovery must be < 1."))
    c=positive("module/permeate ratio",module_to_permeate_ratio)
    c>1/y || throw(ArgumentError("Module/permeate ratio must exceed 1/recovery."))
    r=finite("assumed rejection",rejection)
    0<=r<=1 || throw(ArgumentError("Rejection must be between zero and one."))
    cp=(1-r)*cm
    cr=cm*(c-1+r)/(c-1)
    cf=cm*(1-r*(y*c-1)/(c-1))
    qf=qp/y; qm=qp*c; qr=qm-qp; qreturn=qm-qf; qpurge=qf-qp
    module_residual=finite("module solute residual",qm*cm-qp*cp-qr*cr)
    mixer_residual=finite("mixer solute residual",qf*cf+qreturn*cr-qm*cm)
    overall_residual=finite("overall solute residual",qf*cf-qp*cp-qpurge*cr)
    (; module_inlet_concentration_kg_m3=cm,fresh_concentration_kg_m3=cf,
       permeate_concentration_kg_m3=cp,retentate_and_purge_concentration_kg_m3=cr,
       assumed_rejection=r,module_solute_residual_kg_s=module_residual,
       mixer_solute_residual_kg_s=mixer_residual,overall_solute_residual_kg_s=overall_residual)
end
