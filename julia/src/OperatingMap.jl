"""A bounded literature critical-flux comparison, with SI flux inputs.

Linear interpolation is an explicit scenario; monotonic bounds use only source
points on either side. Neither establishes sustainable operation on a synthetic
membrane. Reading intervals and transfer factors are separate from uncertainty.
"""
function critical_flux_capacity(re_values, flux_nominal, flux_lower, flux_upper,
        query_Re, maximum_flux_m_s, capacity_fraction, transfer_factor)
    n=length(re_values)
    n>=2 && length(flux_nominal)==length(flux_lower)==length(flux_upper)==n ||
        throw(ArgumentError("Require equal source arrays with at least two points."))
    re=positive.("source Reynolds number", re_values)
    f=positive.("source critical flux_m_s", flux_nominal)
    lo=positive.("source lower reading_m_s", flux_lower)
    hi=positive.("source upper reading_m_s", flux_upper)
    all(diff(re).>0) && all(diff(f).>=0) && all(diff(lo).>=0) && all(diff(hi).>=0) ||
        throw(ArgumentError("Source Reynolds numbers must increase and fluxes be nondecreasing."))
    all(lo .<= f) && all(f .<= hi) || throw(ArgumentError("Nominal flux must lie within its reading interval."))
    q=positive("query Reynolds number",query_Re); peak=positive("maximum local flux_m_s",maximum_flux_m_s)
    margin=positive("capacity fraction",capacity_fraction); margin<=1 || throw(ArgumentError("Capacity fraction must be <= 1."))
    transfer=positive("declared transfer factor",transfer_factor)
    re[1]<=q<=re[end] || throw(ArgumentError("Outside source Reynolds range; extrapolation is unsupported."))
    left=searchsortedlast(re,q); right=searchsortedfirst(re,q)
    scale=margin*transfer
    nominal = left==right ? f[left] : f[left]+(f[right]-f[left])*(q-re[left])/(re[right]-re[left])
    lowlinear = left==right ? lo[left] : lo[left]+(lo[right]-lo[left])*(q-re[left])/(re[right]-re[left])
    highlinear = left==right ? hi[left] : hi[left]+(hi[right]-hi[left])*(q-re[left])/(re[right]-re[left])
    lower=scale*lo[left]; upper=scale*hi[right]
    status=peak<=lower ? "admissible" : peak>upper ? "inadmissible" : "indeterminate"
    (; query_Re=q, nominal_capacity_m_s=scale*nominal,
       linear_lower_reading_capacity_m_s=scale*lowlinear, linear_upper_reading_capacity_m_s=scale*highlinear,
       monotonic_lower_capacity_m_s=lower, monotonic_upper_capacity_m_s=upper,
       nominal_admissible=peak<=scale*nominal,
       linear_reading_admissible=peak<=scale*lowlinear,
       monotonic_classification=status, lower_source_index=left,upper_source_index=right)
end

"""Local-profile bound on omitted axial acceleration relative to duct friction.

For q'=-P*J and momentum factor beta, the ratio is 2*beta*rho*C*P*J/(mu*A^2).
The supplied beta bound umax/umean follows from integral(u^2)<=umax*integral(u)
for nonnegative axial flow. This is a reduced-model screen, not CFD validation.
"""
function axial_momentum_screen(C_m4,area_m2,perimeter_m,viscosity_Pa_s,density_kg_m3,
        maximum_flux_m_s,maximum_to_mean_velocity)
    C=positive("geometry coefficient_m4",C_m4);A=positive("open area_m2",area_m2)
    P=positive("perimeter_m",perimeter_m);mu=positive("viscosity_Pa_s",viscosity_Pa_s)
    rho=positive("density_kg_m3",density_kg_m3);j=positive("maximum local flux_m_s",maximum_flux_m_s)
    beta=positive("maximum/mean axial velocity",maximum_to_mean_velocity)
    beta>=1 || throw(ArgumentError("Maximum/mean velocity must be >= 1."))
    finite("momentum/friction screen",2beta*rho*C*P*j/(mu*A^2))
end
