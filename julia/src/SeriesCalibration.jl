"""Fit a single effective non-support impedance Z in Pa s/m from paired
permeances P in m/(Pa s). Only the supplied training vectors enter the fit.
Z=1/Ptotal-1/Psupport = μR avoids inventing an experimental viscosity.
The support coefficient is fixed at one, rather than another fitted parameter.
"""
function series_impedance_fit(support_P::AbstractVector, total_P::AbstractVector)
    n=length(support_P)
    n>0 && length(total_P)==n || throw(ArgumentError("Paired permeance vectors must have equal nonzero length."))
    measured_Z=Float64[]
    for (ps,pt) in zip(support_P,total_P)
        ps,pt=positive("support permeance_m_Pa_s",ps),positive("total permeance_m_Pa_s",pt)
        pt<ps || throw(ArgumentError("Whole-membrane permeance must be smaller than its paired support permeance."))
        push!(measured_Z,positive("non-support impedance_Pa_s_m",1/pt-1/ps))
    end
    Z=positive("fitted non-support impedance_Pa_s_m",sum(measured_Z)/n)
    (; fixed_impedance_Pa_s_m=Z,training_count=n,
       limiting_permeance_m_Pa_s=positive("limiting permeance_m_Pa_s",1/Z),
       impedance_rms_residual_Pa_s_m=finite("impedance RMS residual",norm(measured_Z .- Z)/sqrt(n)))
end

"""Series permeance/flux and resistance shares at a finite signed Δp.
P is m/(Pa s), Z is Pa s/m, Δp is Pa. No water viscosity is assumed.
The effective non-support share need not be intrinsic selective-layer share.
"""
function series_flux(support_P, fixed_Z, delta_p_Pa)
    P=positive("support permeance_m_Pa_s",support_P)
    Z=finite("fixed impedance_Pa_s_m",fixed_Z)
    Z>=0 || throw(ArgumentError("Fixed impedance must be nonnegative."))
    dp=finite("delta_p_Pa",delta_p_Pa)
    Zs=positive("support impedance_Pa_s_m",1/P)
    Zt=positive("total impedance_Pa_s_m",Zs+Z)
    (; support_impedance_Pa_s_m=Zs,total_impedance_Pa_s_m=Zt,
       total_permeance_m_Pa_s=positive("total permeance_m_Pa_s",1/Zt),
       flux_m_s=finite("flux_m_s",dp/Zt),non_support_fraction=Z/Zt,
       support_fraction=Zs/Zt)
end

"""Fixed-flux clean-water pressure/shaft-energy demonstration using measured
hydraulic impedances; same fluid conditions and linear Darcy flow assumed.
"""
function impedance_load(support_P, fixed_Z, flux_m_s, area_m2, efficiency)
    s=series_flux(support_P,fixed_Z,0.)
    J=finite("forward flux_m_s",flux_m_s)
    J>=0 || throw(ArgumentError("Forward filtration flux must be nonnegative."))
    A=positive("membrane area_m2",area_m2)
    eta=positive("pump efficiency",efficiency)
    eta<=1 || throw(ArgumentError("Pump efficiency must be <=1."))
    Q=finite("flow_m3_s",J*A)
    dp=finite("pressure_Pa",J*s.total_impedance_Pa_s_m)
    (; pressure_drop_Pa=dp,flow_m3_s=Q,shaft_power_W=finite("shaft power_W",Q*dp/eta),
       specific_shaft_energy_kWh_m3=finite("specific shaft energy_kWh_m3",dp/(eta*3.6e6)),
       non_support_fraction=s.non_support_fraction,support_fraction=s.support_fraction)
end
