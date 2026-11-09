"""
    filtration_load(k_support_m2, support_thickness_m, selective_resistance_m_inv,
                    viscosity_Pa_s, flux_m_s, membrane_area_m2, efficiency;
                    fixed_resistance_m_inv=0.0)

Clean, saturated, isothermal, incompressible layers in series at common normal
flux. Resistances R=L/k have units m⁻¹; Δp=μJR. The selective resistance is an
explicit supplied property, never inferred from the support permeability.
Efficiency is shaft-to-hydraulic efficiency: motor efficiency is not included.
No crossflow circulation, fouling, osmotic pressure or backwash is modelled.
"""
function filtration_load(k_support_m2, support_thickness_m, selective_resistance_m_inv,
                         viscosity_Pa_s, flux_m_s, membrane_area_m2, efficiency;
                         fixed_resistance_m_inv=0.0)
    k = positive("support permeability_m2", k_support_m2)
    L = positive("support thickness_m", support_thickness_m)
    mu = positive("viscosity_Pa_s", viscosity_Pa_s)
    J = finite("forward filtration flux_m_s", flux_m_s)
    A = positive("membrane area_m2", membrane_area_m2)
    eta = positive("pump efficiency fraction", efficiency)
    J >= 0 || throw(ArgumentError("Forward filtration flux must be >= 0; backwash is outside this model."))
    eta <= 1 || throw(ArgumentError("Pump efficiency fraction must be <= 1."))
    Rl = finite("selective resistance_m_inv", selective_resistance_m_inv)
    Rf = finite("fixed resistance_m_inv", fixed_resistance_m_inv)
    min(Rl,Rf) >= 0 || throw(ArgumentError("Selective and other fixed resistances must be >= 0."))
    Rs = positive("support resistance_m_inv", L/k)
    Rt = positive("total resistance_m_inv", Rs+Rl+Rf)
    Q = finite("filtrate flow_m3_s", J*A)
    dp_s, dp_l, dp_f = mu*J .* (Rs,Rl,Rf)
    dp = finite("total pressure drop_Pa", mu*J*Rt)
    power = finite("shaft power_W", Q*dp/eta)
    # At J=0 the limiting energy per volume is zero in this purely linear model.
    (; support_resistance_m_inv=Rs, selective_resistance_m_inv=Rl,
       fixed_resistance_m_inv=Rf, total_resistance_m_inv=Rt,
       support_pressure_drop_Pa=dp_s, selective_pressure_drop_Pa=dp_l,
       fixed_pressure_drop_Pa=dp_f, total_pressure_drop_Pa=dp,
       flow_m3_s=Q, shaft_power_W=power,
       specific_shaft_energy_kWh_m3=finite("specific shaft energy_kWh_m3",dp/(eta*3.6e6)),
       selective_resistance_fraction=Rl/Rt, support_resistance_fraction=Rs/Rt,
       selective_dominates=Rl>Rs+Rf, selective_at_least_90_percent=Rl>=9*(Rs+Rf))
end

"""Convert pure-water permeance in L m⁻² h⁻¹ bar⁻¹ to m Pa⁻¹ s⁻¹."""
permeance_LMH_bar_to_SI(x) = positive("pure water permeance_LMH_bar",x)*1e-3/(3600*1e5)

"""Resistance from permeance of the SAME specified layer/specimen at stated μ."""
resistance_from_permeance(permeance_m_Pa_s, viscosity_Pa_s) =
    positive("resistance_m_inv",1/(positive("viscosity_Pa_s",viscosity_Pa_s)*positive("permeance_m_Pa_s",permeance_m_Pa_s)))
