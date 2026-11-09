"""
    infer_permeability(length_m, area_m2, viscosity_Pa_s, delta_p_Pa, flow_m3_s)

Infer axial Darcy permeability from measured signed pressure drop and flow in
SI units. The representative area must be supplied explicitly. No geometry
formula is substituted for a missing measurement.
"""
function infer_permeability(length_m, area_m2, viscosity_Pa_s, delta_p_Pa, flow_m3_s)
    L, A, mu = positive("measurement length_m",length_m), positive("representative area_m2",area_m2), positive("viscosity_Pa_s",viscosity_Pa_s)
    dp, Q = finite("measured delta_p_Pa",delta_p_Pa), finite("measured flow_m3_s",flow_m3_s)
    if dp == 0
        Q == 0 || throw(ArgumentError("A nonzero flow with zero measured pressure drop cannot define this passive Darcy permeability."))
        return nothing
    end
    k = finite("inferred permeability_m2",mu*L*Q/(A*dp))
    k >= 0 || throw(ArgumentError("Flow and pressure drop have opposite signs; check the +x orientation convention."))
    k
end

"""Signed Darcy pressure drop for a supplied measured permeability and flow (SI)."""
function predict_column_pressure_drop(permeability_m2, viscosity_Pa_s, length_m, area_m2, flow_m3_s)
    mu, L, A = positive("viscosity_Pa_s",viscosity_Pa_s), positive("column length_m",length_m), positive("column area_m2",area_m2)
    Q = finite("column flow_m3_s",flow_m3_s)
    isnothing(permeability_m2) && return nothing
    k = positive("measured permeability_m2",permeability_m2)
    finite("column pressure drop_Pa",mu*L*Q/(k*A))
end

"""Independent Hagen–Poiseuille oracle for a circular pore inside a declared cell."""
function capillary_reference(radius_m, length_m, viscosity_Pa_s, delta_p_Pa, cell_area_m2)
    r, L, mu, A = positive("radius_m",radius_m), positive("length_m",length_m), positive("viscosity_Pa_s",viscosity_Pa_s), positive("cell_area_m2",cell_area_m2)
    dp = finite("delta_p_Pa",delta_p_Pa)
    pore_area = pi*r^2
    pore_area <= A || throw(ArgumentError("The representative cell area must contain the circular pore area."))
    (; flow_m3_s=finite("reference flow",pi*r^4*dp/(8mu*L)),
       permeability_cell_m2=finite("reference cell permeability",pi*r^4/(8A)),
       permeability_pore_m2=r^2/8, porosity=pore_area/A,
       max_velocity_m_s=finite("reference maximum velocity",r^2*dp/(4mu*L)))
end
