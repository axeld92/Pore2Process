module OperatingPointHydraulics
using ..AdsorptionColumn
export translate

"""Transfer a measured native hydraulic operating point at matched Reynolds.
The apparent K is an algebraic representation of this one point, not a
flow-independent permeability or a nonlinear pressure-law extrapolation.
"""
function translate(;native_gradient_Pa_m,native_velocity_m_s,native_viscosity_Pa_s,
        native_density_kg_m3,carrier_viscosity_Pa_s,carrier_density_kg_m3,
        superficial_velocity_m_s,length_m,area_m2,efficiency_fraction=.75)
    positive=AdsorptionColumn.number
    g=positive("native_gradient_Pa_m",native_gradient_Pa_m)
    un=positive("native_velocity_m_s",native_velocity_m_s)
    mun=positive("native_viscosity_Pa_s",native_viscosity_Pa_s)
    rhon=positive("native_density_kg_m3",native_density_kg_m3)
    mu=positive("carrier_viscosity_Pa_s",carrier_viscosity_Pa_s)
    rho=positive("carrier_density_kg_m3",carrier_density_kg_m3)
    u=positive("superficial_velocity_m_s",superficial_velocity_m_s)
    ratio=(mu/rho)/(mun/rhon)
    matched_velocity=un*ratio
    mismatch=abs(u/matched_velocity-1)
    mismatch<=1e-4 || throw(ArgumentError("Native Reynolds point does not match the requested carrier flow; run a new native point."))
    K=mun*un/g
    load=AdsorptionColumn.pressure_load(K,mu,length_m,u,area_m2,efficiency_fraction)
    (;load...,operating_point_apparent_permeability_m2=K,
      native_to_carrier_velocity_scale=ratio,
      native_to_carrier_pressure_scale=(rho/rhon)*ratio^2,
      requested_to_matched_velocity_ratio=u/matched_velocity,
      Reynolds_relative_mismatch=mismatch)
end
end
