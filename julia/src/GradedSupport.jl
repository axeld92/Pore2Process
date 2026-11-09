"""Independent lubrication reference for piecewise-linear circular radii.

The integral of r^-4 is exact for the declared profile, but reducing 3D Stokes
flow to locally developed sections is approximate for tapered pores. Only a
constant radius supplies an exact Hagen–Poiseuille oracle. Never used by OpenLB.
"""
function graded_pore_reference(x_m::AbstractVector, radii_m::AbstractVector, viscosity_Pa_s, area_m2)
    n=length(x_m)
    n>=2 && length(radii_m)==n || throw(ArgumentError("Require at least two paired axial knots and radii."))
    x=[finite("axial knot_m",v) for v in x_m]
    r=[positive("radius_m",v) for v in radii_m]
    mu,A=positive("viscosity_Pa_s",viscosity_Pa_s),positive("representative area_m2",area_m2)
    all(diff(x).>0) || throw(ArgumentError("Axial knots must be strictly increasing."))
    pi*maximum(r)^2<=A || throw(ArgumentError("Representative area must contain every pore section."))
    integral=0.; volume=0.
    for i in 1:n-1
        L,a,b=x[i+1]-x[i],r[i],r[i+1]
        # Stable form of (a^-3-b^-3)/(3*(b-a)), including equal radii.
        integral+=L*(a*a+a*b+b*b)/(3*a^3*b^3)
        volume+=pi*L*(a*a+a*b+b*b)/3
    end
    L=x[end]-x[1]
    R=positive("lubrication resistance_Pa_s_m3",8mu/pi*integral)
    (; hydraulic_resistance_Pa_s_m3=R,
       permeability_cell_m2=positive("lubrication permeability_m2",mu*L/(A*R)),
       volume_porosity=positive("volume porosity",volume/(A*L)))
end

"""One common baseline calibration carries measured OpenLB resistance ratios
into a literature-informed permeance scenario. No geometry-specific fitting.
Both resistances must use the same body/tap span and representative area.
This scaling is a scenario assumption; it does not reconstruct a specimen.
"""
function graded_support_load(measured_resistance_Pa_s_m3, baseline_resistance_Pa_s_m3,
        baseline_permeance_m_Pa_s, fixed_impedance_Pa_s_m, flux_m_s, area_m2, efficiency)
    R=positive("measured body resistance_Pa_s_m3",measured_resistance_Pa_s_m3)
    R0=positive("baseline body resistance_Pa_s_m3",baseline_resistance_Pa_s_m3)
    P0=positive("baseline permeance_m_Pa_s",baseline_permeance_m_Pa_s)
    ratio=positive("measured body resistance ratio",R/R0)
    P=positive("translated support permeance_m_Pa_s",P0/ratio)
    load=impedance_load(P,fixed_impedance_Pa_s_m,flux_m_s,area_m2,efficiency)
    (; measured_body_resistance_ratio=ratio,support_permeance_m_Pa_s=P,
       support_impedance_Pa_s_m=1/P,non_support_impedance_Pa_s_m=Float64(fixed_impedance_Pa_s_m),load...)
end
