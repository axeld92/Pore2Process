module ColumnContact
using ..FilmParticle
export solve_contact, tracer_step

"""
N equal, perfectly mixed contact sections in series. N is a physical model
choice, NOT a spatial mesh: changing N changes the residence-time distribution.
Each section holds eps*A*L/N liquid and (1-eps)*rho_particle*A*L/N dry solid.
Reuse the frozen, conservative FilmParticle coupling with zero inter-section
backmixing. No axial-dispersion coefficient or inactive mass is introduced.
"""
function solve_contact(;sections=1,length_m=.1,area_m2=pi*.011^2/4,porosity=.427,
        superficial_velocity_m_s=43.8e-6/60/(pi*.011^2/4),
        inlet_concentration_kg_m3=.15,initial_concentration_kg_m3=0.,kwargs...)
    n=FilmParticle.count("sections",sections,1,128)
    L=FilmParticle.positive("length_m",length_m);A=FilmParticle.positive("area_m2",area_m2)
    eps=FilmParticle.positive("porosity",porosity)
    eps<1 || throw(ArgumentError("porosity must be < 1."))
    u=FilmParticle.positive("superficial_velocity_m_s",superficial_velocity_m_s;zero=true)
    cin=FilmParticle.positive("inlet_concentration_kg_m3",inlet_concentration_kg_m3;zero=true)
    c0=FilmParticle.positive("initial_concentration_kg_m3",initial_concentration_kg_m3;zero=true)
    p=FilmParticle.parameters(;kwargs...)
    mass=(1-eps)*p.rho*A*L
    r=FilmParticle.integrate_impl(u;nz=n,L,A,storage=eps,bulk_density=(1-eps)*p.rho,
        physical_mass=mass,active_mass=mass,cin,cinit=c0,dispersion=0.,p...)
    (;r...,sections=n,liquid_residence_time_s=u==0 ? nothing : eps*L/u,
      section_liquid_volume_m3=eps*A*L/n,section_dry_mass_kg=mass/n)
end

"""Independent analytical nonadsorbing step response, no particle integrator.
For equal sections the RTD is Erlang: F=1-exp(-x)*sum(x^k/k!,k=0:N-1),
x=N*t/tau. Evaluate positive Poisson upper-tail terms at small x to avoid
catastrophic cancellation. Intended for verification, not fitting adsorption.
"""
function tracer_step(time_s;sections=1,residence_time_s=1.)
    n=FilmParticle.count("sections",sections,1,128)
    tau=FilmParticle.positive("residence_time_s",residence_time_s)
    t=FilmParticle.positive("time_s",time_s;zero=true)
    x=n*t/tau
    x==0 && return 0.
    x>1000 && return 1.
    if x<n
        term=exp(-x)
        for k in 1:n;term*=x/k;end
        result=term
        for k in n+1:10000
            term*=x/k;result+=term
            term<=eps(Float64)*max(result,1e-300) && break
        end
        return result
    end
    term=exp(-x);survival=term
    for k in 1:n-1;term*=x/k;survival+=term;end
    1-survival
end
end
