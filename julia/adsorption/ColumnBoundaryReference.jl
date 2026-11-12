module ColumnBoundaryReference
using LinearAlgebra
export dirichlet_neumann_moments, danckwerts_moments, solve_dirichlet_neumann

function scales(length_m, velocity_m_s, dispersion_m2_s)
    all(x->x isa Real && !(x isa Bool) && isfinite(x) && x>0,
        (length_m,velocity_m_s,dispersion_m2_s)) ||
        throw(ArgumentError("Require positive finite length, velocity and dispersion in SI."))
    L,v,D=Float64.((length_m,velocity_m_s,dispersion_m2_s))
    L,v,D,v*L/D
end

"""Uniform ADE: fixed inlet concentration, zero outlet gradient.
Analytical moments from -log[cosh(beta*L)+v/(2D*beta)*sinh(beta*L)],
beta^2=(v/(2D))^2+s/D. This is not the Danckwerts flux inlet.
At small Pe evaluate cancellation-prone expressions at higher precision.
"""
function dirichlet_neumann_moments(;length_m,velocity_m_s,dispersion_m2_s)
    L,v,D,Pe=scales(length_m,velocity_m_s,dispersion_m2_s)
    function evaluate(L,v,D)
        p=v*L/D;E=exp(-p)
        mean=L/v-D/v^2*(1-E)
        variance=2D*L/v^3*(1+2*E)-D^2/v^4*(5-4*E-E^2)
        Float64(mean),Float64(variance)
    end
    m,var=Pe<.1 ? setprecision(128) do
        evaluate(BigFloat(L),BigFloat(v),BigFloat(D))
    end : evaluate(L,v,D)
    (;mean_s=m,variance_s2=var,Peclet=Pe)
end

"""Uniform ADE: imposed total inlet solute flux v*Cin; zero outlet gradient.
Independent closed-vessel RTD moments, with actual liquid volume/Q mean.
"""
function danckwerts_moments(;length_m,velocity_m_s,dispersion_m2_s)
    L,v,D,Pe=scales(length_m,velocity_m_s,dispersion_m2_s)
    variance=2D^2/v^4*(Pe+expm1(-Pe))
    (;mean_s=L/v,variance_s2=variance,Peclet=Pe)
end

"""Independent 1D boundary reference, not a second pore-scale solver.
Central finite volumes/backward Euler; inlet concentration=1, outlet diffusive
flux=0. Inlet total flux includes diffusion. Array values are materialized SI.
"""
function solve_dirichlet_neumann(;length_m,velocity_m_s,dispersion_m2_s,
        t_end_s,dt_s=.001,cells=800)
    L,v,D,Pe=scales(length_m,velocity_m_s,dispersion_m2_s)
    cells isa Integer && !(cells isa Bool) && 4<=cells<=4000 ||
        throw(ArgumentError("Reference cells must be integer 4..4000."))
    isfinite(t_end_s) && t_end_s>0 && isfinite(dt_s) && dt_s>0 ||
        throw(ArgumentError("Require positive finite reference duration/time step."))
    h=L/cells
    v*h/D<=1 || throw(ArgumentError("Reference cell Peclet must be <=1; refine cells."))
    steps=ceil(Int,t_end_s/dt_s)
    steps<=200000 || throw(ArgumentError("At most 200000 reference steps."))
    dt=t_end_s/steps
    lower=fill(-dt*(D/h^2+v/(2h)),cells-1)
    upper=fill(dt*(v/(2h)-D/h^2),cells-1)
    diagonal=fill(1+2dt*D/h^2,cells)
    diagonal[1]=1+dt*(3D/h^2+v/(2h))
    diagonal[end]=1+dt*(D/h^2+v/(2h))
    F=lu(Tridiagonal(lower,diagonal,upper))
    c=zeros(cells);time=[0.];exit=[0.];input=[0.];output=[0.];inventory=[0.]
    incoming=0.;outgoing=0.;maxbalance=0.
    for i in 1:steps
        c[1]+=dt*(v/h+2D/h^2);c=F\c
        ji=v+2D/h*(1-c[1]);jo=v*c[end]
        incoming+=dt*ji;outgoing+=dt*jo;mass=h*sum(c)
        maxbalance=max(maxbalance,abs(mass-incoming+outgoing))
        push!(time,i*dt);push!(exit,jo);push!(input,incoming)
        push!(output,outgoing);push!(inventory,mass)
    end
    (;time_s=time,outlet_flux_per_area_m_s=exit,inventory_per_area_m=inventory,
      cumulative_inlet_per_area_m=input,cumulative_outlet_per_area_m=output,
      final_concentration=c,normalized_mass_residual=maxbalance/max(incoming,1e-30))
end
end
