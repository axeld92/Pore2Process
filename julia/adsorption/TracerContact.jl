module TracerContact
using LinearAlgebra
export absorbing_moments,step_moments,derive_contact,solve_ade

positive(x)=x isa Real && !(x isa Bool) && isfinite(x) && x>0 ? Float64(x) : throw(ArgumentError("Require positive finite SI inputs."))
function absorbing_moments(;length_m,velocity_m_s,dispersion_m2_s)
    L=positive(length_m);v=positive(velocity_m_s);D=positive(dispersion_m2_s);Pe=v*L/D
    Pe>=.05 || throw(ArgumentError("Peclet >=0.05 required by this moment evaluator."))
    coth=1+2/expm1(Pe);csch2=Pe>700 ? 0. : 4exp(-Pe)/(1-exp(-Pe))^2
    m=L/v*(coth-2/Pe)
    variance=(L/v)^2*(2/Pe*coth+csch2-8/Pe^2)
    (;mean_s=m,variance_s2=variance,Peclet=Pe,steady_flux_per_area_m_s=v/(-expm1(-Pe)))
end

trap(t,y)=sum(diff(t).*(y[1:end-1]+y[2:end])/2)
"""Moments of the native normalized EXIT-FLUX step, not adsorption data.
No clipping, deconvolution, residual correction or extrapolated tail.
"""
function step_moments(time_s,flux_kg_s)
    t=Float64.(collect(time_s));f=Float64.(collect(flux_kg_s))
    length(t)==length(f)>=20 && all(isfinite,t) && all(isfinite,f) && all(>(0),diff(t)) && t[1]>=0 || throw(ArgumentError("Need >=20 finite ordered times and equal-length exit flux."))
    plateau=sum(f[round(Int,.9length(f)):end])/length(f[round(Int,.9length(f)):end])
    positive(plateau)
    # First native flux is at dt/2; clean origin is declared, not hidden delay.
    tt=vcat(0.,t);F=vcat(0.,f/plateau);S=1 .-F
    m=trap(tt,S);second=2trap(tt,tt.*S);variance=second-m*m
    m>0 && variance>0 || throw(ArgumentError("Nonpositive finite-window RTD moments; extend or inspect the native run."))
    tail=findall(>=(tt[end]/2),tt)
    tail_fraction=abs(trap(tt[tail],S[tail]))/m
    (;mean_s=m,variance_s2=variance,steady_exit_flux_kg_s=plateau,
      final_fraction=F[end],tail_area_fraction=tail_fraction,
      tail_second_moment_fraction=abs(2trap(tt[tail],tt[tail].*S[tail]))/second,
      equivalent_sections=m*m/variance,minimum_fraction=minimum(F),maximum_fraction=maximum(F))
end

"""Two moment-derived, boundary-specific contact coefficients from inert tracer.
Effective velocity is NOT silently substituted for Q/physical pore area.
CV^2>=0.4 cannot be represented by this absorbing-end homogeneous ADE family.
"""
function derive_contact(;length_m,mean_s,variance_s2)
    L=positive(length_m);m=positive(mean_s);var=positive(variance_s2);cv2=var/m^2
    lo=.05;hi=1e6
    cv(p)=begin
        a=absorbing_moments(length_m=1.,velocity_m_s=1.,dispersion_m2_s=1/p)
        a.variance_s2/a.mean_s^2
    end
    cv(hi)<cv2<cv(lo) || throw(ArgumentError("Native CV^2 outside representable absorbing-ADE range; no contact closure accepted."))
    for _ in 1:90
        p=sqrt(lo*hi)
        if cv(p)>cv2;lo=p;else;hi=p;end
    end
    Pe=sqrt(lo*hi);factor=1+2/expm1(Pe)-2/Pe
    v=L/m*factor;D=v*L/Pe
    (;effective_velocity_m_s=v,effective_dispersion_m2_s=D,Peclet=Pe,equivalent_sections=m*m/var)
end

"""Independent central finite-volume ADE, Dirichlet step inlet / absorbing exit.
This Julia reference is not a second pore-scale solver. SI, clean initial state.
"""
function solve_ade(;length_m,velocity_m_s,dispersion_m2_s,t_end_s,dt_s=.01,cells=400)
    L=positive(length_m);v=positive(velocity_m_s);D=positive(dispersion_m2_s)
    tend=positive(t_end_s);dt=positive(dt_s)
    cells isa Integer && !(cells isa Bool) && 4<=cells<=4000 || throw(ArgumentError("cells must be integer 4..4000"))
    dz=L/cells;v*dz/D<=1 || throw(ArgumentError("Refine cells: cell Peclet must be <=1."))
    steps=ceil(Int,tend/dt);steps<=200000 || throw(ArgumentError("At most 200000 reference steps."))
    dt=tend/steps
    lower=fill(-dt*(D/dz^2+v/(2dz)),cells-1)
    upper=fill(dt*(v/(2dz)-D/dz^2),cells-1)
    diag=fill(1+2dt*D/dz^2,cells)
    diag[1]=1+dt*(3D/dz^2+v/(2dz));diag[end]=1+dt*(3D/dz^2-v/(2dz))
    factor=lu(Tridiagonal(lower,diag,upper));rhs=dt*(2D/dz^2+v/dz)
    c=zeros(cells);times=[0.];exit=[0.];inventories=[0.];input=[0.];output=[0.]
    cumulative_in=0.;cumulative_out=0.;maxbalance=0.
    for i in 1:steps
        c[1]+=rhs;c=factor\c
        ji=v+2D/dz*(1-c[1]);jo=2D/dz*c[end]
        cumulative_in+=dt*ji;cumulative_out+=dt*jo;mass=dz*sum(c)
        maxbalance=max(maxbalance,abs(mass-cumulative_in+cumulative_out))
        push!(times,i*dt);push!(exit,jo);push!(inventories,mass);push!(input,cumulative_in);push!(output,cumulative_out)
    end
    (;time_s=times,outlet_flux_per_area_m_s=exit,inventory_per_area_m=inventories,
      cumulative_inlet_per_area_m=input,cumulative_outlet_per_area_m=output,
      normalized_mass_residual=maxbalance/cumulative_in,final_concentration=copy(c),cells,dt_s=dt)
end
end
