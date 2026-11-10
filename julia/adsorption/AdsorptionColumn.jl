module AdsorptionColumn
using LinearAlgebra
export solve_column, equilibrium_loading, ldf_rate, pressure_load

function number(name,x; zero=false)
    x isa Real && !(x isa Bool) && isfinite(x) && (zero ? x>=0 : x>0) ||
        throw(ArgumentError("$name must be finite and $(zero ? "nonnegative" : "positive") in SI units."))
    y=Float64(x);isfinite(y) || throw(ArgumentError("$name exceeds Float64 range."));y
end

"""kg solute/kg dry solid. b is m³/kg solute; C is kg solute/m³ liquid."""
function equilibrium_loading(C,qmax,b; linear_partition=nothing)
    c=number("concentration",C;zero=true)
    if linear_partition!==nothing
        return number("linear partition_m3_kg",linear_partition;zero=true)*c
    end
    q=number("capacity_kg_kg",qmax;zero=true);k=number("affinity_m3_kg",b;zero=true)
    z=k*c
    isfinite(z) || throw(ArgumentError("Affinity*concentration overflows."))
    q*z/(1+z)
end

"""Spherical intraparticle LDF approximation, 15 Deff/R²; not a fitted film law."""
ldf_rate(D,diameter)=number("LDF rate",60number("effective diffusivity_m2_s",D)/number("particle diameter_m",diameter)^2)

"""Darcy load from a separately measured native permeability; no Ergun substitution."""
function pressure_load(K,mu,L,u,A,efficiency)
    k=number("native permeability_m2",K);m=number("viscosity_Pa_s",mu);l=number("length_m",L)
    v=number("superficial velocity_m_s",u);a=number("area_m2",A);e=number("efficiency",efficiency)
    e<=1 || throw(ArgumentError("Pump efficiency must be <= 1."))
    dp=number("pressure drop_Pa",m*l*v/k);Q=number("flow_m3_s",v*a)
    (;pressure_drop_Pa=dp,flow_m3_s=Q,reference_shaft_power_W=number("power_W",Q*dp/e))
end

"""
Conservative cell-centred backward-Euler finite-volume adsorption column.

eps*dC/dt + d(u*C-eps*Dax*dC/dz)/dz = -rho_bulk*dq/dt
dq/dt = k_LDF*(qeq(C)-q).

Danckwerts inlet total flux=u*Cin; outlet zero dispersive flux. q is kg/kg
dry solid; rho_particle is dry mass per particle-envelope volume, distinct
from bed bulk density (1-eps)*rho_particle. Carrier flow is incompressible and
solute is dilute. No thermal, reaction, regeneration or film closure is implied.
Optional linear equilibrium and zero-flow batch cases are numerical oracles.
"""
function solve_column(;length_m=1.,area_m2=.01,porosity=.6,particle_density_kg_m3=1000.,
        superficial_velocity_m_s=1e-4,inlet_concentration_kg_m3=.05,
        capacity_kg_kg=.04169,affinity_m3_kg=94.,ldf_s_inv=1e-4,
        axial_dispersion_m2_s=0.,t_end_s=3e6,cells=160,dt_s=500.,
        initial_concentration_kg_m3=0.,initial_loading_kg_kg=0.,
        linear_partition_m3_kg=nothing,equilibrium=false,breakthrough_fraction=.05)
    L=number("length_m",length_m);A=number("area_m2",area_m2)
    eps=number("porosity",porosity);eps<1 || throw(ArgumentError("Porosity must be in (0,1)."))
    rho=number("particle_density_kg_m3",particle_density_kg_m3)
    u=number("superficial_velocity_m_s",superficial_velocity_m_s;zero=true)
    cin=number("inlet_concentration_kg_m3",inlet_concentration_kg_m3)
    qm=number("capacity_kg_kg",capacity_kg_kg;zero=true);b=number("affinity_m3_kg",affinity_m3_kg;zero=true)
    rate=number("ldf_s_inv",ldf_s_inv;zero=true);D=number("axial_dispersion_m2_s",axial_dispersion_m2_s;zero=true)
    tend=number("t_end_s",t_end_s);dt=number("dt_s",dt_s)
    cinit=number("initial_concentration_kg_m3",initial_concentration_kg_m3;zero=true)
    qinit=number("initial_loading_kg_kg",initial_loading_kg_kg;zero=true)
    threshold=number("breakthrough_fraction",breakthrough_fraction);threshold<1 || throw(ArgumentError("Breakthrough fraction must be < 1."))
    cells isa Integer && !(cells isa Bool) && 2<=cells<=4000 || throw(ArgumentError("Use 2..4000 axial cells."))
    equilibrium isa Bool || throw(ArgumentError("equilibrium must be boolean."))
    K=linear_partition_m3_kg===nothing ? nothing : number("linear_partition_m3_kg",linear_partition_m3_kg;zero=true)
    steps=ceil(Int,tend/dt);steps<=200000 || throw(ArgumentError("Too many time steps; increase dt or reduce duration."))
    dt=tend/steps;n=Int(cells);dz=L/n;rb=(1-eps)*rho
    qe(c)=K===nothing ? qm*(b*c)/(1+b*c) : K*c
    dqe(c)=K===nothing ? qm*b/(1+b*c)^2 : K
    c=fill(cinit,n);q=fill(qinit,n);times=[0.];outlet=[c[end]]
    stored0=A*dz*sum(eps*c+rb*q);input=0.;output=0.;maximum_residual=0.;newton_max=0
    inventories=[stored0];input_history=[0.];output_history=[0.];residual_history=[0.]
    solid_history=[A*dz*rb*sum(q)];liquid_history=[A*dz*eps*sum(c)]
    f=zeros(n+1);res=zeros(n);jacdiag=zeros(n);lower=zeros(n-1);upper=zeros(n-1)
    trial=copy(c);qnew=copy(q)
    candidate=copy(c);candidate_q=copy(q);candidate_res=copy(c)
    function residual!(guess,loading,result,factor)
        loading.=equilibrium ? qe.(guess) : (1-factor)*q+factor*qe.(guess)
        f[1]=u*cin;f[end]=u*guess[end]
        for i in 2:n;f[i]=u*guess[i-1]-eps*D*(guess[i]-guess[i-1])/dz;end
        result.=eps*(guess-c)+dt/dz*(f[2:end]-f[1:end-1])+rb*(loading-q)
    end
    start=time_ns()
    for step in 1:steps
        trial.=c
        factor=equilibrium ? 1. : dt*rate/(1+dt*rate)
        adv=dt*u/dz;diff=dt*eps*D/dz^2
        converged=false
        for iteration in 1:40
            residual!(trial,qnew,res,factor)
            newton_max=max(newton_max,iteration)
            if maximum(abs,res)<=1e-12*max(cin,cinit,1e-8)
                converged=true;break
            end
            jacdiag.=(eps+adv+2diff).+rb*factor*dqe.(trial)
            jacdiag[1]-=diff;jacdiag[end]-=diff
            lower.=-adv-diff;upper.=-diff
            delta=Tridiagonal(lower,jacdiag,upper)\(-res)
            # Bound-constrained Newton search. Project trial ITERATES, then
            # re-evaluate the original conservative residual before acceptance.
            # This handles underflow-size negative increments in a clean-bed
            # tail without clipping a completed solution or correcting its mass.
            alpha=1.;accepted=false;oldnorm=maximum(abs,res)
            for backtrack in 1:40
                candidate.=max.(0.,trial+alpha*delta)
                residual!(candidate,candidate_q,candidate_res,factor)
                if maximum(abs,candidate_res)<=(1-1e-4*alpha)*oldnorm
                    trial.=candidate;accepted=true;break
                end
                alpha/=2
            end
            accepted || throw(ErrorException("Implicit residual line search stalled; reduce dt."))
        end
        converged || throw(ErrorException("Implicit adsorption step did not converge; reduce dt."))
        c.=trial;q.=qnew
        minimum(c)>=0 && minimum(q)>=0 || throw(ErrorException("Negative state."))
        input+=u*A*cin*dt;output+=u*A*c[end]*dt
        stored=A*dz*sum(eps*c+rb*q);balance=stored-stored0-input+output
        maximum_residual=max(maximum_residual,abs(balance))
        push!(times,step*dt);push!(outlet,c[end]);push!(inventories,stored)
        push!(input_history,input);push!(output_history,output);push!(residual_history,balance)
        push!(solid_history,A*dz*rb*sum(q));push!(liquid_history,A*dz*eps*sum(c))
    end
    breakthrough=nothing
    captured=nothing;utilization=nothing
    if u>0
        for i in eachindex(outlet)
            if outlet[i]>=threshold*cin
                weight=i==1 ? 0. : (threshold*cin-outlet[i-1])/(outlet[i]-outlet[i-1])
                breakthrough=i==1 ? 0. : times[i-1]+weight*(times[i]-times[i-1])
                captured=i==1 ? solid_history[1] : solid_history[i-1]+weight*(solid_history[i]-solid_history[i-1])
                utilization=qe(cin)==0 ? nothing : captured/(rb*A*L*qe(cin))
                break
            end
        end
    end
    (;x_m=collect(((1:n).-0.5)*dz),time_s=times,outlet_concentration_kg_m3=outlet,
      final_concentration_kg_m3=copy(c),final_loading_kg_kg=copy(q),inventory_kg=inventories,
      cumulative_inlet_kg=input_history,cumulative_outlet_kg=output_history,mass_residual_kg=residual_history,
      adsorbed_inventory_kg=solid_history,liquid_inventory_kg=liquid_history,
      max_mass_residual_kg=maximum_residual,normalized_mass_residual=maximum_residual/max(stored0,input,1e-30),
      breakthrough_time_s=breakthrough,treated_volume_m3=breakthrough===nothing ? nothing : u*A*breakthrough,
      adsorbed_solute_at_breakthrough_kg=captured,solid_utilization_at_breakthrough=utilization,
      adsorbent_mass_kg=rb*A*L,fluid_holdup_m3=eps*A*L,bulk_density_kg_m3=rb,
      equilibrium_feed_loading_kg_kg=qe(cin),sharp_front_capacity_time_s=u==0 ? nothing : L/u*(eps+rb*qe(cin)/cin),
      dt_s=dt,cells=n,steps=steps,max_Newton_iterations=newton_max,total_s=(time_ns()-start)/1e9)
end
end
