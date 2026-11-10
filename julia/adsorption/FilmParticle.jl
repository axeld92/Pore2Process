module FilmParticle
using LinearAlgebra
export equilibrium, bath_equilibrium, film_coefficient, initial_film_rate, solve_bath, solve_column, solve_fixed_bulk

function positive(name, value; zero=false)
    value isa Real && !(value isa Bool) && isfinite(value) && (zero ? value>=0 : value>0) ||
        throw(ArgumentError("$name must be finite and $(zero ? "nonnegative" : "positive") in SI units."))
    v=Float64(value)
    isfinite(v) || throw(ArgumentError("$name exceeds Float64 range."))
    v
end

function count(name,n,lo,hi)
    n isa Integer && !(n isa Bool) && lo<=n<=hi || throw(ArgumentError("$name must be an integer in $lo..$hi."))
    Int(n)
end

function isotherm(qmax,b,linear)
    qm=positive("capacity_kg_kg",qmax;zero=true)
    k=positive("affinity_m3_kg",b;zero=true)
    K=linear===nothing ? nothing : positive("linear_partition_m3_kg",linear;zero=true)
    qe(c)=K===nothing ? qm*(k*c)/(1+k*c) : K*c
    slope(c)=K===nothing ? qm*k/(1+k*c)^2 : K
    qe,slope
end

function equilibrium(C;capacity_kg_kg=.1115,affinity_m3_kg=139.,linear_partition_m3_kg=nothing)
    c=positive("concentration_kg_m3",C;zero=true)
    qe,_=isotherm(capacity_kg_kg,affinity_m3_kg,linear_partition_m3_kg)
    qe(c)
end

"""Finite-bath endpoint from V*C+m*q(C)=V*C0. No kinetic fit is involved."""
function bath_equilibrium(;initial_concentration_kg_m3=.2,fluid_volume_m3=.00025,
        adsorbent_mass_kg=.001,capacity_kg_kg=.1115,affinity_m3_kg=139.,linear_partition_m3_kg=nothing)
    c0=positive("initial_concentration_kg_m3",initial_concentration_kg_m3;zero=true)
    V=positive("fluid_volume_m3",fluid_volume_m3);m=positive("adsorbent_mass_kg",adsorbent_mass_kg)
    qe,_=isotherm(capacity_kg_kg,affinity_m3_kg,linear_partition_m3_kg)
    lo=0.;hi=c0
    for _ in 1:90
        c=(lo+hi)/2
        if V*c+m*qe(c)>V*c0;hi=c;else;lo=c;end
    end
    c=(lo+hi)/2
    (;concentration_kg_m3=c,loading_kg_kg=qe(c),mass_residual_kg=V*c+m*qe(c)-V*c0)
end

"""Declared Wakao/Funazkri-form scenario; velocity convention is explicit."""
function film_coefficient(;diameter_m=.0008,superficial_velocity_m_s=.00768,
        viscosity_Pa_s=1/1200,density_kg_m3=1000.,molecular_diffusivity_m2_s=1e-9,
        porosity=.427,interstitial=false,multiplier=1.)
    d=positive("diameter_m",diameter_m);u=positive("superficial_velocity_m_s",superficial_velocity_m_s;zero=true)
    mu=positive("viscosity_Pa_s",viscosity_Pa_s);rho=positive("density_kg_m3",density_kg_m3)
    D=positive("molecular_diffusivity_m2_s",molecular_diffusivity_m2_s)
    eps=positive("porosity",porosity);eps<1 || throw(ArgumentError("porosity must be < 1."))
    interstitial isa Bool || throw(ArgumentError("interstitial must be boolean."))
    f=positive("multiplier",multiplier)
    v=interstitial ? u/eps : u
    Re=rho*v*d/mu;Sc=mu/(rho*D);Sh=2+1.1*Re^.6*cbrt(Sc)
    (;film_m_s=f*Sh*D/d,Reynolds=Re,Schmidt=Sc,Sherwood=Sh,velocity_m_s=v)
end

"""Clean-particle t→0 film-supply limit; not an observed initial uptake rate."""
function initial_film_rate(;particle_diameter_m=.0008,particle_density_kg_m3=740/(1-.427),
        film_m_s=1e-5,concentration_kg_m3=.2)
    R=positive("particle_diameter_m",particle_diameter_m)/2
    rho=positive("particle_density_kg_m3",particle_density_kg_m3)
    k=positive("film_m_s",film_m_s;zero=true);c=positive("concentration_kg_m3",concentration_kg_m3;zero=true)
    (;initial_rate_kg_kg_s=3k*c/(rho*R))
end

function operator(n,R,Ds,kf,dt)
    edges=collect(range(0.,R;length=n+1));dr=R/n
    w=diff(edges.^3)/R^3
    g=3Ds.*edges[2:end-1].^2/(R^3*dr)
    outer=Ds>0 && kf>0 ? 3Ds/(R*(dr/2)) : 0.
    diag=ones(n);lower=zeros(n-1);upper=zeros(n-1)
    for i in 1:n-1
        diag[i]+=dt*g[i]/w[i];diag[i+1]+=dt*g[i]/w[i+1]
        upper[i]=-dt*g[i]/w[i];lower[i]=-dt*g[i]/w[i+1]
    end
    diag[end]+=dt*outer/w[end]
    F=lu(Tridiagonal(lower,diag,upper))
    rhs=zeros(n);rhs[end]=dt*outer/w[end]
    b=F\rhs
    (;F,b,w,mean_b=dot(w,b),last_factor=1-b[end],half_width=dr/2,
      radius_m=(edges[1:end-1]+edges[2:end])/2)
end

"""Monotone implicit surface root, in liquid concentration units."""
function surface(c,a_last,op,kf,rho,Ds,qe,slope)
    if kf==0 || Ds==0
        return c,0.
    end
    theta=rho*Ds/(op.half_width*kf)
    lo=0.;hi=c+theta*max(a_last,0.)
    x=min(c,hi)
    tol=2e-13*max(c,1e-6)
    for _ in 1:80
        f=x-c+theta*(op.last_factor*qe(x)-a_last)
        derivative=1+theta*op.last_factor*slope(x)
        if abs(f)<=tol || hi-lo<=tol
            return x,1/derivative
        end
        if f>0;hi=x;else;lo=x;end
        candidate=x-f/derivative
        x=lo<candidate<hi ? candidate : (lo+hi)/2
    end
    throw(ErrorException("Surface concentration root failed to converge."))
end

function local_response!(meanq,derivative,cs,c,a,mean_a,op,kf,rho,Ds,qe,slope)
    for j in eachindex(c)
        cs[j],dcs=surface(c[j],a[end,j],op,kf,rho,Ds,qe,slope)
        meanq[j]=mean_a[j]+op.mean_b*qe(cs[j])
        derivative[j]=op.mean_b*slope(cs[j])*dcs
    end
end

"""
Conservative spherical loading diffusion, coupled to a bath or an axial bed.

q is kg solute/kg dry active solid; intraparticle liquid inventory is neglected.
The surface q_s=qeq(Cs) is connected to the external liquid by the SAME flux
J=kf(C-Cs)=rho_particle*Ds*dq/dr. Shell volume weights sum to one. Particle
states are eliminated from the implicit axial Newton system, not fitted or
advanced with a separate nonconservative source step.
"""
function integrate_impl(cin_flow;nz,nr,L,A,storage,bulk_density,physical_mass,active_mass,
        R,rho,Ds,kf,cin,cinit,qinit,dispersion,tend,dt,qe,slope,
        bath=false,fluid_volume=0.,sample_times=Float64[],sample_volume=0.,fixed_bulk=false)
    dz=L/nz;volume=A*L
    c=fill(cinit,nz);q=fill(qinit,nr,nz)
    first_dt=min(dt,tend)
    ops=Dict(first_dt=>operator(nr,R,Ds,kf,first_dt))
    getop(h)=get!(()->operator(nr,R,Ds,kf,h),ops,h)
    op=getop(min(dt,tend));w=op.w
    meanq=vec(transpose(w)*q)
    initial_mass=bath ? fluid_volume*cinit+active_mass*qinit : volume*(storage*cinit+bulk_density*qinit)
    times=[0.];outlet=[c[end]];loadings=[qinit];inventories=[initial_mass]
    solid_history=[active_mass*qinit];fluid_history=[initial_mass-active_mass*qinit]
    volumes=[bath ? fluid_volume : storage*volume]
    input_history=[0.];output_history=[0.];sample_history=[0.];residuals=[0.]
    input=0.;output=0.;sample_output=0.;maxbalance=0.;maxfluxerror=0.;maxparticleerror=0.
    maxiterations=0;step=0;t=0.;event=1;started=time_ns();cs=similar(c)
    while t<tend-1e-10*max(1.,tend)
        h=min(dt,tend-t)
        if event<=length(sample_times);h=min(h,sample_times[event]-t);end
        h>0 || throw(ErrorException("Sampling events must be strictly increasing and > 0."))
        step+=1;step<=200000 || throw(ArgumentError("More than 200000 steps; increase dt or shorten the window."))
        op=getop(h);a=op.F\q;mean_a=vec(transpose(w)*a)
        old_mean=copy(meanq);old_c=copy(c)
        beta=bath ? active_mass/fluid_volume : bulk_density
        trial=copy(c);qbar=similar(c);derivative=similar(c)
        candidate=similar(c);candidate_cs=similar(c);candidate_qbar=similar(c);candidate_derivative=similar(c)
        res=similar(c);candidate_res=similar(c)
        function residual!(guess,result,loading,deriv,surface_c)
            local_response!(loading,deriv,surface_c,guess,a,mean_a,op,kf,rho,Ds,qe,slope)
            result.=storage*(guess-old_c)+beta*(loading-old_mean)
            if !bath && !fixed_bulk
                for j in 1:nz
                    # Total inlet flux is u*cin; outlet dispersive flux is zero.
                    left=j==1 ? cin_flow*cin : cin_flow*guess[j-1]-storage*dispersion*(guess[j]-guess[j-1])/dz
                    right=j==nz ? cin_flow*guess[j] : cin_flow*guess[j]-storage*dispersion*(guess[j+1]-guess[j])/dz
                    result[j]+=h/dz*(right-left)
                end
            end
        end
        # Fixed-bulk particle oracle deliberately has an unlimited liquid reservoir.
        if fixed_bulk
            local_response!(qbar,derivative,cs,trial,a,mean_a,op,kf,rho,Ds,qe,slope)
        else
            converged=false
            for iteration in 1:50
                residual!(trial,res,qbar,derivative,cs)
                maxiterations=max(maxiterations,iteration)
                norm=maximum(abs,res)
                # Resolve late, small bath transfers as well as early uptake:
                # a looser per-step stop accumulated mass error at smaller dt.
                if norm<=1e-13*max(cin,cinit,1e-6);converged=true;break;end
                adv=bath ? 0. : h*cin_flow/dz
                diff=bath ? 0. : h*storage*dispersion/dz^2
                diag=storage.+beta*derivative.+adv.+2diff
                diag[1]-=diff;diag[end]-=diff
                delta=Tridiagonal(fill(-adv-diff,nz-1),diag,fill(-diff,nz-1))\(-res)
                accepted=false;alpha=1.
                for _ in 1:40
                    candidate.=max.(0.,trial+alpha*delta)
                    residual!(candidate,candidate_res,candidate_qbar,candidate_derivative,candidate_cs)
                    if maximum(abs,candidate_res)<=(1-1e-4*alpha)*norm
                        trial.=candidate;accepted=true;break
                    end
                    alpha/=2
                end
                accepted || throw(ErrorException("Coupled film/particle Newton line search stalled; reduce dt."))
            end
            converged || throw(ErrorException("Coupled film/particle step did not converge; reduce dt."))
        end
        c.=trial
        q.=a+op.b*transpose(qe.(cs))
        meanq.=vec(transpose(w)*q)
        minimum(c)>=0 && minimum(q)>=0 || throw(ErrorException("Negative completed concentration/loading."))
        for j in 1:nz
            film=kf*(c[j]-cs[j])
            radial=Ds==0 || kf==0 ? 0. : rho*Ds/op.half_width*(qe(cs[j])-q[end,j])
            maxfluxerror=max(maxfluxerror,abs(film-radial))
            maxparticleerror=max(maxparticleerror,abs((meanq[j]-old_mean[j])/h-3film/(rho*R)))
        end
        t=min(t+h,tend)
        if !bath && !fixed_bulk
            input+=cin_flow*A*cin*h;output+=cin_flow*A*c[end]*h
        end
        if bath && event<=length(sample_times) && abs(t-sample_times[event])<=1e-9*max(1.,t)
            sample_output+=sample_volume*c[1];fluid_volume-=sample_volume;event+=1
        end
        if fixed_bulk
            input=active_mass*meanq[1]-active_mass*qinit
        end
        solid=active_mass*sum(meanq)/nz
        liquid=bath ? fluid_volume*c[1] : storage*volume*sum(c)/nz
        inventory=solid+liquid
        balance=inventory-initial_mass-input+output+sample_output
        maxbalance=max(maxbalance,abs(balance))
        push!(times,t);push!(outlet,c[end]);push!(loadings,sum(meanq)/nz)
        push!(inventories,inventory);push!(solid_history,solid);push!(fluid_history,liquid)
        push!(input_history,input);push!(output_history,output);push!(sample_history,sample_output)
        push!(residuals,balance);push!(volumes,bath ? fluid_volume : storage*volume)
    end
    (;time_s=times,outlet_concentration_kg_m3=outlet,mean_loading_kg_kg=loadings,
      final_concentration_kg_m3=copy(c),final_shell_loading_kg_kg=permutedims(q),
      final_surface_concentration_kg_m3=copy(cs),radial_centres_m=op.radius_m,
      radial_volume_weights=copy(w),axial_centres_m=collect(((1:nz).-0.5)*dz),
      inventory_kg=inventories,adsorbed_inventory_kg=solid_history,liquid_inventory_kg=fluid_history,
      cumulative_inlet_kg=input_history,cumulative_outlet_kg=output_history,
      cumulative_sampled_solute_kg=sample_history,fluid_volume_m3=volumes,mass_residual_kg=residuals,
      physical_adsorbent_mass_kg=physical_mass,active_adsorbent_mass_kg=active_mass,
      normalized_mass_residual=maxbalance/max(initial_mass,input,1e-30),max_mass_residual_kg=maxbalance,
      max_interface_flux_mismatch_kg_m2_s=maxfluxerror,
      normalized_interface_flux_mismatch=maxfluxerror/max(kf*max(cin,cinit),1e-30),
      max_particle_balance_error_kg_kg_s=maxparticleerror,
      axial_cells=nz,radial_cells=nr,steps=step,max_Newton_iterations=maxiterations,total_s=(time_ns()-started)/1e9)
end

function parameters(;particle_diameter_m=.0008,particle_density_kg_m3=740/(1-.427),
        diffusivity_m2_s=1e-12,film_m_s=1e-5,capacity_kg_kg=.1115,affinity_m3_kg=139.,
        linear_partition_m3_kg=nothing,radial_cells=12,t_end_s=180000.,dt_s=60.,initial_loading_kg_kg=0.)
    R=positive("particle_diameter_m",particle_diameter_m)/2
    rho=positive("particle_density_kg_m3",particle_density_kg_m3)
    Ds=positive("diffusivity_m2_s",diffusivity_m2_s;zero=true)
    kf=positive("film_m_s",film_m_s;zero=true)
    nr=count("radial_cells",radial_cells,1,128)
    tend=positive("t_end_s",t_end_s);dt=positive("dt_s",dt_s)
    tend/dt<=200000 || throw(ArgumentError("More than 200000 time steps; increase dt or shorten duration."))
    qinit=positive("initial_loading_kg_kg",initial_loading_kg_kg;zero=true)
    qe,slope=isotherm(capacity_kg_kg,affinity_m3_kg,linear_partition_m3_kg)
    (;R,rho,Ds,kf,nr,tend,dt,qinit,qe,slope)
end

function solve_bath(;initial_concentration_kg_m3=.2,fluid_volume_m3=.00025,
        adsorbent_mass_kg=.001,sampling_times_s=Float64[],sample_volume_m3=0.,kwargs...)
    p=parameters(;kwargs...)
    c0=positive("initial_concentration_kg_m3",initial_concentration_kg_m3;zero=true)
    V=positive("fluid_volume_m3",fluid_volume_m3);m=positive("adsorbent_mass_kg",adsorbent_mass_kg)
    sample_volume=positive("sample_volume_m3",sample_volume_m3;zero=true)
    events=Float64.(collect(sampling_times_s))
    all(isfinite,events) && all(>(0),events) && all(>(0),diff(events)) && all(<=(p.tend),events) ||
        throw(ArgumentError("Sampling times must be finite, strictly increasing and in (0,t_end_s]."))
    length(events)*sample_volume<V || throw(ArgumentError("Sampling would exhaust the bath."))
    integrate_impl(0.;nz=1,L=1.,A=V,storage=1.,bulk_density=m/V,physical_mass=m,active_mass=m,
        cin=c0,cinit=c0,dispersion=0.,bath=true,fluid_volume=V,sample_times=events,
        sample_volume=sample_volume,p...)
end

function solve_column(;length_m=.1,area_m2=pi*.011^2/4,porosity=.427,
        superficial_velocity_m_s=43.8e-6/60/(pi*.011^2/4),inlet_concentration_kg_m3=.15,
        axial_dispersion_m2_s=0.,axial_cells=80,active_inventory_fraction=1.,
        initial_concentration_kg_m3=0.,kwargs...)
    p=parameters(;kwargs...)
    L=positive("length_m",length_m);A=positive("area_m2",area_m2)
    eps=positive("porosity",porosity);eps<1 || throw(ArgumentError("porosity must be < 1."))
    f=positive("active_inventory_fraction",active_inventory_fraction;zero=true)
    f<=1 || throw(ArgumentError("active_inventory_fraction must be <= 1."))
    u=positive("superficial_velocity_m_s",superficial_velocity_m_s;zero=true)
    cin=positive("inlet_concentration_kg_m3",inlet_concentration_kg_m3;zero=true)
    c0=positive("initial_concentration_kg_m3",initial_concentration_kg_m3;zero=true)
    D=positive("axial_dispersion_m2_s",axial_dispersion_m2_s;zero=true)
    nz=count("axial_cells",axial_cells,2,2000)
    physical=(1-eps)*p.rho*A*L;active=f*physical
    integrate_impl(u;nz,L,A,storage=eps,bulk_density=f*(1-eps)*p.rho,
        physical_mass=physical,active_mass=active,cin,cinit=c0,dispersion=D,p...)
end

"""Particle oracle with a constant external concentration and explicit reservoir uptake."""
function solve_fixed_bulk(;concentration_kg_m3=.2,kwargs...)
    p=parameters(;kwargs...)
    c=positive("concentration_kg_m3",concentration_kg_m3;zero=true)
    integrate_impl(0.;nz=1,L=1.,A=1.,storage=1.,bulk_density=1.,physical_mass=1.,active_mass=1.,
        cin=c,cinit=c,dispersion=0.,fixed_bulk=true,p...)
end

end
