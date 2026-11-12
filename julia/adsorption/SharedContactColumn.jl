module SharedContactColumn

using LinearAlgebra
using ..FilmParticle
using ..ExchangeColumn
export solve_column

"""Invert s + beta*qeq(s) = rhs using the specified batch isotherm."""
@inline function surface_value(rhs,beta,qmax,b,linear)
    if linear!==nothing
        return rhs/(1+beta*linear)
    end
    b==0 && return rhs
    x = 1+b*(beta*qmax-rhs)
    d = sqrt(muladd(x,x,4b*rhs))
    x>=0 ? 2rhs/(x+d) : (d-x)/(2b)
end

function particle_surface!(cs,qs,free,beta,weights,qprev,a_last,outer,cap,qmax,b,linear,qe,slope)
    lo,hi = 0.,cap
    q = clamp(qprev,lo,hi)
    for iteration in 1:40
        mean,derivative = 0.,0.
        @inbounds for i in eachindex(weights)
            weights[i]==0 && continue
            cs[i] = surface_value(free[i]+beta[i]*q,beta[i],qmax,b,linear)
            qs[i] = qe(cs[i])
            mean += weights[i]*qs[i]
            g = beta[i]*slope(cs[i])
            derivative += weights[i]*g/(1+g)
        end
        residual = q-a_last-outer*mean
        if abs(residual)<=2e-13*max(cap,1e-6)
            return mean,iteration
        end
        if residual>0
            hi = q
        else
            lo = q
        end
        candidate = q-residual/(1-outer*derivative)
        q = lo<candidate<hi ? candidate : (lo+hi)/2
    end
    error("Shared particle surface root failed to converge")
end

"""
Each axial slice has one shared radial particle profile and separate path
film concentrations. Wall exposure assigns interfacial area. Liquid paths
mix only through the supplied exchange matrix at each period boundary.
"""
function solve_column(flow_fractions,relative_delays,receiving_mix,contact_fractions;
        periods,subcells=1,length_m,area_m2,porosity,superficial_velocity_m_s,
        inlet_concentration_kg_m3,startup_s=0.,startup_dt_s=.05,
        capacity_kg_kg=.1115,affinity_m3_kg=139.,linear_partition_m3_kg=nothing,kwargs...)
    p = FilmParticle.parameters(;capacity_kg_kg,affinity_m3_kg,linear_partition_m3_kg,kwargs...)
    p.qinit==0 || throw(ArgumentError("Shared-contact column requires initially clean carbon"))
    count = FilmParticle.count("periods",periods,1,2000)
    slices = FilmParticle.count("subcells",subcells,1,64)
    ns = count*slices
    L = FilmParticle.positive("length_m",length_m)
    A = FilmParticle.positive("area_m2",area_m2)
    eps = FilmParticle.positive("porosity",porosity)
    eps<1 || throw(ArgumentError("porosity must be < 1"))
    u = FilmParticle.positive("superficial_velocity_m_s",superficial_velocity_m_s)
    cin = FilmParticle.positive("inlet_concentration_kg_m3",inlet_concentration_kg_m3;zero=true)
    tr = ExchangeColumn.transport(flow_fractions,relative_delays,receiving_mix)
    weights = Vector{Float64}(contact_fractions)
    n = length(tr.flow)
    length(weights)==n && all(isfinite,weights) && all(>=(0),weights) &&
        abs(sum(weights)-1)<=1e-12 || throw(ArgumentError("Contact fractions must be nonnegative and sum to one"))
    solid = (1-eps)*p.rho*A*L
    liquid = eps*A*L
    mass = solid/ns
    area = 3mass/(p.rho*p.R)
    V = (liquid/ns).*tr.flow.*tr.delays
    Q = u*A
    stream = Q.*tr.flow
    cap = p.qe(cin)
    active = p.kf>0 && p.Ds>0 && cap>0
    conductance = active ? area*p.kf.*weights : zeros(n)
    startup = FilmParticle.positive("startup_s",startup_s;zero=true)
    startup<=p.tend || throw(ArgumentError("Startup time must not exceed the observation window"))
    dt0 = FilmParticle.positive("startup_dt_s",startup_dt_s)
    phases = filter(x->x[1]>0,[(startup,dt0),(p.tend-startup,p.dt)])
    steps_by_phase = [ceil(Int,duration/dt) for (duration,dt) in phases]
    steps = sum(steps_by_phase)
    steps<=200000 || throw(ArgumentError("More than 200000 steps with the startup schedule"))

    c = zeros(n,ns)
    q = zeros(p.nr,ns)
    a = similar(q)
    incoming = zeros(n)
    free,cs,qs = (zeros(n) for _ in 1:3)
    time,outlet,liquid_history,solid_history,cin_history,cout_history,residual =
        (zeros(steps+1) for _ in 1:7)
    maxflux = maxparticle = 0.
    minc = maxc = minq = maxq = 0.
    maxroot = step = 0
    absorbed = released = 0.
    started = time_ns()
    for ((duration,_),steps_here) in zip(phases,steps_by_phase)
        h = duration/steps_here
        op = FilmParticle.operator(p.nr,p.R,p.Ds,p.kf,h)
        ks = p.rho*p.Ds/op.half_width
        liquid_denom = V .+ h.*stream
        denom = liquid_denom .+ h.*conductance
        storage = V./liquid_denom
        advection = h.*stream./liquid_denom
        film_response = h.*conductance./denom
        liquid_response = liquid_denom./denom
        beta = active ? (ks/p.kf).*denom./liquid_denom : zeros(n)
        phase_start = time[step+1]
        for phase_step in 1:steps_here
            step += 1
            if active
                copyto!(a,q)
                ldiv!(op.F,a)
            end
            fill!(incoming,cin)
            liquid_inventory = solid_inventory = 0.
            for j in 1:ns
                if j>1 && (j-1)%slices==0
                    mul!(incoming,tr.mix,view(c,:,j-1))
                end
                @inbounds for i in 1:n
                    free[i] = storage[i]*c[i,j]+advection[i]*incoming[i]
                end
                mean_surface = 0.
                old_mean = new_mean = 0.
                if active
                    mean_surface,iterations = particle_surface!(cs,qs,free,beta,weights,
                        q[end,j],a[end,j],op.b[end],cap,capacity_kg_kg,affinity_m3_kg,
                        linear_partition_m3_kg,p.qe,p.slope)
                    maxroot = max(maxroot,iterations)
                    @inbounds for r in 1:p.nr
                        old_mean += op.w[r]*q[r,j]
                        value = a[r,j]+op.b[r]*mean_surface
                        q[r,j] = value
                        minq,maxq = min(minq,value),max(maxq,value)
                        new_mean += op.w[r]*value
                    end
                end
                total_film = 0.
                @inbounds for i in 1:n
                    concentration = liquid_response[i]*free[i]
                    if active && weights[i]>0
                        concentration += film_response[i]*cs[i]
                        film = p.kf*(concentration-cs[i])
                        radial = ks*(qs[i]-q[end,j])
                        maxflux = max(maxflux,abs(film-radial))
                        total_film += weights[i]*film
                        transfer = h*area*weights[i]*film
                        absorbed += max(transfer,0.)
                        released += max(-transfer,0.)
                    end
                    c[i,j] = incoming[i] = concentration
                    minc,maxc = min(minc,concentration),max(maxc,concentration)
                    liquid_inventory += V[i]*concentration
                end
                solid_inventory += mass*new_mean
                maxparticle = max(maxparticle,abs((new_mean-old_mean)/h-3total_film/(p.rho*p.R)))
            end
            j = step+1
            time[j] = phase_start+phase_step*h
            outlet[j] = dot(tr.flow,incoming)
            liquid_history[j] = liquid_inventory
            solid_history[j] = solid_inventory
            cin_history[j] = cin_history[j-1]+h*Q*cin
            cout_history[j] = cout_history[j-1]+h*Q*outlet[j]
            residual[j] = liquid_inventory+solid_inventory-cin_history[j]+cout_history[j]
        end
    end
    time[end] = p.tend
    (;time_s=time,outlet_concentration_kg_m3=outlet,liquid_inventory_kg=liquid_history,
      adsorbed_inventory_kg=solid_history,cumulative_inlet_kg=cin_history,cumulative_outlet_kg=cout_history,
      mass_residual_kg=residual,normalized_mass_residual=maximum(abs,residual)/max(cin_history[end],1e-30),
      normalized_interface_flux_mismatch=maxflux/max(p.kf*cin,1e-30),
      max_particle_balance_error_kg_kg_s=maxparticle,concentration_range_kg_m3=[minc,maxc],
      loading_range_kg_kg=[minq,maxq],liquid_volume_m3=liquid,physical_adsorbent_mass_kg=solid,
      external_particle_area_m2=ns*area,contact_area_by_path_m2=ns*area.*weights,
      total_film_conductance_m3_s=ns*sum(conductance),
      cumulative_film_uptake_kg=absorbed,cumulative_film_release_kg=released,
      nominal_liquid_transit_s=liquid/Q,periods=count,subcells=slices,classes=n,radial_cells=p.nr,
      steps=step,max_surface_iterations=maxroot,total_s=(time_ns()-started)/1e9)
end

end
