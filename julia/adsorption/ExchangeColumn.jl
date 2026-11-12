module ExchangeColumn

using LinearAlgebra
using ..FilmParticle
export solve_column

# Eliminate the liquid concentration from the implicit film and local
# inventory equations. The remaining surface root is monotone.
@inline function surface_root(rhs, factor, guess, qe, slope)
    rhs >= -1e-13 || error("Negative surface root bracket")
    rhs = max(rhs, 0.)
    factor==0 && return rhs
    lo,hi = 0.,rhs
    x = clamp(guess,lo,hi)
    for _ in 1:50
        residual = x + factor*qe(x) - rhs
        abs(residual)<=1e-13*max(rhs,1e-6) && return x
        if residual>0
            hi = x
        else
            lo = x
        end
        next = x - residual/(1 + factor*slope(x))
        x = lo<next<hi ? next : (lo+hi)/2
    end
    error("Exchange-column surface root failed to converge")
end

function transport(flow_fractions, relative_delays, receiving_mix)
    flow = Vector{Float64}(flow_fractions)
    delays = Vector{Float64}(relative_delays)
    mix = Matrix{Float64}(receiving_mix)
    n = length(flow)
    n>=1 && length(delays)==n && size(mix)==(n,n) || throw(ArgumentError("Inconsistent path dimensions"))
    all(isfinite,flow) && all(isfinite,delays) && all(isfinite,mix) &&
        all(>(0),flow) && all(>(0),delays) && all(>=(0),mix) || throw(ArgumentError("Invalid path flow, time or mixing"))
    abs(sum(flow)-1)<=1e-12 && abs(dot(flow,delays)-1)<=1e-12 ||
        throw(ArgumentError("Path fractions or mean relative transit do not sum to one"))
    maximum(abs,vec(sum(mix;dims=2)).-1)<=1e-12 &&
        maximum(abs,transpose(mix)*flow-flow)<=1e-12 ||
        throw(ArgumentError("Mixing must preserve uniform concentration and path throughputs"))
    (;flow,delays,mix)
end

"""
Forward path classes within each axial period. Adjacent periods exchange
flow through the receiving-class mixing matrix. Liquid and solid inventory
are fixed. The local particle operator and film/isotherm are those in
FilmParticle, and its local flux is shared with the liquid balance.
"""
function solve_column(flow_fractions,relative_delays,receiving_mix;
        periods,subcells=1,length_m,area_m2,porosity,superficial_velocity_m_s,
        inlet_concentration_kg_m3,startup_s=0.,startup_dt_s=.05,kwargs...)
    p = FilmParticle.parameters(;kwargs...)
    p.qinit==0 || throw(ArgumentError("Exchange column requires initially clean particles"))
    section_count = FilmParticle.count("periods",periods,1,2000)
    slices = FilmParticle.count("subcells",subcells,1,64)
    L = FilmParticle.positive("length_m",length_m)
    A = FilmParticle.positive("area_m2",area_m2)
    eps = FilmParticle.positive("porosity",porosity)
    eps<1 || throw(ArgumentError("porosity must be < 1"))
    u = FilmParticle.positive("superficial_velocity_m_s",superficial_velocity_m_s)
    cin = FilmParticle.positive("inlet_concentration_kg_m3",inlet_concentration_kg_m3;zero=true)
    tr = transport(flow_fractions,relative_delays,receiving_mix)
    n = length(tr.flow)
    nv = n*slices*section_count
    solid = (1-eps)*p.rho*A*L
    liquid = eps*A*L
    share = tr.flow.*tr.delays/(slices*section_count)
    V = liquid.*share
    M = solid.*share
    Q = u*A
    stream = Q.*tr.flow
    startup = FilmParticle.positive("startup_s",startup_s;zero=true)
    startup<=p.tend || throw(ArgumentError("Startup time must not exceed the observation window"))
    initial_dt = FilmParticle.positive("startup_dt_s",startup_dt_s)
    phases = filter(phase->phase[1]>0,[(startup,initial_dt),(p.tend-startup,p.dt)])
    counts = [ceil(Int,duration/dt) for (duration,dt) in phases]
    steps = sum(counts)
    steps<=200000 || throw(ArgumentError("More than 200000 steps with the startup schedule"))

    c = zeros(nv)
    q = zeros(p.nr,nv)
    a = similar(q)
    meanq = zeros(nv)
    prior_out = zeros(n)
    incoming = similar(prior_out)
    time = zeros(steps+1)
    outlet = zeros(steps+1)
    cin_mass = zeros(steps+1)
    cout_mass = zeros(steps+1)
    liquid_mass = zeros(steps+1)
    solid_mass = zeros(steps+1)
    residual = zeros(steps+1)
    max_flux_mismatch = 0.
    max_particle_mismatch = 0.
    concentration_min,concentration_max,loading_min,loading_max = 0.,0.,0.,0.
    started = time_ns()
    step = 0
    for ((duration,_), count) in zip(phases,counts)
        h = duration/count
        op = FilmParticle.operator(p.nr,p.R,p.Ds,p.kf,h)
        theta = p.kf==0 || p.Ds==0 ? 0. : p.rho*p.Ds/(op.half_width*p.kf)
        denominator = V+h*stream
        storage = V./denominator
        advection = h*stream./denominator
        solid_ratio = M./denominator
        factor = theta*op.last_factor .+ solid_ratio.*op.mean_b
        phase_start = time[step+1]
        for phase_step in 1:count
            step += 1
            copyto!(a,q)
            ldiv!(op.F,a)
            fill!(incoming,cin)
            liquid_inventory,solid_inventory = 0.,0.
            for period in 1:section_count
                period>1 && mul!(incoming,tr.mix,prior_out)
                @inbounds for path in 1:n
                    upstream = incoming[path]
                    for slice in 1:slices
                        index = slice+slices*(path-1+n*(period-1))
                        mean_a = 0.
                        for r in 1:p.nr
                            mean_a += op.w[r]*a[r,index]
                        end
                        oldmean = meanq[index]
                        rhs = storage[path]*c[index] + advection[path]*upstream +
                              solid_ratio[path]*(oldmean-mean_a) + theta*a[end,index]
                        cs = surface_root(rhs,factor[path],c[index],p.qe,p.slope)
                        qs = p.qe(cs)
                        concentration = cs + theta*(op.last_factor*qs-a[end,index])
                        concentration >= -1e-12*max(cin,1e-6) || error("Negative liquid concentration")
                        concentration = max(concentration,0.)
                        c[index] = concentration
                        upstream = concentration
                        meanq[index] = mean_a + op.mean_b*qs
                        for r in 1:p.nr
                            loading = a[r,index]+op.b[r]*qs
                            q[r,index] = loading
                            loading_min = min(loading_min,loading)
                            loading_max = max(loading_max,loading)
                        end
                        film = p.kf*(concentration-cs)
                        radial = theta==0 ? 0. : p.rho*p.Ds/op.half_width*(qs-q[end,index])
                        max_flux_mismatch = max(max_flux_mismatch,abs(film-radial))
                        max_particle_mismatch = max(max_particle_mismatch,
                            abs((meanq[index]-oldmean)/h - 3film/(p.rho*p.R)))
                        liquid_inventory += V[path]*concentration
                        solid_inventory += M[path]*meanq[index]
                        concentration_min = min(concentration_min,concentration)
                        concentration_max = max(concentration_max,concentration)
                    end
                    prior_out[path] = upstream
                end
            end
            j = step+1
            time[j] = phase_start + phase_step*h
            outlet[j] = dot(tr.flow,prior_out)
            cin_mass[j] = cin_mass[j-1]+h*Q*cin
            cout_mass[j] = cout_mass[j-1]+h*Q*outlet[j]
            liquid_mass[j] = liquid_inventory
            solid_mass[j] = solid_inventory
            residual[j] = liquid_inventory+solid_inventory-cin_mass[j]+cout_mass[j]
        end
    end
    time[end] = p.tend
    (;time_s=time,outlet_concentration_kg_m3=outlet,liquid_inventory_kg=liquid_mass,
      adsorbed_inventory_kg=solid_mass,cumulative_inlet_kg=cin_mass,cumulative_outlet_kg=cout_mass,
      mass_residual_kg=residual,normalized_mass_residual=maximum(abs,residual)/max(cin_mass[end],1e-30),
      normalized_interface_flux_mismatch=max_flux_mismatch/max(p.kf*cin,1e-30),
      max_particle_balance_error_kg_kg_s=max_particle_mismatch,
      concentration_range_kg_m3=[concentration_min,concentration_max],
      loading_range_kg_kg=[loading_min,loading_max],
      liquid_volume_m3=liquid,physical_adsorbent_mass_kg=solid,
      nominal_liquid_transit_s=liquid/Q,periods=section_count,subcells=slices,
      classes=n,radial_cells=p.nr,steps,total_s=(time_ns()-started)/1e9)
end

end
