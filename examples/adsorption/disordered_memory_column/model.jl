module MemoryColumn

using LinearAlgebra
using ..FilmParticle
using ..SharedContactColumn: particle_surface!

struct Delays
    donor::Vector{Int}
    receiver::Vector{Int}
    share::Vector{Float64}
    seconds::Vector{Float64}
end

function delays(T, first, second, donor_flow)
    result = Delays(Int[],Int[],Float64[],Float64[])
    for ij in CartesianIndices(T)
        p = T[ij]
        p == 0 && continue
        i,j = Tuple(ij)
        m = first[ij]/p
        v = second[ij]/p-m^2
        m > 0 && v >= -1e-12*m^2 || error("Invalid conditional transit-time moments")
        weight = donor_flow[i]*p
        if v <= 4eps(m^2)
            push!(result.donor,i); push!(result.receiver,j)
            push!(result.share,weight); push!(result.seconds,m)
        else
            offset = min(m/2,sqrt(v))
            high = offset^2/(v+offset^2)
            append!(result.donor,(i,i)); append!(result.receiver,(j,j))
            append!(result.share,(weight*(1-high),weight*high))
            append!(result.seconds,(m-offset,m+v/offset))
        end
    end
    isapprox(sum(result.share),sum(donor_flow);rtol=1e-12,atol=1e-15) &&
        isapprox(dot(result.share,result.seconds),dot(donor_flow,vec(sum(first;dims=2)));
                 rtol=1e-12,atol=0.) &&
        isapprox(dot(result.share,result.seconds.^2),dot(donor_flow,vec(sum(second;dims=2)));
                 rtol=1e-12,atol=0.) || error("Delay quadrature failed the supplied moments")
    result
end

function schedule(tend,dt,startup,early)
    0 <= startup <= tend && dt > 0 && early > 0 || throw(ArgumentError("Invalid time schedule"))
    time = [0.]
    for (duration,h) in ((startup,early),(tend-startup,dt))
        duration == 0 && continue
        count = ceil(Int,duration/h)
        count + length(time) <= 200001 || throw(ArgumentError("Too many time steps"))
        append!(time,range(last(time),last(time)+duration;length=count+1)[2:end])
    end
    time
end

function prepare_intervals!(left,right,alpha,delay,time,next,slot,scale)
    @inbounds for atom in eachindex(delay)
        target = max(0.,time[next]-scale*delay[atom])
        k = searchsortedlast(time,target,1,next-1,Base.Order.Forward)
        left[atom],right[atom] = slot[k],slot[k+1]
        alpha[atom] = (target-time[k])/(time[k+1]-time[k])
    end
end

"""Conserve each delayed stream by differencing integrated inlet histories."""
function advect!(bulk,history,delivered,rule,phase,current,left,right,alpha,h,flow)
    fill!(bulk,0.)
    volume_concentration = 0.
    @inbounds for a in eachindex(rule.share)
        i,j = rule.donor[a],rule.receiver[a]
        k,l,f = left[a],right[a],alpha[a]
        x = history[i,k,phase]
        lagged = muladd(f,history[i,l,phase]-x,x)
        bulk[j] += rule.share[a]*(lagged-delivered[a])
        delivered[a] = lagged
        volume_concentration += rule.share[a]*(history[i,current,phase]-lagged)
    end
    @inbounds for j in eachindex(flow)
        bulk[j] /= h*flow[j]
    end
    volume_concentration
end

function solve_column(rules,flow,contact; blocks, liquid_volume_m3, dry_mass_kg,
        flow_m3_s, inlet_concentration_kg_m3, startup_s=0., startup_dt_s=.01,
        capacity_kg_kg=.1115,affinity_m3_kg=139.,linear_partition_m3_kg=nothing,
        seam_history=false,kwargs...)
    p = FilmParticle.parameters(;capacity_kg_kg,affinity_m3_kg,linear_partition_m3_kg,kwargs...)
    p.qinit==0 || throw(ArgumentError("The column requires initially clean carbon"))
    nsample,n = length(rules),length(flow)
    size(contact)==(n,nsample) && all(>=(0),contact) &&
        all(x->isapprox(x,1.;atol=1e-13,rtol=0.),sum(contact;dims=1)) ||
        throw(ArgumentError("Expected normalized nonnegative contact shares"))
    all(>(0),flow) && isapprox(sum(flow),1.;atol=1e-13,rtol=0.) ||
        throw(ArgumentError("Expected normalized positive state flow shares"))
    all(k->isfinite(k) && k>0,(flow_m3_s,liquid_volume_m3,dry_mass_kg)) && blocks>=1 ||
        throw(ArgumentError("Column inventories, flow and pattern count must be positive"))
    Q, cin = flow_m3_s, inlet_concentration_kg_m3
    tau = sum(dot(r.share,r.seconds) for r in rules)
    scale = liquid_volume_m3/(Q*blocks*tau)
    max_delay = scale*maximum(maximum(r.seconds) for r in rules)
    for (j,rule) in enumerate(rules)
        got, sent = zeros(n), zeros(n)
        for a in eachindex(rule.share)
            got[rule.receiver[a]] += rule.share[a]
            sent[rule.donor[a]] += rule.share[a]
        end
        target = j==1 ? [1.;zeros(n-1)] : flow
        isapprox(got,flow;rtol=1e-12,atol=1e-15) &&
            isapprox(sent,target;rtol=1e-12,atol=1e-15) || error("Delay network does not conserve water")
    end
    time = schedule(p.tend,p.dt,startup_s,startup_dt_s)
    nstep = length(time)-1
    window = min(length(time),ceil(Int,max_delay/minimum(diff(time)))+3)
    slot = [mod1(k,window) for k in eachindex(time)]
    ns = nsample*blocks
    kinds = repeat(1:nsample,blocks)
    mass = dry_mass_kg/ns
    area = 3mass/(p.rho*p.R)
    cap = p.qe(cin)
    active = p.kf>0 && p.Ds>0 && cap>0
    weights = Matrix{Float64}(contact)
    particle = zeros(p.nr,ns)
    interior = similar(particle)
    c,cs,qs,rhs,incoming = (zeros(n) for _ in 1:5)
    history = zeros(n,window,ns)
    delivered = [zeros(length(rules[k].share)) for k in kinds]
    left = [zeros(Int,length(r.share)) for r in rules]
    right,alpha = deepcopy(left),[zeros(length(r.share)) for r in rules]
    out,passed,liquid,solid,residual = (zeros(length(time)) for _ in 1:5)
    seams = seam_history ? zeros(ns,length(time)) : zeros(0,0)
    maxflux,uptake,release = 0.,0.,0.
    minc,maxc,minq,maxq = 0.,0.,0.,0.
    old_h = -1.
    op = FilmParticle.operator(p.nr,p.R,p.Ds,p.kf,time[2])
    ks = p.rho*p.Ds/op.half_width
    ntu = area*p.kf.*weights./(Q.*flow)
    beta = active ? (ks/p.kf).*(1 .+ ntu) : zeros(n,nsample)
    t0 = time_ns()
    for step in 1:nstep
        current,prev = slot[step+1],slot[step]
        h = time[step+1]-time[step]
        if active
            if !isapprox(h,old_h;rtol=0.,atol=8eps(p.tend))
                op = FilmParticle.operator(p.nr,p.R,p.Ds,p.kf,h)
                old_h = h
            end
            copyto!(interior,particle)
            ldiv!(op.F,interior)
        end
        for k in 1:nsample
            prepare_intervals!(left[k],right[k],alpha[k],rules[k].seconds,time,step+1,slot,scale)
        end
        incoming[1] = cin
        for s in 1:ns
            k = kinds[s]
            if k==1
                history[1,current,s] = muladd(h,incoming[1],history[1,prev,s])
            else
                @inbounds for i in 1:n
                    history[i,current,s] = muladd(h,incoming[i],history[i,prev,s])
                end
            end
            liquid[step+1] += Q*advect!(c,history,delivered[s],rules[k],s,current,
                left[k],right[k],alpha[k],h,flow)
            if active
                @inbounds for i in 1:n
                    rhs[i] = clamp(c[i],0.,cin)
                end
                mean_surface,_ = particle_surface!(cs,qs,rhs,view(beta,:,k),view(weights,:,k),
                    particle[end,s],interior[end,s],op.b[end],cap,capacity_kg_kg,affinity_m3_kg,
                    linear_partition_m3_kg,p.qe,p.slope)
                mean_loading = 0.
                @inbounds for r in 1:p.nr
                    q = interior[r,s]+op.b[r]*mean_surface
                    particle[r,s] = q
                    mean_loading += op.w[r]*q
                    minq,maxq = min(minq,q),max(maxq,q)
                end
                solid[step+1] += mass*mean_loading
                @inbounds for i in 1:n
                    if weights[i,k]>0
                        cout = (c[i]+ntu[i,k]*cs[i])/(1+ntu[i,k])
                        film = p.kf*(cout-cs[i])
                        radial = ks*(qs[i]-particle[end,s])
                        maxflux = max(maxflux,abs(film-radial))
                        transferred = h*area*weights[i,k]*film
                        uptake += max(transferred,0.)
                        release += max(-transferred,0.)
                        c[i] = cout
                    end
                end
            end
            @inbounds for i in 1:n
                minc,maxc = min(minc,c[i]),max(maxc,c[i])
                incoming[i] = c[i]
            end
            if seam_history
                seams[s,step+1] = dot(flow,c)
            end
            if k==nsample
                incoming[1] = dot(flow,c)
            end
        end
        out[step+1] = incoming[1]
        passed[step+1] = passed[step]+h*Q*out[step+1]
        residual[step+1] = liquid[step+1]+solid[step+1]+passed[step+1]-Q*cin*time[step+1]
    end
    (;time_s=time,outlet_concentration_kg_m3=out,cumulative_outlet_kg=passed,
      liquid_inventory_kg=liquid,adsorbed_inventory_kg=solid,mass_residual_kg=residual,
      normalized_mass_residual=maximum(abs,residual)/max(Q*cin*time[end],1e-30),
      normalized_interface_flux_mismatch=maxflux/max(p.kf*cin,1e-30),
      concentration_range_kg_m3=[minc,maxc],loading_range_kg_kg=[minq,maxq],
      liquid_volume_m3=Q*blocks*tau*scale,physical_adsorbent_mass_kg=ns*mass,
      external_particle_area_m2=ns*area,total_film_conductance_m3_s=ns*area*p.kf,
      cumulative_film_uptake_kg=uptake,cumulative_film_release_kg=release,
      delay_scale=scale,maximum_delay_s=max_delay,sections=ns,seam_concentration=seams,
      wall_s=(time_ns()-t0)/1e9)
end

end
