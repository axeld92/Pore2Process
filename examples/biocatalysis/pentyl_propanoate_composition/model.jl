module CompositionKinetics

using LinearAlgebra: Diagonal, norm, svdvals
import ..PentylPropanoate
import ..CatalystLoading

export Parameters, Curve, make_curves, fit_equilibrium, equilibrium_at, equilibrium_conversion,
       static_terms, net_mass_rate, acid_hazard, activity_bounds, batch_progress,
       predict_curves, fit_models, parameter_values, fit_model, specification

const MA = PentylPropanoate.ACID_MOLAR_MASS_KG_MOL
const MB = 0.088148
const CREF = inv(MA + 2MB)          # mol acid/kg nominal dry liquid in a 2:1 feed
const H_SCALE = 0.01              # min^-1 scale for the nonnegative optimization coordinate
const MODELS = ("ratio_reference","reversible_shared","acid_analogue","acid","composition",
                "acid_thermal","acid_decay","acid_activation","acid_both","protected_decay","protected_both",
                "residual_decay","protected_residual_decay")

struct Parameters
    k2::Float64                    # mol/(kg immobilised preparation s), apparent reference
    k2_slope::Float64
    k1::Float64                    # distinct empirical feed-rate reference only
    k1_slope::Float64
    acid_b::Float64                # dimensionless at reference acid molality
    alcohol_b::Float64
    product_slope::Float64         # lumped E/W response, not a water activity
    acid_heat_slope::Float64
    hazard_ref_min_inv::Float64
    hazard_temperature_slope::Float64
    initial_deficit::Float64
    recovery_min::Float64
    acid_exponent::Float64
    hazard_exponent::Float64
    protection_exponent::Float64    # relative acid/alcohol contact proxy
    residual_active_fraction::Float64 # hypothesized stable activity after damage
end

parameter_values(p::Parameters) = NamedTuple{fieldnames(Parameters)}(
    Tuple(getfield(p,f) for f in fieldnames(Parameters)))
parameter_vector(p::Parameters) = [getfield(p,f) for f in fieldnames(Parameters)]

struct Curve
    assay::Int
    feed::Int
    temperature_C::Float64
    ratio::Float64
    loading_pct::Float64
    Kx::Float64
    time_min::Vector{Float64}
    observed::Vector{Float64}
    weight::Float64
end

function make_curves(ids, feeds, temperatures, ratios, loadings, equilibria, counts, times, observed;
                     equal_curve_weights=true)
    ids, feeds, counts = (Int.(collect(v)) for v in (ids,feeds,counts))
    T,R,w,K,t,x = (Float64.(collect(v)) for v in (temperatures,ratios,loadings,equilibria,times,observed))
    n=length(ids)
    n>0 && length(unique(ids))==n && all(v -> length(v)==n,(feeds,T,R,w,K,counts)) &&
        sum(counts)==length(t)==length(x) && all(>(0),counts) &&
        all(v->all(isfinite,v),(T,R,w,K,t,x)) && all(v->all(>(0),v),(R,w,K,t)) &&
        all(>(-273.15),T) && all(v->v in (1,2),feeds) && all(y->0<=y<=1,x) ||
        throw(ArgumentError("Use unique assays, physical positive inputs and finite positive-time conversions."))
    curves=Curve[]
    index=1
    for i in 1:n
        range=index:index+counts[i]-1
        all(>(0),diff(t[range])) || throw(ArgumentError("Batch times must strictly increase."))
        weight=equal_curve_weights ? inv(sqrt(counts[i])) : sqrt(n/length(t))
        push!(curves,Curve(ids[i],feeds[i],T[i],R[i],w[i],K[i],t[range],x[range],weight))
        index+=counts[i]
    end
    curves
end

temperature_coordinate(T) = PentylPropanoate.temperature_coordinate(T)

"""Fit a separate apparent Kx(T) from training equilibrium observations only."""
function fit_equilibrium(temperatures, constants)
    t,K = (Float64.(collect(v)) for v in (temperatures,constants))
    length(t)==length(K)>=2 && length(unique(t))>=2 && all(>(0),K) &&
        all(isfinite,K) && all(x->isfinite(x) && x > -273.15,t) ||
        throw(ArgumentError("At least two distinct temperatures and positive training Kx values are required."))
    A=hcat(ones(length(t)),temperature_coordinate.(t))
    c=A\log.(K)
    (; log_K_reference=c[1], inverse_temperature_slope=c[2],
       log_K_rmse=sqrt(sum(abs2,A*c-log.(K))/length(t)))
end
equilibrium_at(log_K,slope,T)=exp(log_K+slope*temperature_coordinate(T))
equilibrium_conversion(R,K)=2R/(1+R+sqrt((R-1)^2+4R/K))

struct StaticTerms
    k::Float64
    Kx::Float64
    ratio::Float64
    mass_ratio::Float64             # c_A0/c_ref (stoichiometry and conserved mass)
    acid_b::Float64
    alcohol_b::Float64
    product_slope::Float64
    acid_exponent::Float64
    normalization::Float64
    equilibrium::Float64
    hazard_at_ref_min_inv::Float64
    hazard_exponent::Float64
    protection_exponent::Float64
    residual_active_fraction::Float64
    deficit::Float64
    tau_min::Float64
end

function static_terms(p::Parameters,T,R,K;feed=2)
    all(isfinite,(T,R,K)) && T > -273.15 && R>0 && K>0 && feed in (1,2) ||
        throw(ArgumentError("Physical temperature, feed and equilibrium required."))
    th=temperature_coordinate(T)
    k=feed==1 ? p.k1*exp(p.k1_slope*th) : p.k2*exp(p.k2_slope*th)
    b=p.acid_b*exp(p.acid_heat_slope*th)
    StaticTerms(k,K,R,inv(MA+R*MB)/CREF,b,p.alcohol_b,p.product_slope,p.acid_exponent,
        1+b+p.alcohol_b,equilibrium_conversion(R,K),p.hazard_ref_min_inv*exp(p.hazard_temperature_slope*th),
        p.hazard_exponent,p.protection_exponent,p.residual_active_fraction,p.initial_deficit,p.recovery_min)
end

"""Signed net reaction rate per kg of fully active immobilised preparation in mol/(kg s)."""
@inline function net_mass_rate(s::StaticTerms,x)
    # Local mole fractions at fixed original liquid mass and total moles.
    a=s.mass_ratio*(1-x)
    b=s.mass_ratio*(s.ratio-x)/2
    product=s.mass_ratio*x
    D=1+s.alcohol_b*b*b+(s.acid_b==0 ? 0. : s.acid_b*max(a,0.)^s.acid_exponent)
    factor=s.normalization/D
    if s.product_slope!=0
        factor*=exp(s.product_slope*product)
    end
    s.k*factor*((1-x)*(s.ratio-x)-x*x/s.Kx)/(1+s.ratio)^2
end

@inline function acid_hazard(s::StaticTerms,x)
    s.hazard_at_ref_min_inv==0 && return 0.
    base=s.hazard_at_ref_min_inv*(max(0.,s.mass_ratio*(1-x)))^s.hazard_exponent
    if s.protection_exponent!=0
        # At reference 2:1 dry feed, 2*A/B=1. Excess alcohol may
        # reduce acid-mediated damage more than acid molality alone.
        # This is a clock-law scenario, not the source's fitted Kd.
        base*=(2max(0.,1-x)/max(s.ratio-x,1e-12))^s.protection_exponent
    end
    base
end
@inline recovery(s::StaticTerms,t)=s.deficit==0 ? 1. : 1-s.deficit*exp(-t/s.tau_min)
@inline retained(s::StaticTerms,H)=s.residual_active_fraction+(1-s.residual_active_fraction)*exp(-max(0.,H))
@inline active_fraction(s::StaticTerms,t,H)=recovery(s,t)*retained(s,H)

"""Conditional activity envelope for stationary beads at elapsed age, not fluid age."""
function activity_bounds(s::StaticTerms,age_min)
    isfinite(age_min) && age_min>=0 || throw(ArgumentError("Use a finite nonnegative catalyst age."))
    max_h=acid_hazard(s,0.)
    min_h=acid_hazard(s,s.equilibrium)
    g=recovery(s,age_min)
    (; inlet_hazard_min_inv=max_h,equilibrium_hazard_min_inv=min_h,recovery_fraction=g,
       lower=g*retained(s,max_h*age_min),upper=g*retained(s,min_h*age_min))
end

@inline function rhs(s,exposure_rate,t,x,H)
    # Stage extrapolations must not invent negative concentrations. The
    # accepted state stays below the thermodynamic root; no clipping of data.
    state=clamp(x,0.,s.equilibrium)
    (exposure_rate*max(0.,net_mass_rate(s,state))*active_fraction(s,t,H),
     acid_hazard(s,state))
end

@inline function rk4(s,E,t,x,h,dt)
    ax,ah=rhs(s,E,t,x,h)
    bx,bh=rhs(s,E,t+dt/2,x+dt*ax/2,h+dt*ah/2)
    cx,ch=rhs(s,E,t+dt/2,x+dt*bx/2,h+dt*bh/2)
    dx,dh=rhs(s,E,t+dt,x+dt*cx,h+dt*ch)
    (x+dt*(ax+2bx+2cx+dx)/6,h+dt*(ah+2bh+2ch+dh)/6)
end

"""Progress and catalyst state at the requested elapsed minutes.

Step-doubling controls error in both X and integrated hazard. The unclipped
physical root is checked on each accepted state. The proportional,
clock-free reversible limit uses the previous independently tested formula.
"""
function batch_progress(p::Parameters,T,R,w,K,times;feed=2,rtol=2e-7,atol=2e-9,max_step_min=20.)
    t=Float64.(collect(times))
    all(isfinite,t) && all(>=(0),t) && all(>(0),diff(t)) &&
        isfinite(w) && w>=0 && rtol>0 && atol>0 && max_step_min>0 ||
        throw(ArgumentError("Use ordered nonnegative times, physical mass and integration tolerances."))
    s=static_terms(p,T,R,K;feed)
    E=60MA*w/100
    result=zeros(length(t)); activity=ones(length(t)); hazards=zeros(length(t))
    # Conditioning-only integrates exactly into the reaction exposure for
    # the un-inhibited reversible case; otherwise integrate the same ODE.
    if s.acid_b==s.alcohol_b==s.product_slope==s.hazard_at_ref_min_inv==0
        for (i,ti) in enumerate(t)
            effective_time=ti+s.deficit*s.tau_min*expm1(-ti/s.tau_min)
            result[i]=CatalystLoading.progress_conversion(s.k,s.Kx,s.ratio,E*max(0.,effective_time))
            activity[i]=recovery(s,ti)
        end
        return (; conversion=result,activity,hazards,accepted_steps=0,rejected_steps=0,
                  equilibrium=s.equilibrium)
    end
    x,H,clock=0.,0.,0.
    dt=min(1.,max_step_min)
    accepted,rejected=0,0
    for (i,goal) in enumerate(t)
        while clock<goal
            remaining=goal-clock
            if remaining <= 4eps(max(goal,1.))
                clock=goal
                break
            end
            step=min(dt,remaining,max_step_min)
            xc,hc=rk4(s,E,clock,x,H,step)
            xm,hm=rk4(s,E,clock,x,H,step/2)
            xf,hf=rk4(s,E,clock+step/2,xm,hm,step/2)
            err=max(abs(xf-xc)/(15*(atol+rtol*max(abs(x),abs(xf)))),
                    abs(hf-hc)/(15*(atol+rtol*max(abs(H),abs(hf)))))
            all(isfinite,(xf,hf,err)) || error("Nonfinite reaction integration.")
            physical=xf>=x-10atol && xf<=s.equilibrium+10atol && hf>=H-10atol
            if err<=1 && physical
                x=clamp(xf,x,s.equilibrium);H=max(H,hf)
                clock=step==remaining ? goal : clock+step
                accepted+=1
                dt=min(max_step_min,step*clamp(err==0 ? 3. : .9err^(-.2),.3,3.))
            else
                rejected+=1
                dt=step*(physical ? clamp(.8max(err,1.)^(-.2),.1,.5) : .5)
            end
            accepted+rejected<=100_000 && dt>1e-12 ||
                error("Reaction integrator exceeded its step limit.")
        end
        result[i]=x;hazards[i]=H;activity[i]=active_fraction(s,goal,H)
    end
    (; conversion=result,activity,hazards,accepted_steps=accepted,rejected_steps=rejected,
      equilibrium=s.equilibrium)
end

function predict_curves(p,curves;rtol=2e-7,atol=2e-9,max_step_min=20.)
    predictions=[batch_progress(p,c.temperature_C,c.ratio,c.loading_pct,c.Kx,c.time_min;
        feed=c.feed,rtol,atol,max_step_min) for c in curves]
    (; conversion=vcat([r.conversion for r in predictions]...),
       activity=vcat([r.activity for r in predictions]...),hazards=vcat([r.hazards for r in predictions]...),
       max_accepted_steps=maximum(r.accepted_steps for r in predictions),
       max_rejected_steps=maximum(r.rejected_steps for r in predictions))
end

# Optimization coordinates preserve exact nested zero terms and are scaled
# to [0,1] for finite differences. Positive constants use log transforms.
function specification(name;fixed_protection=NaN)
    name in MODELS || throw(ArgumentError("Unknown model $name"))
    names=["log_k2","k2_slope"]
    lower=[log(1e-4),name=="ratio_reference" ? -4. : 0.]
    upper=[log(10.),4.]
    initial=[log(.14),.35]
    function add(n,lo,hi,value)
        push!(names,n);push!(lower,lo);push!(upper,hi);push!(initial,value)
    end
    if name=="ratio_reference"
        add("log_k1",log(1e-4),log(10.),log(.13));add("k1_slope",-4.,4.,-.4)
    elseif name in ("acid","composition","acid_thermal","acid_decay","acid_activation","acid_both",
                    "protected_decay","protected_both","residual_decay","protected_residual_decay")
        add("log1p_acid_b",0.,log1p(100.),log1p(1.))
    end
    if name=="composition"
        add("log1p_alcohol_b",0.,log1p(100.),log1p(.4))
        add("product_slope",-4.,4.,0.)
    end
    name=="acid_thermal" && add("acid_heat_slope",0.,8.,2.)
    if name in ("acid_decay","acid_both","protected_decay","protected_both",
                "residual_decay","protected_residual_decay")
        add("log1p_hazard_ref_scaled",0.,log1p(100.),log1p(.25))
        add("hazard_temperature_slope",0.,10.,4.)
    end
    if name in ("acid_activation","acid_both","protected_both")
        add("initial_deficit",0.,1.,.35)
        add("log_recovery_min",log(.05),log(300.),log(8.))
    end
    if name in ("protected_decay","protected_both","protected_residual_decay") && !isfinite(fixed_protection)
        add("protection_exponent",0.,10.,4.)
    end
    if name in ("residual_decay","protected_residual_decay")
        add("residual_active_fraction",0.,1.,.08)
    end
    (; names,lower,upper,initial)
end

function decode(name,spec,values,n,m,analogue_KI,fixed_protection)
    v=Dict(zip(spec.names,values))
    k2=exp(v["log_k2"]); slope=v["k2_slope"]
    ba=name=="acid_analogue" ? (CREF/analogue_KI)^n : expm1(get(v,"log1p_acid_b",0.))
    Parameters(k2,slope,exp(get(v,"log_k1",v["log_k2"])),get(v,"k1_slope",slope),
        ba,expm1(get(v,"log1p_alcohol_b",0.)),get(v,"product_slope",0.),get(v,"acid_heat_slope",0.),
        H_SCALE*expm1(get(v,"log1p_hazard_ref_scaled",0.)),get(v,"hazard_temperature_slope",0.),
        get(v,"initial_deficit",0.),exp(get(v,"log_recovery_min",0.)),n,m,
        name in ("protected_decay","protected_both","protected_residual_decay") ?
        (isfinite(fixed_protection) ? fixed_protection : get(v,"protection_exponent",0.)) : 0.,
        get(v,"residual_active_fraction",0.))
end

function jacobian(residual,q;shift=2e-4)
    J=zeros(length(residual(q)),length(q))
    for j in eachindex(q)
        a=copy(q);b=copy(q)
        a[j]=max(0.,q[j]-shift);b[j]=min(1.,q[j]+shift)
        J[:,j]=(residual(b)-residual(a))/(b[j]-a[j])
    end
    J
end

function projected_gradient(q,g)
    [((q[j]<=1e-9 && g[j]>0) || (q[j]>=1-1e-9 && g[j]<0)) ? 0. : g[j]
     for j in eachindex(q)]
end

"""Deterministic projected damped least squares in bounded unit coordinates."""
function bounded_fit(residual,start;max_iterations=90)
    q=clamp.(copy(start),0.,1.)
    r=residual(q);sse=sum(abs2,r);damping=1e-3
    iterations=0
    for iteration in 1:max_iterations
        iterations=iteration
        J=jacobian(residual,q)
        g=J'*r
        pg=projected_gradient(q,g)
        norm(pg,Inf)<=2e-5 && break
        free=findall(!=(0.),pg)
        isempty(free) && break
        H=J[:,free]'*J[:,free]
        diagH=max.([H[i,i] for i in axes(H,1)],1e-4)
        accepted=false
        for _ in 1:16
            step=-(H+damping*Diagonal(diagH))\g[free]
            trial=copy(q)
            trial[free]=clamp.(q[free]+step,0.,1.)
            newr=residual(trial);newcost=sum(abs2,newr)
            if newcost<sse-1e-14
                q,r,sse=trial,newr,newcost
                damping=max(damping*.3,1e-10);accepted=true;break
            end
            damping*=5
        end
        if !accepted
            # Projection can interfere with the Gauss-Newton direction near
            # active bounds. Try a scaled descent step before declaring a stall.
            direction=-pg/max(norm(pg,Inf),1e-20)
            for j in 0:22
                trial=clamp.(q+(.15*0.5^j)*direction,0.,1.)
                newr=residual(trial);newcost=sum(abs2,newr)
                if newcost<sse-1e-14
                    q,r,sse=trial,newr,newcost;accepted=true;damping=1e-2;break
                end
            end
        end
        !accepted && break
    end
    J=jacobian(residual,q)
    g=J'*r
    pg=norm(projected_gradient(q,g),Inf)
    sv=svdvals(J)
    numerical_rank=count(>(maximum(sv)*1e-5),sv)
    (; q,sse,iterations,projected_gradient_inf=pg,numerical_rank,
       singular_values=sv,converged=pg<=1e-4)
end

function seed_from(spec,fit;positive_new_terms=false)
    old=Dict(zip(fit.names,fit.coordinates))
    [haskey(old,n) ? old[n] :
      (positive_new_terms ? spec.initial[i] :
        (startswith(n,"log1p_") || n in ("initial_deficit","acid_heat_slope","product_slope",
                                         "protection_exponent","residual_active_fraction") ? 0. : spec.initial[i]))
     for (i,n) in enumerate(spec.names)]
end

function fit_model(name,curves;starts_from=[],acid_exponent=3.5,hazard_exponent=3.5,analogue_KI=1.70,
                   fixed_protection=NaN,
                   rtol=2e-7,atol=2e-9,max_step_min=20.,max_iterations=90)
    (!isfinite(fixed_protection) || 0<=fixed_protection<=10) ||
        throw(ArgumentError("Protection exponent must be in [0,10]."))
    spec=specification(name;fixed_protection)
    observed=vcat([c.observed for c in curves]...)
    weights=vcat([fill(c.weight,length(c.observed)) for c in curves]...)
    decode_q(q)=decode(name,spec,spec.lower+q.*(spec.upper-spec.lower),
                       acid_exponent,hazard_exponent,analogue_KI,fixed_protection)
    residual(q)=weights.*(predict_curves(decode_q(q),curves;rtol,atol,max_step_min).conversion-observed)
    seeds=[spec.initial]
    for f in starts_from
        push!(seeds,seed_from(spec,f))
        push!(seeds,seed_from(spec,f;positive_new_terms=true))
    end
    if isempty(starts_from)
        alt=copy(spec.initial)
        alt[1]=log(.08);alt[2]=.8
        push!(seeds,alt)
    end
    # Unique starting coordinates, all derived solely from training fits
    # or predeclared physical-scale initial guesses.
    seeds=unique(seeds)
    runs=[bounded_fit(residual,(seed-spec.lower)./(spec.upper-spec.lower);max_iterations)
          for seed in seeds]
    best=runs[argmin([r.sse for r in runs])]
    raw=spec.lower+best.q.*(spec.upper-spec.lower)
    p=decode_q(best.q)
    active_bounds=[best.q[i]<1e-6 ? "lower" : best.q[i]>1-1e-6 ? "upper" : "interior"
                   for i in eachindex(best.q)]
    output=predict_curves(p,curves;rtol,atol,max_step_min)
    (; name,names=spec.names,coordinates=raw,parameters=p,weighted_sse=best.sse,
       best.iterations,best.projected_gradient_inf,best.numerical_rank,best.singular_values,
       best.converged,active_bounds,start_sse=[r.sse for r in runs],
       max_accepted_steps=output.max_accepted_steps,max_rejected_steps=output.max_rejected_steps)
end

"""Fit candidate family using only training curves and nested training-only starts."""
function fit_models(curves;names=collect(MODELS),acid_exponent=3.5,hazard_exponent=3.5,
                    analogue_KI=1.70,fixed_protection=NaN,rtol=2e-7,atol=2e-9,max_step_min=20.,max_iterations=90)
    results=Dict{String,Any}()
    for name in names
        name=String(name)
        parent_names = name=="acid" || name=="acid_analogue" ? ("reversible_shared",) :
            name in ("composition","acid_thermal","acid_decay","acid_activation") ? ("acid",) :
            name=="acid_both" ? ("acid_decay","acid_activation") :
            name=="protected_decay" ? ("acid_decay",) :
            name=="protected_both" ? ("protected_decay","acid_both") :
            name=="residual_decay" ? ("acid_decay",) :
            name=="protected_residual_decay" ? ("protected_decay","residual_decay") : ()
        parents=[results[p] for p in parent_names if haskey(results,p)]
        results[name]=fit_model(name,curves;starts_from=parents,acid_exponent,hazard_exponent,
            analogue_KI,fixed_protection,rtol,atol,max_step_min,max_iterations)
    end
    [results[String(name)] for name in names]
end

end
