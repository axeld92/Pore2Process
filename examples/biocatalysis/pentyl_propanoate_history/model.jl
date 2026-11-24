module SpatialCatalystHistory

using LinearAlgebra: Tridiagonal, norm, SingularException
import ..CompositionKinetics
import ..BeadFilm

const CK = CompositionKinetics
const GAMMA = 1 - inv(sqrt(2.0))
const STAGE_TOL = 5e-12

export simulate, film_slope_screen

"""Frozen local rate and stationary exposure, expressed per litre of liquid."""
struct LocalLaw
    state::CK.StaticTerms
    scale::Float64              # eta_I * .06 * C_cat [g/L]; multiplies mol/(kg s)
    acid0::Float64
    capacity::Float64           # mol/(L min) per unit difference in progress
    use_film::Bool
    surface_hazard::Bool
    recovered::Bool
end

@inline activity(l::LocalLaw,t,H) = l.recovered ? CK.retained(l.state,H) : CK.active_fraction(l.state,t,H)
@inline hazard(l::LocalLaw,x) = CK.acid_hazard(l.state,x)
@inline reaction(l::LocalLaw,a,x) = l.scale*a*CK.net_mass_rate(l.state,x)

@inline function derivative(fun,x,eq)
    lo=max(0.,x-1e-6); hi=min(eq,x+1e-6)
    (fun(hi)-fun(lo))/(hi-lo)
end

struct StepFailure <: Exception
    message::String
end
Base.showerror(io::IO, e::StepFailure) = print(io,e.message)

"""Evaluate a coupled implicit reaction/age stage at a given bulk progress.

At this stage H = H_base + alpha * h(X_local). The existing film root
computes surface progress. Differentiate that scalar balance implicitly for
the tridiagonal liquid Newton solve: if R(s,b) is the local reaction,
dR/db = capacity*(R_s+R_b)/(capacity-R_s). All derivatives use the
unchanged mass-specific and hazard functions, not a refitted rate law.
"""
@inline function film_result(l,bulk,surface,a,H,hh,source,dy,db)
    dy/l.capacity<.99 || throw(StepFailure("Local implicit film slope approaches its capacity."))
    d=l.capacity*(dy+db)/(l.capacity-dy)/l.acid0
    err=abs(source-l.capacity*(surface-bulk))
    (;surface,activity=a,H,hazard=hh,source=source/l.acid0,slope=d,
      film_error=err,positive_film_slope_ratio=max(0.,dy/l.capacity))
end

@inline function local_bare(l,alpha,Hbase,t,bulk)
    hh=hazard(l,bulk); H=Hbase+alpha*hh; a=activity(l,t,H)
    rf=x->reaction(l,activity(l,t,Hbase+alpha*hazard(l,x)),x)
    src=reaction(l,a,bulk)
    d=derivative(rf,bulk,l.state.equilibrium)/l.acid0
    (;surface=bulk,activity=a,H,hazard=hh,source=src/l.acid0,
      slope=d,film_error=0.,positive_film_slope_ratio=0.)
end

@inline function local_surface(l,alpha,Hbase,t,bulk)
    # H at the stage depends on the same surface composition as rate.
    rf=y->reaction(l,activity(l,t,Hbase+alpha*hazard(l,y)),y)
    surface=BeadFilm.surface_progress(rf,bulk,l.state.equilibrium,l.capacity)
    hh=hazard(l,surface); H=Hbase+alpha*hh; a=activity(l,t,H)
    source=reaction(l,a,surface)
    dy=derivative(rf,surface,l.state.equilibrium)
    film_result(l,bulk,surface,a,H,hh,source,dy,0.)
end

@inline function local_bulk(l,alpha,Hbase,t,bulk)
    # For the calibration-basis exposure rule H depends on bulk only.
    hh=hazard(l,bulk); H=Hbase+alpha*hh; a=activity(l,t,H)
    rf=y->reaction(l,a,y)
    surface=BeadFilm.surface_progress(rf,bulk,l.state.equilibrium,l.capacity)
    source=rf(surface)
    dy=derivative(rf,surface,l.state.equilibrium)
    a_at=x->activity(l,t,Hbase+alpha*hazard(l,x))
    db=alpha==0 || l.state.hazard_at_ref_min_inv==0 ? 0. :
        l.scale*CK.net_mass_rate(l.state,surface)*derivative(a_at,bulk,l.state.equilibrium)
    film_result(l,bulk,surface,a,H,hh,source,dy,db)
end

@noinline function local_state(l::LocalLaw,alpha,Hbase,t,bulk)
    !l.use_film ? local_bare(l,alpha,Hbase,t,bulk) :
        l.surface_hazard ? local_surface(l,alpha,Hbase,t,bulk) : local_bulk(l,alpha,Hbase,t,bulk)
end

"""Sample the film bracket monotonicity over the whole physical interval.

For surface-age coupling the stage reaction derivative is affine in the
surviving labile fraction at the start of the stage. Its largest positive
value over past damage/recovery occurs at a recovered, initially undamaged
or fully depleted labile population. Sample both limiting derivatives and
the stage-substep range, separately from the Newton iteration's local check.
"""
function film_slope_screen(p::CK.Parameters,T,feed,K;catalyst_g_L,eta_internal,capacity,
                           max_step_min=2.,surface_hazard=false)
    s=CK.static_terms(p,T,Float64(feed),K;feed)
    scale=eta_internal*.06*catalyst_g_L
    q=y->scale*CK.net_mass_rate(s,y)
    function sampled(nx,nt)
        # Includes fully active and stable-only limits. For a bulk hazard
        # the film source varies only through the stationary activity factor.
        x=range(0.,s.equilibrium;length=nx)
        .0+max(0.,maximum(derivative(q,y,s.equilibrium) for y in x)/capacity,
            surface_hazard ? maximum(derivative(z->CK.retained(s,alpha*CK.acid_hazard(s,z))*q(z),
                       y,s.equilibrium)/capacity for alpha in range(0.,GAMMA*max_step_min;length=nt),y in x) : 0.)
    end
    a=sampled(65,9); b=sampled(257,17)
    (;max_sampled_slope_ratio=max(a,b),refinement_change=abs(b-a))
end

"""Same exponentially fitted face flux as the existing steady bed solver."""
struct AxialGrid
    n::Int
    tau_min::Float64
    down::Float64
    inlet::Float64
    lower::Vector{Float64}
    diagonal::Vector{Float64}
    upper::Vector{Float64}
end

function AxialGrid(n::Int,tau,Pe;inlet_flux=false)
    n>=2 && isfinite(tau) && tau>0 && isfinite(Pe) && Pe>0 ||
        throw(ArgumentError("Require positive finite residence and Pe and at least two cells."))
    down=inv(expm1(Pe/n)); inlet=inlet_flux ? 0. : inv(expm1(Pe/(2n)))
    a=fill(1+2down,n); a[1]=1+down+inlet; a[end]=1+down
    AxialGrid(n,tau,down,inlet,fill(-(1+down),n-1),a,fill(-down,n-1))
end

function face_fluxes(g::AxialGrid,x)
    j=Vector{Float64}(undef,g.n+1)
    j[1]=-g.inlet*x[1]
    @inbounds for i in 1:g.n-1
        j[i+1]=x[i]+g.down*(x[i]-x[i+1])
    end
    j[end]=x[end]
    j
end

"""Evaluate physical stage residual and fluxes independently of the Newton step."""
function evaluate(g,l,b,Hbase,t,alpha,x)
    n=g.n
    f=zeros(n); H=zeros(n); h=zeros(n); a=zeros(n); surf=zeros(n)
    q=zeros(n); dq=zeros(n)
    j=face_fluxes(g,x)
    max_film=0.; max_slope=0.
    @inbounds for i in 1:n
        local_result=local_state(l,alpha,Hbase[i],t,x[i])
        surf[i]=local_result.surface; a[i]=local_result.activity
        H[i]=local_result.H; h[i]=local_result.hazard
        q[i]=local_result.source; dq[i]=local_result.slope
        f[i]=q[i]-(n/g.tau_min)*(j[i+1]-j[i])
        max_film=max(max_film,local_result.film_error)
        max_slope=max(max_slope,local_result.positive_film_slope_ratio)
    end
    r=x.-b.-alpha.*f
    valid=all(isfinite,r) && all(isfinite,H) && all(isfinite,dq) &&
        all(isfinite,q) && all(v->0<=v<=1,a) && all(>=(0.),h) &&
        all(v->-1e-12<=v<=l.state.equilibrium+1e-12,surf) &&
        all(v->v>=-1e-12,q)
    valid || throw(StepFailure("Nonphysical or nonfinite implicit stage."))
    (;x,H,h,a,surf,q,dq,f,j,r,generated=g.tau_min*sum(q)/n,
      max_film,max_slope)
end

@inline function stage_measure(r,alpha)
    # A small componentwise residual alone can accumulate systematically
    # over a long history. Budget the spatial mean residual per step too.
    max(norm(r,Inf)/STAGE_TOL,abs(sum(r)/length(r))/max(5e-15,2e-12*alpha))
end

function solve_stage(g,l,b,Hbase,t,alpha,guess)
    # No clipping of an accepted update. Clamped trial states are evaluated
    # by the original residual and accepted only on actual convergence.
    x=copy(guess); state=evaluate(g,l,b,Hbase,t,alpha,x)
    lambda=alpha*g.n/g.tau_min
    lower=lambda.*g.lower; upper=lambda.*g.upper
    diag_base=1 .+ lambda.*g.diagonal
    for k in 0:18
        res=stage_measure(state.r,alpha)
        res<=1. && return state,k
        k==18 && break
        step=Tridiagonal(lower,diag_base-alpha.*state.dq,upper)\state.r
        accepted=false
        for halve in 0:15
            trial=clamp.(state.x.-step./2.0^halve,0.,l.state.equilibrium)
            try
                candidate=evaluate(g,l,b,Hbase,t,alpha,trial)
                if stage_measure(candidate.r,alpha)<res
                    state=candidate; accepted=true; break
                end
            catch e
                e isa StepFailure || rethrow()
            end
        end
        accepted || throw(StepFailure("Implicit stage did not reduce the physical residual."))
    end
    throw(StepFailure("Implicit stage exceeded the Newton iteration limit."))
end

"""Two-stage SDIRK2; both local exposure and liquid progress are at stage time."""
function advance(g,l,x,H,t,dt)
    alpha=GAMMA*dt
    a,ni=solve_stage(g,l,x,H,t+alpha,alpha,x)
    bbase=x.+(1-GAMMA)*dt.*a.f
    Hbase=H.+(1-GAMMA)*dt.*a.h
    b,nj=solve_stage(g,l,bbase,Hbase,t+dt,alpha,a.x)
    all(b.H.>=H) || throw(StepFailure("Stationary cumulative exposure decreased."))
    # Stiffly accurate final stage. Physical source/face integrals use its
    # actual fluxes and the first stage, not a difference of inventories.
    integral(first,second)=dt*((1-GAMMA)*first+GAMMA*second)
    (; state=b,generated=integral(a.generated,b.generated),
      effluent=integral(a.j[end],b.j[end]),backdiffusion=integral(-a.j[1],-b.j[1]),
      iterations=max(ni,nj),stage_error=max(norm(a.r,Inf),norm(b.r,Inf)),
      max_film=max(a.max_film,b.max_film),max_slope=max(a.max_slope,b.max_slope))
end

mutable struct HistoryLedger
    t::Float64
    dt::Float64
    generated::Float64
    effluent::Float64
    backdiffusion::Float64
    accepted::Int
    rejected::Int
    max_iters::Int
    min_step::Float64
    max_step::Float64
    max_stage::Float64
    max_balance::Float64
    max_balance_mol::Float64
    max_bound::Float64
    max_film::Float64
    max_slope::Float64
end

@inline inventory_norm(g::AxialGrid,x)=g.tau_min*sum(x)/g.n

@noinline function accept_step!(led,g,result,step,newt,hmin,hmax,initial_inventory,molfactor,max_step_min)
    # Exposure stays between the inlet and fully reacted hazard bounds at
    # every cell; it is not transported with the liquid.
    bound=max(0.,hmin*newt-minimum(result.state.H),maximum(result.state.H)-hmax*newt)/max(1.,hmax*newt)
    bound<=1e-10 || throw(StepFailure("Integrated exposure violates the physical age bounds."))
    next_gen=led.generated+result.generated
    next_out=led.effluent+result.effluent
    next_back=led.backdiffusion+result.backdiffusion
    balance=abs(inventory_norm(g,result.state.x)-initial_inventory+next_out+next_back-next_gen)
    norm_balance=balance/max(1.,g.tau_min,newt)
    norm_balance<=1e-9 || throw(StepFailure("Cumulative species balance exceeds tolerance."))
    led.t=newt
    led.generated=next_gen
    led.effluent=next_out
    led.backdiffusion=next_back
    led.accepted+=1
    led.max_iters=max(led.max_iters,result.iterations)
    led.min_step=min(led.min_step,step)
    led.max_step=max(led.max_step,step)
    led.max_stage=max(led.max_stage,result.stage_error)
    led.max_film=max(led.max_film,result.max_film)
    led.max_slope=max(led.max_slope,result.max_slope)
    led.max_balance=max(led.max_balance,norm_balance)
    led.max_balance_mol=max(led.max_balance_mol,balance*molfactor)
    led.max_bound=max(led.max_bound,bound)
    # Keep the nominal step when the sample interval clipped it.
    led.dt=min(max_step_min,max(led.dt,1.5step))
    return nothing
end

@noinline function integrate_to!(led,g,l,state,goal,hmin,hmax,initial_inventory,molfactor,max_step_min)
    while led.t<goal
        remaining=goal-led.t
        if remaining<=8eps(max(goal,1.))
            led.t=goal; break
        end
        # Startup can contain rapid liquid displacement and catalyst
        # damage. Independently refine that period without unnecessarily
        # restricting the whole 180-minute run to the startup step size.
        stage_cap,phase_stop = led.t<10. ? (max_step_min/4,10.) :
            led.t<30. ? (max_step_min/2,30.) : (max_step_min,Inf)
        step=min(led.dt,remaining,stage_cap,phase_stop-led.t)
        try
            result=advance(g,l,state.x,state.H,led.t,step)
            newt=step==remaining ? goal : step==phase_stop-led.t ? phase_stop : led.t+step
            accept_step!(led,g,result,step,newt,hmin,hmax,initial_inventory,molfactor,stage_cap)
            state=result.state
        catch e
            if e isa StepFailure || e isa SingularException
                led.rejected+=1; led.dt=step/2
                led.dt>=1e-7 || error("Unable to resolve transient at t=$(led.t) min: $(sprint(showerror,e))")
            else
                rethrow()
            end
        end
        led.accepted+led.rejected<=1_000_000 || error("Transient integration exceeded step limit.")
    end
    return state
end

"""Transient axial liquid balance with separate stationary catalyst exposure.

Inputs are kg catalyst, L of interbead liquid, L/min flow, mol/L feed,
minutes, and the frozen batch Parameters. Return snapshots (cells × times)
of local liquid/surface progress, integrated hazard, activity and rate, as
well as product inventory/generation/outlet/backdiffusion in mol. The
exposure age is elapsed time at the fixed cell, never liquid residence.
"""
function simulate(p::CK.Parameters,T,feed,K,times;
                  catalyst_kg,liquid_L,flow_L_min,acid0,alcohol0,capacity,
                  eta_internal,peclet,cells=128,max_step_min=1.,
                  surface_hazard=false,recovered=false,equilibrium_prefill=false,
                  use_film=true,inlet_flux=false)
    ts=Float64.(collect(times)); n=Int(cells)
    !isempty(ts) && first(ts)==0 && all(isfinite,ts) && all(>(0.),diff(ts)) ||
        throw(ArgumentError("Sample times must start at zero and increase strictly."))
    all(isfinite,(T,K,catalyst_kg,liquid_L,flow_L_min,acid0,alcohol0,capacity,eta_internal,peclet,max_step_min)) &&
        T>-273.15 && feed in (1,2) && K>0 && catalyst_kg>=0 && liquid_L>0 && flow_L_min>0 &&
        acid0>0 && isapprox(alcohol0/acid0,feed;atol=1e-12,rtol=0) && capacity>0 &&
        0<eta_internal<=1 && max_step_min>0 || throw(ArgumentError("Invalid physical reactor inputs."))
    s=CK.static_terms(p,T,alcohol0/acid0,K;feed)
    all(isfinite,(s.k,s.equilibrium,s.hazard_at_ref_min_inv,s.deficit,s.tau_min,s.residual_active_fraction)) &&
        s.k>=0 && 0<s.equilibrium<1 && s.hazard_at_ref_min_inv>=0 && s.hazard_exponent>=0 &&
        s.protection_exponent>=0 && 0<=s.deficit<=1 && s.tau_min>0 && 0<=s.residual_active_fraction<=1 ||
        throw(ArgumentError("Invalid frozen rate/activity parameters."))
    tau=liquid_L/flow_L_min
    cat=1000*catalyst_kg/liquid_L
    g=AxialGrid(n,tau,peclet;inlet_flux)
    l=LocalLaw(s,eta_internal*.06*cat,acid0,capacity,use_film,surface_hazard,recovered)
    hmax=hazard(l,0.); hmin=hazard(l,s.equilibrium)
    all(isfinite,(hmax,hmin)) && 0<=hmin<=hmax || throw(ArgumentError("Invalid exposure bounds."))

    x=fill(equilibrium_prefill ? s.equilibrium : 0.,n); H=zeros(n)
    state=evaluate(g,l,x,H,0.,0.,x)
    columns=length(ts)
    bulk=zeros(n,columns); surface=zeros(n,columns); hazards=zeros(n,columns)
    activity_profile=zeros(n,columns); local_rate=zeros(n,columns); local_hazard=zeros(n,columns)
    gen=zeros(columns); out=zeros(columns); back=zeros(columns); inventory=zeros(columns)
    # tau*mean(X) is product moles in the void liquid / (acid0*Q); multiply
    # every ledger quantity by acid0*Q to export physical mol.
    initial_inventory=inventory_norm(g,x)
    molfactor=acid0*flow_L_min
    dt=min(max_step_min/8,.1tau,.2/max(hmax,1e-12))
    led=HistoryLedger(0.,dt,0.,0.,0.,0,0,0,Inf,0.,0.,0.,0.,0.,state.max_film,state.max_slope)
    function snapshot!(i,state)
        bulk[:,i]=state.x; surface[:,i]=state.surf; hazards[:,i]=state.H
        activity_profile[:,i]=state.a; local_rate[:,i]=acid0.*state.q; local_hazard[:,i]=state.h
        gen[i]=molfactor*led.generated; out[i]=molfactor*led.effluent; back[i]=molfactor*led.backdiffusion
        inventory[i]=molfactor*inventory_norm(g,state.x)
    end
    snapshot!(1,state)
    for i in 2:columns
        state=integrate_to!(led,g,l,state,ts[i],hmin,hmax,initial_inventory,molfactor,max_step_min)
        snapshot!(i,state)
    end
    (; times_min=ts,cell_center_fraction=(collect(1:n).-0.5)./n,bulk_progress=bulk,surface_progress=surface,
       cumulative_hazard=hazards,activity=activity_profile,local_rate_mol_L_min=local_rate,
       local_hazard_min_inv=local_hazard,outlet_conversion=vec(bulk[end,:]),
       generated_mol=gen,effluent_product_mol=out,backdiffused_product_mol=back,
       liquid_product_mol=inventory,initial_liquid_product_mol=initial_inventory*molfactor,
       catalyst_mass_kg=catalyst_kg,catalyst_mass_per_cell_kg=catalyst_kg/n,
       catalyst_g_L=cat,residence_min=tau,equilibrium_conversion=s.equilibrium,
       inlet_hazard_min_inv=hmax,equilibrium_hazard_min_inv=hmin,accepted_steps=led.accepted,
       rejected_steps=led.rejected,min_step_min=led.accepted==0 ? 0. : led.min_step,max_step_min=led.max_step,
       max_newton_iterations=led.max_iters,max_stage_residual=led.max_stage,max_normalized_balance_error=led.max_balance,
       max_absolute_balance_mol=led.max_balance_mol,max_normalized_exposure_bound_error=led.max_bound,
       max_film_error_mol_L_min=led.max_film,max_positive_film_slope_ratio=led.max_slope)
end

end
