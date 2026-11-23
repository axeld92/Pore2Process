module DualFeedBed

import ..CompositionKinetics
import ..LorenzoPackedBed: rate, rate_terms, packed_bed
import ..BeadFilm

export FrozenKinetics, WithFilm, prediction

"""Frozen batch law, with explicit nominal feed and separate stationary activity."""
struct FrozenKinetics
    parameters::CompositionKinetics.Parameters
    Kx::Float64
    activity::Float64
    feed::Int
    function FrozenKinetics(p,Kx,activity,feed)
        isfinite(Kx) && Kx>0 && isfinite(activity) && 0<=activity<=1 && feed in (1,2) ||
            throw(ArgumentError("Require positive Kx, activity in [0,1] and feed ratio 1 or 2."))
        new(p,Kx,activity,feed)
    end
end

function rate_terms(p::FrozenKinetics,T,acid0,alcohol0,temperature_dependent_ester_term)
    isapprox(alcohol0/acid0,p.feed;atol=1e-12,rtol=0) ||
        throw(ArgumentError("Feed concentrations do not match the nominal calibration ratio."))
    st=CompositionKinetics.static_terms(p.parameters,T,alcohol0/acid0,p.Kx;feed=p.feed)
    # The existing axial solver obtains its equilibrium bracket from k1/k2.
    (;k1=st.k,k2=st.k/st.Kx,state=st)
end

function rate(p::FrozenKinetics,terms,acid0,alcohol0,catalyst_g_L,x)
    # mol/(kg_cat s) * g_cat/L_liquid * 60/1000 -> mol/(L_liquid min).
    .06*catalyst_g_L*p.activity*CompositionKinetics.net_mass_rate(terms.state,x)
end

"""Use the existing generic BeadFilm local surface state and root solver."""
struct WithFilm{K}
    kinetics::K
    capacity_mol_L_min::Float64
    eta_internal::Float64
    function WithFilm(kinetics::K,capacity,eta) where K
        isfinite(capacity) && capacity>0 && isfinite(eta) && 0<eta<=1 ||
            throw(ArgumentError("Require positive film capacity and eta in (0,1]."))
        new{K}(kinetics,capacity,eta)
    end
end

rate_terms(f::WithFilm,args...)=rate_terms(f.kinetics,args...)

function rate(f::WithFilm,terms,acid0,alcohol0,catalyst_g_L,x)
    # packed_bed multiplies the returned rate by eta, once. The same eta
    # appears in the local film root; no separate printed eta_external.
    BeadFilm.surface_state(f,terms,acid0,alcohol0,catalyst_g_L,x).reaction/f.eta_internal
end

equilibrium(terms,R)=CompositionKinetics.equilibrium_conversion(R,terms.k1/terms.k2)

function film_audit(f::WithFilm,terms,acid0,alcohol0,catalyst_g_L)
    eq=equilibrium(terms,alcohol0/acid0)
    reaction(x)=f.eta_internal*rate(f.kinetics,terms,acid0,alcohol0,catalyst_g_L,x)
    # Check the same sufficient sampled monotonicity condition as the
    # previous 2:1 example, refined over the entire physical interval.
    slopes(n)=maximum((reaction(x+1e-6)-reaction(x-1e-6))/(2e-6)
                      for x in range(0.,eq;length=n))
    coarse,fine=slopes(65),slopes(257)
    bulk=range(0.,eq;length=33)
    states=[BeadFilm.surface_state(f,terms,acid0,alcohol0,catalyst_g_L,x) for x in bulk]
    residuals=[abs(s.reaction-f.capacity_mol_L_min*(s.surface-x)) for (s,x) in zip(states,bulk)]
    (;surface_progress_at_inlet=first(states).surface,
      max_film_balance_error_mol_L_min=maximum(residuals),
      max_positive_rate_slope_over_film_capacity=max(0.,fine/f.capacity_mol_L_min),
      film_slope_sampling_refinement_change=abs(fine-coarse)/f.capacity_mol_L_min)
end

function prediction(kinetics,T,residence,catalyst_g_L,capacity,eta,peclet;
                    acid0,alcohol0,cells=1024,use_film=true)
    model=use_film ? WithFilm(kinetics,capacity,eta) : kinetics
    terms=rate_terms(kinetics,T,acid0,alcohol0,false)
    details=use_film ? film_audit(model,terms,acid0,alcohol0,catalyst_g_L) :
        (;surface_progress_at_inlet=0.,max_film_balance_error_mol_L_min=0.,
          max_positive_rate_slope_over_film_capacity=0.,film_slope_sampling_refinement_change=0.)
    result=packed_bed(model,T,residence,catalyst_g_L,eta,peclet;acid0,alcohol0,cells)
    (;result...,equilibrium_conversion=equilibrium(terms,alcohol0/acid0),details...)
end

end
