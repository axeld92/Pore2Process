module CompositionBedTransfer

import ..CompositionKinetics
import ..LorenzoPackedBed: rate, rate_terms, packed_bed
import ..BeadFilm

export BedKinetics, CompositionFilm, inlet_diagnostics, frozen_prediction

struct BedKinetics
    parameters::CompositionKinetics.Parameters
    Kx::Float64
    activity::Float64
    function BedKinetics(p,Kx,activity)
        isfinite(Kx) && Kx>0 && isfinite(activity) && 0<=activity<=1 ||
            throw(ArgumentError("Use positive equilibrium and fractional stationary activity."))
        new(p,Kx,activity)
    end
end

function rate_terms(p::BedKinetics,T,acid0,alcohol0,temperature_dependent_ester_term)
    st=CompositionKinetics.static_terms(p.parameters,T,alcohol0/acid0,p.Kx;feed=2)
    # The existing finite-volume reactor solver locates equilibrium by these
    # two constants. It is the same reversible root, independent of activity.
    (; k1=st.k,k2=st.k/st.Kx,state=st)
end

function rate(p::BedKinetics,terms,acid0,alcohol0,catalyst_g_L,x)
    # Mass-specific rate [mol/kg_cat/s] -> mol/(L_liquid min).
    # The weighed inventory is supplied once by the existing axial balance.
    .06*catalyst_g_L*p.activity*CompositionKinetics.net_mass_rate(terms.state,x)
end

struct CompositionFilm
    kinetics::BedKinetics
    capacity_mol_L_min::Float64
    eta_internal::Float64
    function CompositionFilm(kinetics,capacity,eta)
        isfinite(capacity) && capacity>0 && isfinite(eta) && 0<eta<=1 ||
            throw(ArgumentError("Require positive film capacity and internal effectiveness."))
        new(kinetics,capacity,eta)
    end
end

rate_terms(f::CompositionFilm,args...)=rate_terms(f.kinetics,args...)

function surface_state(f::CompositionFilm,terms,acid0,alcohol0,catalyst_g_L,bulk)
    reaction(x)=f.eta_internal*rate(f.kinetics,terms,acid0,alcohol0,catalyst_g_L,x)
    surface=BeadFilm.surface_progress(reaction,bulk,terms.state.equilibrium,f.capacity_mol_L_min)
    (; surface,reaction=reaction(surface))
end

function rate(f::CompositionFilm,terms,acid0,alcohol0,catalyst_g_L,x)
    # Reuse existing film root; existing axial solver multiplies eta once.
    surface_state(f,terms,acid0,alcohol0,catalyst_g_L,x).reaction/f.eta_internal
end

"""Film root audit, including possible nonmonotonic reaction response."""
function inlet_diagnostics(f::CompositionFilm,T,catalyst_g_L;acid0=3.3,alcohol0=6.6)
    terms=rate_terms(f,T,acid0,alcohol0,false)
    eq=terms.state.equilibrium
    r(x)=f.eta_internal*rate(f.kinetics,terms,acid0,alcohol0,catalyst_g_L,x)
    # A sufficient sampled check for locally single-valued film roots:
    # d(k_f*X-r)/dX>0. Check along the whole physical interval,
    # not merely at the inlet; refine the stationary sampling.
    slopes(n) = maximum((r(x+1e-6)-r(x-1e-6))/(2e-6)
                        for x in range(0.,eq;length=n))
    c0,c1=slopes(65),slopes(257)
    balances=[surface_state(f,terms,acid0,alcohol0,catalyst_g_L,x)
              for x in range(0.,eq;length=33)]
    errors=[abs(b.reaction-f.capacity_mol_L_min*(b.surface-x))
            for (b,x) in zip(balances,range(0.,eq;length=33))]
    (; surface_progress_at_inlet=first(balances).surface,
       maximum_sampled_film_residual_mol_L_min=maximum(errors),
       maximum_positive_rate_slope_over_film_capacity=max(0.,c1/f.capacity_mol_L_min),
       film_slope_refinement_ratio_change=abs(c1-c0)/f.capacity_mol_L_min)
end

function frozen_prediction(p,T,Kx,activity,residence,catalyst_g_L,capacity,eta,peclet;
                           acid0=3.3,alcohol0=6.6,cells=1024)
    kinetics=BedKinetics(p,Kx,activity)
    film=CompositionFilm(kinetics,capacity,eta)
    result=packed_bed(film,T,residence,catalyst_g_L,eta,peclet;acid0,alcohol0,cells)
    diagnostics=inlet_diagnostics(film,T,catalyst_g_L;acid0,alcohol0)
    (; result.conversion,result.balance_error,diagnostics...)
end

end
