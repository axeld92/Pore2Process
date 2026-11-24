using Test
include("../pentyl_propanoate/model.jl")
include("../pentyl_propanoate_loading/model.jl")
include("../pentyl_propanoate_composition/model.jl")
include("../pentyl_propanoate_transfer/published/model.jl")
include("../pentyl_propanoate_film_water/film/model.jl")
include("../pentyl_propanoate_equimolar/model.jl")
include("model.jl")

const S=SpatialCatalystHistory
const CK=CompositionKinetics

parameters(;k=.14,slope=0.,kb=.21,acid=0.,alcohol=0.,product=0.,h=0.,hs=0.,
            d=0.,tau=15.,m=3.5,protect=0.,floor=0.) =
    CK.Parameters(k,slope,kb,slope,acid,alcohol,product,0.,h,hs,d,tau,3.5,m,protect,floor)

function mass_balance(r,acid0,flow)
    @test r.max_normalized_balance_error<1e-9
    @test r.max_stage_residual<1e-10
    @test r.max_film_error_mol_L_min<1e-10
    @test maximum(abs.(r.liquid_product_mol .- r.initial_liquid_product_mol .+
                       r.effluent_product_mol .+ r.backdiffused_product_mol .- r.generated_mol))<1e-8*max(1.,acid0*flow*last(r.times_min))
    @test all(diff(r.cumulative_hazard;dims=2).>=-1e-12)
    @test all(0 .<= r.activity .<= 1)
    @test all(0 .<= r.bulk_progress .<= r.equilibrium_conversion)
end

@testset "Stationary constant exposure and blank reaction" begin
    p=parameters(k=0.,kb=0.,h=.03,m=0.,d=.45,tau=8.,floor=.05)
    times=[0.,1.,6.,10.,30.,180.]
    kw=(;catalyst_kg=.006,liquid_L=.018,flow_L_min=.003,acid0=5.3,alcohol0=5.3,
        eta_internal=.893,capacity=40.,peclet=115.,cells=32,max_step_min=2.)
    for recovered in (false,true)
        r=S.simulate(p,60.,1,12.,times;kw...,recovered)
        # Analytic and independent of age of the liquid. Identical hazard
        # at every fixed catalyst position for a constant imposed damage rate.
        expected_H=repeat(reshape(.03 .* times,1,:),32,1)
        expected_a=[(recovered ? 1. : 1-.45exp(-t/8.))*(.05+.95exp(-.03t)) for t in times]
        @test maximum(abs.(r.cumulative_hazard.-expected_H))<1e-11
        @test maximum(abs.(r.activity.-reshape(expected_a,1,:)))<1e-12
        @test all(iszero,r.bulk_progress)
        @test all(iszero,r.generated_mol)
        @test r.catalyst_mass_per_cell_kg*32≈.006 rtol=1e-15
        mass_balance(r,5.3,.003)
    end
end

@testset "Analytic reversible batch progress with exposure, closed-liquid limit" begin
    # With h constant the time integral of activity has an independent
    # elementary expression, including retained population and recovery.
    k=.19;h=.027;f=.09;delta=.35;tr=11.;R=1.;K=12.;C=40.;A=5.3
    p=parameters(k=k,kb=k,h=h,m=0.,d=delta,tau=tr,floor=f)
    effective(t)=f*t+(1-f)*(-expm1(-h*t))/h-delta*(f*tr*(-expm1(-t/tr))+
        (1-f)*(-expm1(-(h+inv(tr))*t))/(h+inv(tr)))
    time=[0.,.1,1.,5.,15.,60.]
    kw=(;catalyst_kg=C/1000,liquid_L=1.,flow_L_min=1e-12,acid0=A,alcohol0=A,
        capacity=100.,eta_internal=1.,peclet=115.,cells=8,use_film=false)
    coarse=S.simulate(p,60.,1,K,time;kw...,max_step_min=.2)
    fine=S.simulate(p,60.,1,K,time;kw...,max_step_min=.1)
    oracle=[CatalystLoading.progress_conversion(k,K,R,.06*C/A*effective(t)) for t in time]
    ec=maximum(abs.(coarse.outlet_conversion.-oracle)); ef=maximum(abs.(fine.outlet_conversion.-oracle))
    @test ef<ec/2
    @test ef<3e-6
    @test maximum(abs.(fine.cumulative_hazard .- reshape(h.*time,1,:)))<1e-11
    @test maximum(abs.(fine.bulk_progress.-reshape(oracle,1,:)))<3e-6
    mass_balance(fine,A,1e-12)
end

@testset "Independent plug-flow startup and unchanged steady finite-Pe calculation" begin
    p=parameters(k=.15,kb=.18);R=2.;K=6.;A=3.3;C=333.;tau=6.
    kw=(;catalyst_kg=C*.018/1000,liquid_L=.018,flow_L_min=.018/tau,acid0=A,alcohol0=6.6,
        capacity=35.,eta_internal=1.)
    # Advection without dispersion: fluid age at the outlet is min(t,tau)
    # when fresh mixture fills the entire bed at t=0.
    t=[0.,1.,2.,3.,10.,18.]
    r=S.simulate(p,60.,2,K,t;kw...,peclet=1e6,cells=512,max_step_min=.04,use_film=false)
    want=[CatalystLoading.progress_conversion(.15,K,R,.06*C/A*min(u,tau)) for u in t]
    @test r.outlet_conversion≈want atol=.001
    mass_balance(r,A,.018/tau)
    # With dispersion, retain the identical finite-volume inlet convention
    # and film kinetics; the no-aging transient should reach the old steady answer.
    for flux in (false,true)
        t=[0.,10.,30.,60.,120.,180.]
        r=S.simulate(p,60.,2,K,t;kw...,peclet=115.,cells=128,max_step_min=.5,inlet_flux=flux)
        m=DualFeedBed.FrozenKinetics(p,K,1.,2)
        f=DualFeedBed.WithFilm(m,kw.capacity,1.)
        prev=LorenzoPackedBed.packed_bed(f,60.,tau,C,1.,115.;acid0=A,alcohol0=6.6,cells=128,inlet_flux=flux)
        @test last(r.outlet_conversion)≈prev.conversion atol=2e-9
        @test all(==(1.),r.activity)
        @test r.backdiffused_product_mol[end]>=0.
        @test (flux ? r.backdiffused_product_mol[end]==0. : r.backdiffused_product_mol[end]>0.)
        mass_balance(r,A,.018/tau)
    end
end

@testset "Implicit film source derivative and surface exposure" begin
    p=parameters(h=.015,m=3.5,protect=2.,acid=.6,alcohol=.3,product=-.2,d=.2,floor=.05)
    st=CK.static_terms(p,70.,1.,12.;feed=1)
    for surface in (false,true), x in (.1,.4,.65)
        l=S.LocalLaw(st,18.,5.3,40.,true,surface,false)
        local_state=S.local_state(l,.3,.4,6.,x)
        lo=S.local_state(l,.3,.4,6.,x-1e-5); hi=S.local_state(l,.3,.4,6.,x+1e-5)
        fd=(hi.source-lo.source)/(2e-5)
        @test local_state.slope≈fd rtol=5e-5 atol=1e-7
        @test local_state.film_error<1e-10
        @test local_state.H≈.4+.3*CK.acid_hazard(st,surface ? local_state.surface : x) rtol=1e-13
        @test local_state.surface>=x
    end
    screen=S.film_slope_screen(p,70.,1,12.;catalyst_g_L=333.,capacity=40.,eta_internal=.99,surface_hazard=true)
    @test screen.max_sampled_slope_ratio<.99
end

@testset "Initial product ledger and invalid cases" begin
    p=parameters(k=0.,kb=0.)
    kw=(;catalyst_kg=0.,liquid_L=.018,flow_L_min=.003,acid0=3.3,alcohol0=6.6,
        capacity=50.,eta_internal=1.,peclet=115.,cells=64,max_step_min=.5,equilibrium_prefill=true)
    r=S.simulate(p,60.,2,6.,[0.,1.,6.,30.];kw...)
    @test r.initial_liquid_product_mol≈.018*3.3*r.equilibrium_conversion rtol=1e-14
    @test all(iszero,r.generated_mol)
    @test last(r.liquid_product_mol)<1e-5*r.initial_liquid_product_mol
    @test r.effluent_product_mol[end]+r.backdiffused_product_mol[end]≈r.initial_liquid_product_mol atol=1e-6
    mass_balance(r,3.3,.003)
    @test_throws ArgumentError S.simulate(p,60.,2,6.,[1.,2.];kw...)
    @test_throws ArgumentError S.simulate(p,60.,1,6.,[0.,1.];kw...)
    @test_throws ArgumentError S.simulate(p,60.,2,-6.,[0.,1.];kw...)
end
