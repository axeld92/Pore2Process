using Test
include("../pentyl_propanoate/model.jl")
include("../pentyl_propanoate_loading/model.jl")
include("../pentyl_propanoate_composition/model.jl")
include("../pentyl_propanoate_transfer/published/model.jl")
include("../pentyl_propanoate_film_water/film/model.jl")
include("../pentyl_propanoate_composition/bed_transfer.jl")
include("model.jl")
using .DualFeedBed

const CK=CompositionKinetics
const Bed=LorenzoPackedBed

@testset "Separate nominal feeds and consistent mass-specific rate" begin
    # Deliberately unequal reference coefficients so using the 2:1 rate
    # silently for the equimolar feed cannot pass.
    p=CK.Parameters(.14,0.,.35,0.,0.,0.,0.,0.,0.,0.,0.,8.,3.5,3.5,0.,0.)
    C=333.05
    for (R,A,K,k) in ((1,5.3,12.,.35),(2,3.3,6.,.14))
        kin=FrozenKinetics(p,K,1.,R)
        terms=Bed.rate_terms(kin,60.,A,R*A,false)
        @test terms.k1==k
        @test terms.k1/terms.k2≈K rtol=1e-14
        exact=.06*C*k*R/(1+R)^2
        @test Bed.rate(kin,terms,A,R*A,C,0.)≈exact rtol=1e-14
        @test abs(Bed.rate(kin,terms,A,R*A,C,CK.equilibrium_conversion(R,K)))<1e-12
        inactive=FrozenKinetics(p,K,0.,R)
        @test Bed.rate(inactive,terms,A,R*A,C,0.)==0.
        # The existing axial solver at large Pe tends to the independent
        # analytic plug-flow reversible solution: no redesign of solver.
        t=6.
        result=prediction(kin,60.,t,C,10.,1.,1e5;acid0=A,alcohol0=R*A,cells=1024,use_film=false)
        oracle=CatalystLoading.progress_conversion(k,K,R,.06*C*t/A)
        @test result.conversion≈oracle atol=1e-3
        @test result.balance_error<1e-9
        @test result.equilibrium_conversion≈CK.equilibrium_conversion(R,K) rtol=1e-14
        @test result.integrated_source_fraction≈result.conversion+result.inlet_backdiffusion_fraction atol=1e-9
        zero=prediction(inactive,60.,t,C,10.,1.,115.;acid0=A,alcohol0=R*A,cells=256)
        @test zero.conversion==0.
    end
    @test_throws ArgumentError Bed.rate_terms(FrozenKinetics(p,12.,1.,1),60.,3.3,6.6,false)
    @test_throws ArgumentError FrozenKinetics(p,0.,1.,1)
    @test_throws ArgumentError FrozenKinetics(p,12.,1.01,1)
    @test_throws ArgumentError FrozenKinetics(p,12.,1.,3)
end

@testset "Actual feed film root and unchanged prior 2:1 path" begin
    p=CK.Parameters(.14,.3,.35,.2,.8,.2,-.4,0.,0.,0.,0.,8.,3.5,3.5,0.,0.)
    C=333.05; cap=10.; eta=.893; pe=115.; t=9.
    m=FrozenKinetics(p,12.,1.,1)
    f=WithFilm(m,cap,eta)
    A=5.3; B=A
    terms=Bed.rate_terms(f,60.,A,B,false)
    for x in (0.,.2,.4,.65)
        s=BeadFilm.surface_state(f,terms,A,B,C,x)
        @test s.surface>=x
        @test abs(s.reaction-cap*(s.surface-x))<1e-10
        @test eta*Bed.rate(f,terms,A,B,C,x)≈s.reaction rtol=1e-14
    end
    lo=prediction(m,60.,t,C,cap/2,eta,pe;acid0=A,alcohol0=B,cells=512)
    hi=prediction(m,60.,t,C,cap*2,eta,pe;acid0=A,alcohol0=B,cells=512)
    free=prediction(m,60.,t,C,cap,eta,pe;acid0=A,alcohol0=B,cells=512,use_film=false)
    @test lo.conversion < hi.conversion < free.conversion
    @test hi.max_film_balance_error_mol_L_min<1e-10
    @test hi.max_positive_rate_slope_over_film_capacity<.99
    @test_throws ArgumentError WithFilm(m,0.,eta)
    @test_throws ArgumentError WithFilm(m,cap,1.1)
    # Use the same 2:1 feed, model coefficients, effectiveness and Pe as
    # the old frozen wrapper; unchanged physical calculation.
    prev=CompositionBedTransfer.frozen_prediction(p,60.,6.,.9,t,C,cap,.869,pe;cells=512)
    current=prediction(FrozenKinetics(p,6.,.9,2),60.,t,C,cap,.869,pe;
        acid0=3.3,alcohol0=6.6,cells=512)
    @test prev.conversion≈current.conversion atol=1e-12
    @test current.balance_error<1e-9
end

@testset "Correct dimensionless film transfer and stationary clock" begin
    u=2e-4; d=593e-6; D=1e-9; eps=.64
    a=BeadFilm.film_coefficients(u,d,800.,.001,D,eps,"wilson_geankoplis")
    b=BeadFilm.film_coefficients(u,d,860.,.002,D,eps,"wilson_geankoplis")
    @test a.reynolds*a.schmidt≈u*d/D rtol=1e-14
    @test a.coefficient_m_s≈b.coefficient_m_s rtol=1e-14
    c=BeadFilm.film_coefficients(u,d,860.,.002,2D,eps,"wilson_geankoplis")
    @test c.coefficient_m_s/a.coefficient_m_s≈2^(2/3) rtol=1e-14
    p=CK.Parameters(.14,0.,.14,0.,0.,0.,0.,0.,.001,.0,0.,8.,3.5,3.5,4.,.05)
    s1=CK.static_terms(p,60.,1.,12.;feed=1)
    s2=CK.static_terms(p,60.,2.,6.;feed=2)
    a1=CK.activity_bounds(s1,30.); a2=CK.activity_bounds(s2,30.)
    @test a1.lower<a2.lower
    @test a1.lower>p.residual_active_fraction
    @test a1.lower<CK.activity_bounds(s1,0.).lower
    @test CK.activity_bounds(s1,180.).lower<a1.lower
end
