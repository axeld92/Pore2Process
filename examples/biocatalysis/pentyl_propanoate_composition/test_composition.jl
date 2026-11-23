using Test
include("../pentyl_propanoate/model.jl")
include("../pentyl_propanoate_loading/model.jl")
include("model.jl")
include("../pentyl_propanoate_transfer/published/model.jl")
include("../pentyl_propanoate_film_water/film/model.jl")
include("bed_transfer.jl")
using .CompositionKinetics
using .CompositionBedTransfer

parameters(;k=.14,slope=.3,acid=.0,alcohol=.0,product=.0,heat=.0,h=.0,hs=.0,delta=.0,tau=8.,n=3.5,m=3.5,protect=.0,floor=.0) =
    Parameters(k,slope,k,slope,acid,alcohol,product,heat,h,hs,delta,tau,n,m,protect,floor)

@testset "Stoichiometric reversible state and loading invariant" begin
    for R in (1.,2.,.98)
        K=5.7
        p=parameters(acid=.8,alcohol=.2,product=-.4)
        state=static_terms(p,60.,R,K)
        eq=equilibrium_conversion(R,K)
        @test abs(net_mass_rate(state,eq)) < 1e-14
        @test net_mass_rate(state,eq*.9)>0
        @test net_mass_rate(state,eq+min(.03,(min(1.,R)-eq)/2))<0
        low=batch_progress(p,60.,R,5.,K,[20.,60.,100.,150.])
        high=batch_progress(p,60.,R,20.,K,[5.,15.,25.,37.5])
        @test low.conversion ≈ high.conversion atol=2e-6
        @test all(==(1.),high.activity)
    end
    s=static_terms(parameters(acid=20.,alcohol=3.),60.,2.,5.)
    @test net_mass_rate(s,0.) ≈ .14*2/9 rtol=1e-14
    @test batch_progress(parameters(),60.,2.,0.,5.,[60.]).conversion==[0.]
    @test_throws ArgumentError batch_progress(parameters(),60.,2.,10.,5.,[20.,10.])
end

@testset "Independent constant-hazard and recovery oracle" begin
    # With a constant hazard (m=0) and no composition inhibition,
    # integrate the product of exponential decay and recovery exactly.
    k=.13; h=.012; delta=.45; tau=15.; w=10.; R=1.; K=10.
    p=parameters(k=k,slope=0.,h=h,m=0.,delta=delta,tau=tau)
    t=[0.,.1,5.,15.,40.,120.,300.]
    progress=batch_progress(p,60.,R,w,K,t)
    effective(u)=-expm1(-h*u)/h+delta*expm1(-(h+inv(tau))*u)/(h+inv(tau))
    exact=[CatalystLoading.progress_conversion(k,K,R,60*CompositionKinetics.MA*w/100*effective(u)) for u in t]
    @test progress.conversion ≈ exact atol=1e-6
    exact_a=[(1-delta*exp(-u/tau))*exp(-h*u) for u in t]
    @test progress.activity ≈ exact_a atol=1e-10
    @test progress.hazards ≈ h.*t atol=1e-10
    # Same catalyst mass multiplied by elapsed time is no longer the
    # whole history when stationary/moving catalyst clock is involved.
    a=batch_progress(p,60.,R,5.,K,[80.])
    b=batch_progress(p,60.,R,20.,K,[20.])
    @test abs(only(a.conversion)-only(b.conversion))>1e-3
    # Pure recovery with zero damage also has analytic exposure.
    q=parameters(k=k,slope=0.,delta=delta,tau=tau)
    test=batch_progress(q,60.,R,w,K,t)
    @test test.conversion ≈ [CatalystLoading.progress_conversion(k,K,R,
        60*CompositionKinetics.MA*w/100*(u+delta*tau*expm1(-u/tau))) for u in t] atol=1e-12
    terms=static_terms(parameters(h=.01),60.,2.,5.)
    envelope=activity_bounds(terms,180.)
    @test envelope.lower≈exp(-.01*180.)
    @test 0<envelope.lower<envelope.upper<1
    protected=parameters(h=.01,protect=4.)
    st2=static_terms(protected,60.,2.,5.)
    st1=static_terms(protected,60.,1.,12.;feed=1)
    unprotected=static_terms(parameters(h=.01),60.,1.,12.;feed=1)
    @test acid_hazard(st2,0.)≈.01
    @test acid_hazard(st1,0.)/acid_hazard(unprotected,0.)≈16.
    @test acid_hazard(st2,.5)<acid_hazard(static_terms(parameters(h=.01),60.,2.,5.),.5)
    # A stable residual population is tested separately with an exact
    # integral of a constant-hazard activity: f*t+(1-f)*(1-exp(-h*t))/h.
    stable=.08
    u=parameters(k=k,slope=0.,h=h,m=0.,floor=stable)
    y=batch_progress(u,60.,R,w,K,t)
    E=60*CompositionKinetics.MA*w/100
    @test y.conversion≈[CatalystLoading.progress_conversion(k,K,R,E*(stable*z-(1-stable)*expm1(-h*z)/h)) for z in t] atol=1e-6
    @test last(y.activity)>stable
    @test activity_bounds(static_terms(u,60.,2.,5.),1e8).lower≈stable
end

@testset "Equilibrium regression and bounded shared parameter recovery" begin
    temps=[40.,50.,60.]
    exact=[exp(1.6+.4*PentylPropanoate.temperature_coordinate(t)) for t in temps]
    eq=fit_equilibrium(temps,exact)
    @test eq.log_K_reference≈1.6 atol=1e-12
    @test eq.inverse_temperature_slope≈.4 atol=1e-12
    @test_throws ArgumentError fit_equilibrium([40.],[4.])
    # Exact bounded optimum against a synthetic least-squares problem.
    f=CompositionKinetics.bounded_fit(q->[q[1]-.4,q[2]+.2],[.8,.7])
    @test f.converged
    @test f.q≈[.4,0.] atol=1e-5
    # Numerical recovery across feeds with no held-out loading in fitting.
    p=parameters(k=.16,slope=.4,acid=.55)
    ts=[4.,12.,25.,55.,100.,180.,300.]
    ids=[1,2,3,4];feeds=[1,1,2,2];T=[40.,60.,60.,80.]
    ratios=Float64.(feeds);loading=[5.,10.,5.,10.];K=[11.,17.,6.,7.]
    X=[batch_progress(p,T[i],ratios[i],loading[i],K[i],ts;feed=feeds[i],
        rtol=2e-10,atol=2e-12,max_step_min=5.).conversion for i in 1:4]
    curves=make_curves(ids,feeds,T,ratios,loading,K,fill(length(ts),4),repeat(ts,4),vcat(X...))
    fit=fit_model("acid",curves)
    @test fit.converged
    @test fit.parameters.k2≈.16 rtol=1e-3
    @test fit.parameters.k2_slope≈.4 atol=1e-3
    @test fit.parameters.acid_b≈.55 atol=1e-3
    got=batch_progress(fit.parameters,60.,1.,20.,17.,ts;feed=1)
    wanted=batch_progress(p,60.,1.,20.,17.,ts;feed=1,rtol=2e-10,atol=2e-12)
    @test maximum(abs.(got.conversion-wanted.conversion))<1e-4
end

@testset "Existing axial and film balances with new mass-specific kinetics" begin
    T,K,C,t,p=60.,6.,333.,6.,parameters()
    acid,alcohol,Pe=3.3,6.6,1e4
    model=CompositionBedTransfer.BedKinetics(p,K,1.)
    q=LorenzoPackedBed.packed_bed(model,T,t,C,1.,Pe;acid0=acid,alcohol0=alcohol,cells=1024)
    exact=CatalystLoading.progress_conversion(.14,K,2.,60*C/1000*t/acid)
    @test q.conversion≈exact atol=1e-3
    @test q.balance_error<1e-9
    # Do not alter the mass to change activity.
    a=parameters(acid=.2)
    cases=[frozen_prediction(a,60.,K,f,t,C,3.0,.869,115.;cells=512) for f in (0.,.5,1.)]
    @test abs(cases[1].conversion)<1e-12
    @test 0<cases[2].conversion<cases[3].conversion<equilibrium_conversion(2.,K)
    @test all(r->r.balance_error<1e-9,cases)
    @test all(r->r.maximum_sampled_film_residual_mol_L_min<1e-10,cases)
end
