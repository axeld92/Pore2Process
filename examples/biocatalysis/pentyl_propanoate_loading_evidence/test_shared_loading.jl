using Test
include("../pentyl_propanoate/model.jl")
include("../pentyl_propanoate_loading/model.jl")
include("shared_loading.jl")
using .SharedLoadingEvidence

@testset "Shared shape, separate chemical equilibria" begin
    ref = 25.
    C = vcat(repeat([6.25,12.5,25.,50.],inner=6),repeat([10.,20.,40.,80.],inner=6))
    t = repeat([10.,25.,50.,100.,200.,300.],8)
    g = repeat([2,1],inner=24)
    # Exposure (kg catalyst*s/mol acid) differs at the two feed ratios
    # even at the same concentration. Generate the data with independent
    # integration of the mole-fraction balance rather than predict_series.
    E = t .* C .* ifelse.(g .== 2, .017, .011)
    R = Float64.(g)
    k = Dict(1=>0.167,2=>0.13)
    K = Dict(1=>10.3,2=>4.7)
    shape = .19
    function independent(C,E,R,g)
        dt = E/2048
        rhs(x) = k[g]/(1+shape*(C/ref-1))*((1-x)*(R-x)-x^2/K[g])/(1+R)^2
        x = 0.
        for _ in 1:2048
            a = rhs(x); b = rhs(x+dt*a/2); c = rhs(x+dt*b/2); d = rhs(x+dt*c)
            x += dt*(a+2b+2c+d)/6
        end
        x
    end
    X = [independent(c,e,r,b) for (c,e,r,b) in zip(C,E,R,g)]
    train = C .!= 80
    result = fit_shared("crowding",C[train],E[train],R[train],X[train],ones(count(train)),g[train];
        reference_concentration=ref)
    @test result.shape ≈ shape atol=2e-5
    for (group,f) in zip(result.series,result.fits)
        @test f.k_ref ≈ k[group] rtol=2e-5
        @test f.Kx ≈ K[group] rtol=2e-5
    end
    f = result.fits[findfirst(==(1),result.series)]
    y = predict_series("crowding",result.shape,f.k_ref,f.Kx,C[.!train],E[.!train],R[.!train];
        reference_concentration=ref)
    @test maximum(abs.(y-X[.!train])) < 1e-6
    @test_throws ArgumentError fit_shared("crowding",C,E,R,X,ones(length(C)),g;
        reference_concentration=0)
end

@testset "Unidentified equilibrium reported as a numerical limit" begin
    exposures=[1.,4.,12.,25.,50.,100.,300.]
    k=.14
    # Independent solution to dX/dE=k*(1-X)^2/4 (no reverse rate).
    x=k.*exposures ./ (4 .+ k.*exposures)
    fit=fit_shared("proportional",fill(25.,7),exposures,ones(7),x,ones(7),ones(Int,7);
        reference_concentration=25.)
    p=only(fit.fits)
    @test p.solver == "upper_equilibrium_limit"
    @test p.Kx == 1e9 || p.Kx ≈ 1e9
    @test p.k_ref ≈ k rtol=1e-5
    @test p.limit_probe_change < 1e-5
    predictions=predict_series("proportional",1.,p.k_ref,p.Kx,fill(25.,7),exposures,ones(7);
        reference_concentration=25.)
    @test maximum(abs.(predictions-x)) < 1e-6
end

@testset "Matched-progress laws across independent chemistries" begin
    u = repeat([.1,.5,1.,1.5,2.],3)
    blocks = repeat([1,2,3],inner=5)
    anchors = [3.,16.,38.]
    for (law,shape) in (("proportional",1.),("power",.71),("crowding",.23))
        # Independent closed-form crossing times with three arbitrary
        # matched-progress times. Hold out the largest loading everywhere.
        times = [anchors[b]*(law=="crowding" ? (1+shape*(c-1))/c : c^(-shape))
                 for (b,c) in zip(blocks,u)]
        train = u .!= 2.
        fit = fit_crossing_law(law,u[train],times[train],ones(count(train)),blocks[train])
        @test fit.shape ≈ shape atol=2e-6
        prediction = [exp(fit.intercepts[b]-log_total_rate(law,fit.shape,c))
                      for (b,c) in zip(blocks[.!train],u[.!train])]
        @test prediction ≈ times[.!train] rtol=1e-6
    end
    @test_throws ArgumentError fit_crossing_law("power",[1.,1.],[2.,2.],[1.,1.],[1,1])
end
