using Test
include("../pentyl_propanoate/model.jl")
include("model.jl")
include("../pentyl_propanoate_transfer/published/model.jl")
include("../pentyl_propanoate_film_water/film/model.jl")
include("film_transfer.jl")
using .CatalystLoading: rate_scale, progress_conversion, predict_batch, fit_loading

@testset "Mass-specific scaling and reversible progress" begin
    @test rate_scale("proportional",1,4) == 1
    @test rate_scale("power",1,4) == 1
    @test rate_scale("crowding",0,4) == 1
    @test rate_scale("power",0.8,4) ≈ 4^(-0.2)
    @test 4rate_scale("crowding",1,4) == 1
    @test rate_scale("crowding",0.2,2) ≈ 1/1.2
    @test_throws ArgumentError rate_scale("crowding",-0.2,2)
    @test_throws ArgumentError rate_scale("proportional",1,0)
    k,K = 0.13,4.53
    prior = PentylPropanoate.Kinetics(log(k),0,log(K),0)
    for exposure in (0.,1.,10.,40.,300.)
        @test progress_conversion(k,K,2.,exposure) ≈
            PentylPropanoate.conversion(prior,40.,exposure) rtol=1e-12
    end
    @test CatalystLoading.batch_conversion("proportional",1,k,K,2.5,80,2.) ≈
        CatalystLoading.batch_conversion("proportional",1,k,K,20.,10,2.) rtol=1e-14

    # Integrate the original reversible mole-fraction balance in exposure
    # independently of the analytical root/exponential expression.
    R, exposure = 1.907, 43.0
    rhs(x) = k*((1-x)*(R-x)-x^2/K)/(1+R)^2
    dt = exposure/4096
    x = 0.
    for _ in 1:4096
        a = rhs(x)
        b = rhs(x+dt*a/2)
        c = rhs(x+dt*b/2)
        d = rhs(x+dt*c)
        x += dt*(a+2b+2c+d)/6
    end
    @test progress_conversion(k,K,R,exposure) ≈ x atol=1e-12
    @test progress_conversion(k,K,R,1e6) ≈ 2R/(1+R+sqrt((R-1)^2+4R/K))
end

@testset "Recover independent synthetic loading curves" begin
    t = [5.,12.,24.,48.,96.,180.,300.]
    loading = repeat([2.5,5.,10.,20.],inner=length(t))
    times = repeat(t,4)
    ratios = fill(2.,length(loading))
    mask = loading .!= 5.
    holdout = .!mask
    k,K = 0.142,4.7
    for (law,shape) in (("proportional",1.0),("power",0.82),("crowding",0.21))
        # Generate a noiseless independent experimental pattern, then hide
        # an entire catalyst loading during the kinetic parameter recovery.
        observed = predict_batch(law,shape,k,K,loading,times,ratios)
        train = fit_loading(law,loading[mask],times[mask],ratios[mask],
                            observed[mask],ones(count(mask)))
        @test train.k_ref ≈ k rtol=2e-5
        @test train.Kx ≈ K rtol=2e-5
        @test train.shape ≈ shape atol=2e-5
        predicted = predict_batch(law,train.shape,train.k_ref,train.Kx,
                                  loading[holdout],times[holdout],ratios[holdout])
        @test maximum(abs.(predicted-observed[holdout])) < 1e-6
        if law == "proportional"
            for alternative in ("power","crowding")
                null = fit_loading(alternative,loading[mask],times[mask],ratios[mask],
                                   observed[mask],ones(count(mask)))
                @test null.shape == (alternative=="power" ? 1. : 0.)
            end
        end
    end
end

@testset "Rate multiplier acts inside existing film balance at fixed mass" begin
    p = LorenzoPackedBed.PingPong([3520.6,1746.6,260.76,1913.3,42250,1.2764e-6,103740,776560,
        848920,4.4007e-6,16898,12813,6.6093e-6,12.292,3.3629,-28091,313.15])
    T, residence, concentration, eta, pe = 60.,6.,333.,0.869,115.
    film = BeadFilm.FilmKinetics(p,2.,eta)
    solve(model) = LorenzoPackedBed.packed_bed(model,T,residence,concentration,eta,pe;cells=256)
    baseline = solve(film)
    results = [solve(LoadingFilmTransfer.LoadingAdjustedFilm(film,f)) for f in (0.7,1.,1.3)]
    @test results[2].conversion ≈ baseline.conversion atol=1e-12
    @test 0 < results[1].conversion < baseline.conversion < results[3].conversion < 1
    @test all(r -> r.balance_error < 1e-9,results)
    local_result = LoadingFilmTransfer.inlet_balance(
        LoadingFilmTransfer.LoadingAdjustedFilm(film,1.3),T,concentration)
    @test local_result.film_balance_error_mol_L_min < 1e-10
    @test_throws ArgumentError LoadingFilmTransfer.LoadingAdjustedFilm(film,0.)
end
