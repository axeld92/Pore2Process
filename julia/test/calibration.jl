@testset "Support closure, analytic log fit and held-out diagnostics" begin
    # Independent hand-computable CK case: d=6e-7, e=1/2 => k=d²/360.
    d,e = 6e-7,0.5
    k = d^2/360
    @test support_permeability(d,e) ≈ k rtol=1e-12
    @test support_permeability(2d,e) ≈ 4k rtol=1e-12
    @test support_permeability(d,e;multiplier=3.) ≈ 3k rtol=1e-12
    # Ratios 2 and 8 give geometric mean 4 (not arithmetic mean 5).
    fit=fit_support_multiplier([d,d],[e,e],[2k,8k])
    @test fit.multiplier ≈ 4. rtol=1e-12
    @test fit.calibration_log_rmse ≈ log(2.) rtol=1e-12
    @test fit.calibration_count == 2
    @test all(x -> isapprox(x,4k;rtol=1e-12),fit.predictions_m2)
    @test fit_support_multiplier([d],[e],[7k]).multiplier ≈ 7. rtol=1e-12
    exact=fit_support_multiplier([d,2d],[e,e],[3k,12k])
    @test exact.multiplier ≈ 3. rtol=1e-12
    @test exact.calibration_log_rmse < 1e-12
    err=permeability_errors([2k,k/2],[k,k])
    @test err.mean_absolute_relative_error ≈ 0.75
    @test err.max_absolute_relative_error ≈ 1.
    @test err.mean_relative_error ≈ 0.25
    @test err.log_rmse ≈ log(2.)
    @test err.count == 2
    for (d,e,b) in ((0.,.3,1.),(-1.,.3,1.),(Inf,.3,1.),(1e-6,0.,1.),(1e-6,1.,1.),
                    (1e-6,NaN,1.),(1e-6,true,1.),(1e-6,.3,-1.),(1e-6,.3,Inf),
                    (1e-6,.3,true),(1e300,.3,1.),(1e-300,.3,1.))
        @test_throws ArgumentError support_permeability(d,e;multiplier=b)
    end
    for (ds,es,ks) in (([],[],[]),([d],[e,e],[k]),([d],[e],[0.]),([d],[e],[NaN]))
        @test_throws ArgumentError fit_support_multiplier(ds,es,ks)
    end
    @test_throws ArgumentError permeability_errors([],[])
    @test_throws ArgumentError permeability_errors([k],[k,k])
    @test_throws ArgumentError permeability_errors([k],[Inf])
    @test_throws ArgumentError permeability_errors([1e300],[1e-300])
end
