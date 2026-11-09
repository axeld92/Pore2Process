@testset "Paired measured permeance, series impedance and energy" begin
    # Independent oracle: support Z 1e8 and fixed Z 2e8 -> Ptotal=1/3e8.
    fit=series_impedance_fit([1e-8,5e-9],[1/3e8,1/4e8])
    @test fit.fixed_impedance_Pa_s_m ≈ 2e8
    @test fit.training_count==2
    @test fit.limiting_permeance_m_Pa_s ≈ 1/2e8
    @test fit.impedance_rms_residual_Pa_s_m < 1e-7
    # Arithmetic residual fit uses the mean of per-case impedances.
    @test series_impedance_fit([1e-8,1e-8],[1/2e8,1/4e8]).fixed_impedance_Pa_s_m ≈ 2e8
    s=series_flux(1e-8,2e8,1e5)
    @test s.total_permeance_m_Pa_s ≈ 1/3e8
    @test s.flux_m_s ≈ 1/3000
    @test s.non_support_fraction ≈ 2/3
    @test s.support_fraction ≈ 1/3
    @test series_flux(1e-8,2e8,-1e5).flux_m_s ≈ -1/3000
    @test series_flux(1e-8,2e8,0.).flux_m_s==0.
    @test series_flux(1e-8,0.,1e5).total_permeance_m_Pa_s ≈ 1e-8
    load=impedance_load(1e-8,2e8,1e-5,2.,.5)
    @test load.pressure_drop_Pa ≈ 3000.
    @test load.flow_m3_s ≈ 2e-5
    @test load.shaft_power_W ≈ .12
    @test load.specific_shaft_energy_kWh_m3 ≈ 1/600
    @test impedance_load(1e-8,2e8,0.,2.,.5).shaft_power_W==0.
    @test impedance_load(1e-8,2e8,2e-5,2.,.5).shaft_power_W ≈ 4load.shaft_power_W
    for (ps,pt) in (([],[]),([1e-8],[]),([1e-8],[1e-8]),([1e-8],[2e-8]),([Inf],[1e-8]),([1e-8],[NaN]))
        @test_throws ArgumentError series_impedance_fit(ps,pt)
    end
    for args in ((0.,1.,1.),(1e-8,-1.,1.),(1e-8,Inf,1.),(1e-8,1.,NaN),(true,1.,1.))
        @test_throws ArgumentError series_flux(args...)
    end
    for args in ((1e-8,2e8,-1.,1.,.5),(1e-8,2e8,1.,0.,.5),(1e-8,2e8,1.,1.,1.1))
        @test_throws ArgumentError impedance_load(args...)
    end
end
