@testset "Qualified operating-window source relation" begin
    re=[1000.,2000.,4000.]; f=[20.,50.,80.]./3.6e6
    v=critical_flux_capacity(re,f,.9f,1.1f,1500.,30/3.6e6,.95,1.)
    @test v.nominal_capacity_m_s ≈ .95*35/3.6e6
    @test v.monotonic_lower_capacity_m_s ≈ .95*18/3.6e6
    @test v.monotonic_upper_capacity_m_s ≈ .95*55/3.6e6
    @test v.nominal_admissible && !v.linear_reading_admissible
    @test v.monotonic_classification == "indeterminate"
    anchor=critical_flux_capacity(re,f,.9f,1.1f,2000.,40/3.6e6,1.,1.)
    @test anchor.monotonic_classification == "admissible"
    @test anchor.lower_source_index == anchor.upper_source_index == 2
    too_high=critical_flux_capacity(re,f,.9f,1.1f,2000.,60/3.6e6,1.,1.)
    @test too_high.monotonic_classification == "inadmissible"
    reduced=critical_flux_capacity(re,f,.9f,1.1f,1500.,30/3.6e6,.8,.8)
    @test !reduced.nominal_admissible
    @test reduced.nominal_capacity_m_s ≈ .64*35/3.6e6
    for query in (999.,4001.,NaN,Inf,-1.)
        @test_throws ArgumentError critical_flux_capacity(re,f,.9f,1.1f,query,30/3.6e6,1.,1.)
    end
    @test_throws ArgumentError critical_flux_capacity(re,f,.9f,1.1f,1500.,NaN,1.,1.)
    @test_throws ArgumentError critical_flux_capacity(re,f,.9f,1.1f,1500.,30/3.6e6,1.01,1.)
    @test_throws ArgumentError critical_flux_capacity(re,f,1.01f,1.1f,1500.,30/3.6e6,1.,1.)
    @test_throws ArgumentError critical_flux_capacity(reverse(re),f,.9f,1.1f,1500.,30/3.6e6,1.,1.)
end

@testset "Axial momentum screen, independent circular profile" begin
    D=.0023;mu=.001;rho=1000.;J=3e-6
    A=pi*D^2/4;P=pi*D;C=pi*D^4/128
    # Exact parabolic circular profile has beta=4/3. Its momentum/friction
    # ratio reduces independently to beta*rho*D*J/(4mu).
    exact=(4/3)*rho*D*J/(4mu)
    @test axial_momentum_screen(C,A,P,mu,rho,J,4/3) ≈ exact
    @test axial_momentum_screen(C,A,P,mu,rho,J,2.) ≈ 1.5exact
    @test axial_momentum_screen(C,A,P,2mu,rho,J,4/3) ≈ exact/2
    @test_throws ArgumentError axial_momentum_screen(C,A,P,mu,rho,J,.9)
    @test_throws ArgumentError axial_momentum_screen(C,A,P,mu,rho,NaN,2.)
end
