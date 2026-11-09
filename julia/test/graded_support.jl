@testset "Graded support independent limits and scale boundary" begin
    r,L,mu,A=4e-6,40e-6,1e-3,20e-6^2
    s=graded_pore_reference([0.,L],[r,r],mu,A)
    @test s.hydraulic_resistance_Pa_s_m3 ≈ 8mu*L/(pi*r^4) rtol=1e-14
    @test s.permeability_cell_m2 ≈ pi*r^4/(8A) rtol=1e-14
    @test s.volume_porosity ≈ pi*r^2/A rtol=1e-14
    # r increases linearly from r to 2r: integral = 7L/(24r^4).
    t=graded_pore_reference([0.,L],[r,2r],mu,A)
    @test t.hydraulic_resistance_Pa_s_m3/s.hydraulic_resistance_Pa_s_m3 ≈ 7/24
    rev=graded_pore_reference([0.,L],[2r,r],mu,A)
    @test rev==t
    split=graded_pore_reference([0.,L/2,L],[r,1.5r,2r],mu,A)
    @test split.hydraulic_resistance_Pa_s_m3 ≈ t.hydraulic_resistance_Pa_s_m3
    base=graded_support_load(8e14,8e14,1e-8,2e8,1e-5,2.,.5)
    better=graded_support_load(4e14,8e14,1e-8,2e8,1e-5,2.,.5)
    @test base.support_permeance_m_Pa_s==1e-8
    @test better.support_permeance_m_Pa_s==2e-8
    @test better.non_support_impedance_Pa_s_m==base.non_support_impedance_Pa_s_m
    @test better.pressure_drop_Pa/base.pressure_drop_Pa ≈ 2.5/3
    @test better.shaft_power_W/base.shaft_power_W ≈ 2.5/3
    @test better.non_support_fraction ≈ .8
    @test better.flow_m3_s==2e-5
    @test better.specific_shaft_energy_kWh_m3 ≈ better.pressure_drop_Pa/.5/3.6e6
    @test graded_support_load(4e14,8e14,1e-8,2e8,0.,2.,.5).shaft_power_W==0
    for (x,radii,mu0,a) in (([0.,0.],[r,r],mu,A),([0.,L],[r],mu,A),
            ([0.,Inf],[r,r],mu,A),([0.,L],[r,-r],mu,A),([0.,L],[r,r],0.,A),
            ([0.,L],[r,r],mu,r^2))
        @test_throws ArgumentError graded_pore_reference(x,radii,mu0,a)
    end
    for args in ((0.,1.,1e-8,2e8,1e-5,1.,.75), (1.,-1.,1e-8,2e8,1e-5,1.,.75),
                 (1.,1.,NaN,2e8,1e-5,1.,.75),(1.,1.,1e-8,-2e8,1e-5,1.,.75))
        @test_throws ArgumentError graded_support_load(args...)
    end
end
