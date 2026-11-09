@testset "Fixed-flux layers in series: independent resistance and energy oracle" begin
    # Independent dimensional example: Rs=2e10, Rl=3e10, Rf=1e10 m^-1.
    k,L,Rl,mu,J,A,eta = 5e-14,1e-3,3e10,1e-3,1e-5,2.,0.5
    f=filtration_load(k,L,Rl,mu,J,A,eta;fixed_resistance_m_inv=1e10)
    @test f.support_resistance_m_inv ≈ 2e10
    @test f.total_resistance_m_inv ≈ 6e10
    @test f.support_pressure_drop_Pa ≈ 200.
    @test f.selective_pressure_drop_Pa ≈ 300.
    @test f.fixed_pressure_drop_Pa ≈ 100.
    @test f.total_pressure_drop_Pa ≈ 600.
    @test f.flow_m3_s ≈ 2e-5
    @test f.shaft_power_W ≈ 0.024
    @test f.specific_shaft_energy_kWh_m3 ≈ 1/3000
    @test f.selective_resistance_fraction ≈ 0.5
    @test !f.selective_dominates
    @test !f.selective_at_least_90_percent
    @test f.support_pressure_drop_Pa+f.selective_pressure_drop_Pa+f.fixed_pressure_drop_Pa ≈ f.total_pressure_drop_Pa
    @test f.shaft_power_W/(f.flow_m3_s*3.6e6) ≈ f.specific_shaft_energy_kWh_m3
    # At fixed selective layer, a fourfold support-k improvement removes 3/4
    # of the support contribution, not 3/4 of the total energy.
    for alpha in (0.,0.1,1.,10.,100.)
        base=filtration_load(k,L,alpha*L/k,mu,J,A,eta)
        better=filtration_load(4k,L,alpha*L/k,mu,J,A,eta)
        @test 1-better.shaft_power_W/base.shaft_power_W ≈ 0.75/(1+alpha)
        @test better.selective_resistance_m_inv == base.selective_resistance_m_inv
        @test better.flow_m3_s == base.flow_m3_s
        @test better.shaft_power_W < base.shaft_power_W
    end
    for (alpha,majority,ninety) in ((1.,false,false),(1.01,true,false),(8.99,true,false),(9.,true,true))
        v=filtration_load(k,L,alpha*L/k,mu,J,A,eta)
        @test v.selective_dominates == majority
        @test v.selective_at_least_90_percent == ninety
    end
    z=filtration_load(k,L,Rl,mu,0.,A,eta)
    @test z.flow_m3_s == z.total_pressure_drop_Pa == z.shaft_power_W == z.specific_shaft_energy_kWh_m3 == 0.
    doubled=filtration_load(k,L,Rl,mu,2J,A,eta;fixed_resistance_m_inv=1e10)
    @test doubled.total_pressure_drop_Pa ≈ 2f.total_pressure_drop_Pa
    @test doubled.shaft_power_W ≈ 4f.shaft_power_W
    @test doubled.specific_shaft_energy_kWh_m3 ≈ 2f.specific_shaft_energy_kWh_m3
    @test permeance_LMH_bar_to_SI(360.) ≈ 1e-9
    @test resistance_from_permeance(1e-9,1e-3) ≈ 1e12
    @test_throws ArgumentError permeance_LMH_bar_to_SI(0.)
    @test_throws ArgumentError resistance_from_permeance(Inf,1e-3)
    defaults=Any[k,L,Rl,mu,J,A,eta]
    for (index,value) in ((1,0.),(1,nothing),(2,-1.),(3,-1.),(4,0.),(5,-1.),(6,0.),(7,0.),(7,1.1),(1,NaN),(3,Inf),(5,true),(7,true))
        bad=copy(defaults); bad[index]=value
        @test_throws ArgumentError filtration_load(bad...)
    end
    @test_throws ArgumentError filtration_load(defaults...;fixed_resistance_m_inv=-1.)
    @test_throws ArgumentError filtration_load(defaults...;fixed_resistance_m_inv=NaN)
end
