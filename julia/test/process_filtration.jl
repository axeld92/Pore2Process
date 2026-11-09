@testset "Once-through filtration: independent hydraulic balances" begin
    tmp,Qp,eta=1e4,1e-4,.75
    f=once_through_filtration(tmp,Qp,.5,eta)
    @test f.feed_reference_flow_m3_s==2Qp
    @test f.retentate_reference_flow_m3_s==Qp
    @test f.shaft_power_W ≈ 2Qp*tmp/eta
    @test f.membrane_dissipation_W ≈ Qp*tmp
    @test f.retentate_throttle_dissipation_W ≈ Qp*tmp
    @test f.pump_loss_W ≈ f.shaft_power_W*(1-eta)
    for y in (.25,.5,.75), c in (0.,2e4)
        s=once_through_filtration(tmp,Qp,y,eta,c)
        @test s.feed_reference_flow_m3_s ≈ s.permeate_reference_flow_m3_s+s.retentate_reference_flow_m3_s
        @test s.hydraulic_power_W ≈ s.membrane_dissipation_W+s.retentate_throttle_dissipation_W+s.channel_dissipation_W
        @test s.shaft_power_W ≈ s.hydraulic_power_W+s.pump_loss_W
        @test s.specific_shaft_energy_kWh_m3 ≈ (tmp+c)/(eta*y*3.6e6)
        better=once_through_filtration(.8tmp,Qp,y,eta,c)
        @test 1-better.shaft_power_W/s.shaft_power_W ≈ .2tmp/(tmp+c)
    end
    @test once_through_filtration(tmp,Qp,.25,eta).shaft_power_W ≈ 2f.shaft_power_W
    @test once_through_filtration(tmp,2Qp,.5,eta).shaft_power_W ≈ 2f.shaft_power_W
    for args in ((0.,Qp,.5,eta,0.),(tmp,0.,.5,eta,0.),(tmp,Qp,0.,eta,0.),
                 (tmp,Qp,1.,eta,0.),(tmp,Qp,1.1,eta,0.),(tmp,Qp,.5,1.1,0.),
                 (tmp,Qp,.5,eta,-1.),(NaN,Qp,.5,eta,0.),(tmp,Qp,.5,eta,Inf),
                 (tmp,Qp,true,eta,0.),(1e308,Qp,.5,eta,1e308),
                 (1e300,1e-300,1e-100,1.,0.))
        @test_throws ArgumentError once_through_filtration(args...)
    end
end
