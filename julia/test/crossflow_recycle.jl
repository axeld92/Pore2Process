@testset "Crossflow recycle: independent flow, work and scaling reference" begin
    tmp,Qp,y,C,R,ef,ec=1e4,1e-4,.5,5.,1e8,.75,.8
    a=crossflow_recycle_load(tmp,Qp,y,C,R,ef,ec)
    @test a.fresh_feed_reference_flow_m3_s ≈ 2Qp
    @test a.module_reference_flow_m3_s ≈ 5Qp
    @test a.recycle_reference_flow_m3_s ≈ 3Qp
    @test a.purge_reference_flow_m3_s ≈ Qp
    @test a.membrane_permeate_split_fraction ≈ .2
    @test a.purge_retentate_split_fraction ≈ .25
    @test a.loop_return_fraction ≈ .6
    @test a.channel_pressure_drop_Pa ≈ 5Qp*R
    @test a.feed_shaft_power_W ≈ 2Qp*tmp/ef
    @test a.circulation_shaft_power_W ≈ (5Qp)^2*R/ec
    for c in (5.,10.,15.), recovery in (.25,.5,.75)
        s=crossflow_recycle_load(tmp,Qp,recovery,c,R,ef,ec)
        @test s.fresh_feed_reference_flow_m3_s ≈ s.permeate_reference_flow_m3_s+s.purge_reference_flow_m3_s
        @test s.module_reference_flow_m3_s ≈ s.fresh_feed_reference_flow_m3_s+s.recycle_reference_flow_m3_s
        @test s.module_retentate_reference_flow_m3_s ≈ s.purge_reference_flow_m3_s+s.recycle_reference_flow_m3_s
        @test s.module_reference_flow_m3_s*s.membrane_permeate_split_fraction ≈ Qp
        @test s.module_retentate_reference_flow_m3_s*s.purge_retentate_split_fraction ≈ s.purge_reference_flow_m3_s
        @test s.total_shaft_power_W ≈ s.membrane_hydraulic_dissipation_W+s.purge_throttle_dissipation_W+s.channel_hydraulic_dissipation_W+s.pump_loss_W
        # Geometric-series recurrence m_next = beta*m + beta*m_f.
        @test s.loop_return_fraction*s.fresh_feed_reference_flow_m3_s/(1-s.loop_return_fraction) ≈ s.recycle_reference_flow_m3_s
    end
    twice=crossflow_recycle_load(tmp,Qp,y,2C,R,ef,ec)
    @test twice.circulation_shaft_power_W ≈ 4a.circulation_shaft_power_W
    @test twice.feed_shaft_power_W ≈ a.feed_shaft_power_W
    resistance=crossflow_recycle_load(tmp,Qp,y,C,2R,ef,ec)
    @test resistance.circulation_shaft_power_W ≈ 2a.circulation_shaft_power_W
    changed=crossflow_recycle_load(.8tmp,Qp,y,C,R,ef,ec)
    @test changed.circulation_shaft_power_W==a.circulation_shaft_power_W
    @test a.total_shaft_power_W-changed.total_shaft_power_W ≈ .2a.feed_shaft_power_W
    no_loss=crossflow_recycle_load(tmp,Qp,y,C,0.,ef,ec)
    @test no_loss.circulation_shaft_power_W==0
    @test no_loss.total_shaft_power_W ≈ once_through_filtration(tmp,Qp,y,ef).shaft_power_W
    @test a.specific_shaft_energy_kWh_m3 ≈ a.total_shaft_power_W/Qp/3.6e6
    for args in ((0.,Qp,y,C,R,ef,ec),(tmp,0.,y,C,R,ef,ec),(tmp,Qp,1.,C,R,ef,ec),
                 (tmp,Qp,0.,C,R,ef,ec),(tmp,Qp,y,2.,R,ef,ec),(tmp,Qp,y,1.,R,ef,ec),
                 (tmp,Qp,y,C,-R,ef,ec),(tmp,Qp,y,C,R,0.,ec),(tmp,Qp,y,C,R,ef,1.1),
                 (NaN,Qp,y,C,R,ef,ec),(tmp,Qp,y,Inf,R,ef,ec),(tmp,Qp,true,C,R,ef,ec),
                 (tmp,1.,y,C,1e308,ef,ec),(1e300,1e-300,1e-100,1e101,R,ef,ec),
                 (tmp,1e308,y,C,R,ef,ec))
        @test_throws ArgumentError crossflow_recycle_load(args...)
    end
end
