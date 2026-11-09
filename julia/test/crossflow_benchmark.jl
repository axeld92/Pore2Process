@testset "Crossflow source energy, target flux and flow convention" begin
    u=[.6,.9,1.4,1.8,2.2,3.1]
    j=[65.,445.,515.,555.,575.,760.]./3.6e6
    dp=[39.,69.,147.,215.,333.,637.]
    q=[15.,25.,35.,45.,55.,80.]./3.6e6
    reported=[300.,130.,330.,580.,1060.,2230.]
    for i in eachindex(u)
        a=crossflow_benchmark_point(u[i],j[i],dp[i],q[i],.03,1.,100/3.6e6,1.,1.)
        b=crossflow_benchmark_point(u[i],j[i],dp[i],q[i],.03,1.,100/3.6e6,7.,1.)
        @test abs(a.source_specific_axial_energy_J_m3/reported[i]-1)<.01
        @test a.source_specific_axial_energy_J_m3 ≈ dp[i]*q[i]/(j[i]*.03)
        @test a.nominal_admissible == (i>1)
        @test b.scenario_specific_axial_energy_J_m3 ≈ 7a.scenario_specific_axial_energy_J_m3
        @test b.source_specific_axial_energy_J_m3==a.source_specific_axial_energy_J_m3
        @test b.channel_resistance_Pa_s_m3 ≈ a.channel_resistance_Pa_s_m3/7
        @test a.module_reference_flow_m3_s/a.module_to_permeate_ratio ≈ 100/3.6e6
    end
    a=crossflow_benchmark_point(.9,445/3.6e6,69.,25/3.6e6,.03,1.,100/3.6e6,1.,1.)
    @test a.source_specific_axial_energy_J_m3 ≈ 1725/13.35
    @test a.scenario_specific_axial_energy_J_m3 ≈ 575.
    c=crossflow_recycle_load(1e4,a.permeate_reference_flow_m3_s,.5,a.module_to_permeate_ratio,a.channel_resistance_Pa_s_m3,.75,.75)
    @test c.channel_pressure_drop_Pa ≈ 69.
    @test c.circulation_shaft_power_W/c.permeate_reference_flow_m3_s ≈ 575/.75
    @test !crossflow_benchmark_point(.9,445/3.6e6,69.,25/3.6e6,.03,1.,400/3.6e6,1.,.8).nominal_admissible
    defaults=Any[.9,445/3.6e6,69.,25/3.6e6,.03,1.,100/3.6e6,1.,1.]
    for (idx,value) in ((1,NaN),(2,0.),(3,-1.),(4,Inf),(5,0.),(6,true),(7,0.),(8,-7.),(9,1.1),(8,1e308))
        args=copy(defaults);args[idx]=value
        @test_throws ArgumentError crossflow_benchmark_point(args...)
    end
end

@testset "Bounded monotonic flux constraints and transition uncertainty" begin
    vl=[.58,.97,1.37];vh=[.61,1.01,1.41]
    jl=[60.,450.,510.]./3.6e6;jh=[69.,460.,520.]./3.6e6
    @test critical_flux_bracket(vl,vh,jl,jh,.9,100/3.6e6).classification=="indeterminate"
    @test critical_flux_bracket(vl,vh,jl,jh,1.4,100/3.6e6).classification=="admissible"
    @test critical_flux_bracket(vl,vh,jl,jh,.8,500/3.6e6).classification=="inadmissible"
    @test isnothing(critical_flux_bracket(vl,vh,jl,jh,.59,100/3.6e6).lower_flux_m_s)
    @test isnothing(critical_flux_bracket(vl,vh,jl,jh,1.4,100/3.6e6).upper_flux_m_s)
    @test critical_flux_bracket(vl,vh,jl,jh,1.4,400/3.6e6,.8).classification=="indeterminate"
    for u in (.1,2.,NaN)
        @test_throws ArgumentError critical_flux_bracket(vl,vh,jl,jh,u,1e-5)
    end
    @test_throws ArgumentError critical_flux_bracket(vl[1:2],vh,jl,jh,.9,1e-5)
    @test_throws ArgumentError critical_flux_bracket(reverse(vl),vh,jl,jh,.9,1e-5)
    @test_throws ArgumentError critical_flux_bracket(vl,vh,reverse(jl),jh,.9,1e-5)
    @test_throws ArgumentError critical_flux_bracket(vl,vh,jl,jh,.9,1e-5,1.1)
end

@testset "Finite-pass passive solute conservation" begin
    for c in (5.,8.333333333333,60.), y in (.25,.5,.75), r in (0.,.5,1.)
        s=crossflow_solute_balance(.01,100/3.6e6,y,c,r)
        @test s.fresh_concentration_kg_m3>=0
        @test s.permeate_concentration_kg_m3 ≈ (1-r)*.01
        @test s.retentate_and_purge_concentration_kg_m3>=.01*(1-1e-14)
        @test abs(s.module_solute_residual_kg_s)<1e-18
        @test abs(s.mixer_solute_residual_kg_s)<1e-18
        @test abs(s.overall_solute_residual_kg_s)<1e-18
        cf=s.fresh_concentration_kg_m3;cp=s.permeate_concentration_kg_m3;cr=s.retentate_and_purge_concentration_kg_m3
        @test cf ≈ y*cp+(1-y)*cr
        @test cr ≈ .01*(c-1+r)/(c-1)
        if r==0
            @test cf ≈ .01
            @test cp ≈ .01
            @test cr ≈ .01
        elseif r==1
            @test cp==0
        end
    end
    @test crossflow_solute_balance(0.,1e-4,.5,5.,1.).fresh_concentration_kg_m3==0
    for args in ((-.01,1e-4,.5,5.,1.),(.01,1e-4,1.,5.,1.),(.01,1e-4,.5,2.,1.),(.01,1e-4,.5,5.,1.1),(.01,Inf,.5,5.,1.))
        @test_throws ArgumentError crossflow_solute_balance(args...)
    end
end
