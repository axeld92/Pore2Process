@testset "Independent polygon duct reference" begin
    v=[0. 0.; 1. 0.; 1. 1.; 0. 1.]
    exact=(1-192/pi^5*sum(tanh(n*pi/2)/n^5 for n in 1:2:501))/12
    a=duct_stokes_reference(v,32); b=duct_stokes_reference(v,64)
    @test abs(b.conductance_geometry_m4/exact-1)<.002
    @test abs(b.conductance_geometry_m4/exact-1)<abs(a.conductance_geometry_m4/exact-1)
    @test b.residual<1e-9
    @test_throws ArgumentError duct_stokes_reference(v,2)
end

@testset "Distributed star membrane and biological tracer" begin
    C=2.67e-12; A=9.66e-6; P=.015; L=.3; mu=.001; rho=1000.; Re=1388.; J=15/3.6e6
    zs=3.6e8; zf=8e8
    a=star_membrane_load(C,A,P,L,mu,rho,Re,J,1.,zs,zf)
    @test a.channel_inlet_flow_m3_s-a.channel_outlet_flow_m3_s≈a.channel_permeate_flow_m3_s
    @test a.pressure_Pa[1]≈a.inlet_TMP_Pa
    @test a.pressure_Pa[end]≈a.outlet_TMP_Pa
    @test a.channel_flow_m3_s[end]≈a.channel_outlet_flow_m3_s rtol=1e-12
    @test a.inlet_TMP_Pa>a.mean_TMP_Pa>a.outlet_TMP_Pa>0
    @test maximum(a.local_flux_m_s)>J>minimum(a.local_flux_m_s)
    dx=L/100
    meanJ=sum((a.local_flux_m_s[1:end-1]+a.local_flux_m_s[2:end])/2)*dx/L
    @test meanJ≈J rtol=1e-5
    @test maximum(abs.(diff(a.pressure_Pa)/dx+mu/C*(a.channel_flow_m3_s[1:end-1]+a.channel_flow_m3_s[2:end])/2))/a.inlet_TMP_Pa<1e-5
    @test maximum(abs.(diff(a.channel_flow_m3_s)/dx+P/(zs+zf)*(a.pressure_Pa[1:end-1]+a.pressure_Pa[2:end])/2))/a.channel_inlet_flow_m3_s<1e-5
    b=star_membrane_load(C,A,P,L,mu,rho,Re,2J,1.,zs,zf)
    @test b.mean_TMP_Pa≈2a.mean_TMP_Pa
    @test b.channel_drop_Pa<a.channel_drop_Pa
    for rejection in (0.,.5,.65,1.)
        s=star_feed_balance(.0687,rejection,a,.5)
        @test all(abs.([s.module_residual_kg_s,s.mixer_residual_kg_s,s.overall_residual_kg_s]).<1e-18)
        @test s.permeate_COD_kg_m3>=-1e-12
        @test s.retentate_COD_kg_m3+1e-14>=s.module_inlet_COD_kg_m3>=s.fresh_COD_kg_m3-1e-14
        if rejection==0.; @test s.permeate_COD_kg_m3≈.0687; end
        if rejection==1.; @test abs(s.permeate_COD_kg_m3)<1e-12; end
    end
    @test_throws ArgumentError star_feed_balance(.0687,1.1,a,.5)
    @test_throws ArgumentError star_membrane_load(-C,A,P,L,mu,rho,Re,J,1.,zs,zf)
    @test_throws ArgumentError star_membrane_load(C,A,P,L,mu,rho,Re,1000J,1.,zs,zf)
end
