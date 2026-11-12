using Test
include("TracerContact.jl")
using .TracerContact
@testset "Independent absorbing-end tracer reference" begin
    for Pe in (2.,10.,40.,100.)
        L=.008;v=.0005;D=v*L/Pe
        exact=absorbing_moments(length_m=L,velocity_m_s=v,dispersion_m2_s=D)
        r=solve_ade(length_m=L,velocity_m_s=v,dispersion_m2_s=D,t_end_s=100.,dt_s=.001,cells=800)
        m=step_moments(r.time_s,r.outlet_flux_per_area_m_s)
        @test isapprox(m.mean_s,exact.mean_s;rtol=.002)
        @test isapprox(m.variance_s2,exact.variance_s2;rtol=.01)
        @test isapprox(r.outlet_flux_per_area_m_s[end],exact.steady_flux_per_area_m_s;rtol=.002)
        @test r.normalized_mass_residual<1e-9
        @test minimum(r.final_concentration)>=0
        @test maximum(r.final_concentration)<=1
        c=derive_contact(length_m=L,mean_s=exact.mean_s,variance_s2=exact.variance_s2)
        @test isapprox(c.effective_velocity_m_s,v;rtol=1e-8)
        @test isapprox(c.effective_dispersion_m2_s,D;rtol=1e-8)
    end
    @test_throws ArgumentError derive_contact(length_m=1.,mean_s=1.,variance_s2=.5)
    @test_throws ArgumentError solve_ade(length_m=1.,velocity_m_s=1.,dispersion_m2_s=.00001,t_end_s=1.)
    @test_throws ArgumentError step_moments([0.,1.],[0.,1.])
    @test_throws ArgumentError absorbing_moments(length_m=NaN,velocity_m_s=1.,dispersion_m2_s=1.)
end
