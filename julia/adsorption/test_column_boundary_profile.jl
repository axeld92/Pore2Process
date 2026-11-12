using Test
include("ColumnBoundaryReference.jl")
include("TracerContact.jl")
using .ColumnBoundaryReference
@testset "Conservative outlet profile reference" begin
    p=(length_m=1.,velocity_m_s=1.,dispersion_m2_s=.05)
    r=solve_dirichlet_neumann(;p...,t_end_s=12.,dt_s=.00025,cells=400)
    exact=dirichlet_neumann_moments(;p...)
    moments=TracerContact.step_moments(r.time_s[2:end],r.outlet_flux_per_area_m_s[2:end])
    @test isapprox(moments.mean_s,exact.mean_s;rtol=.001)
    @test isapprox(moments.variance_s2,exact.variance_s2;rtol=.01)
    @test r.normalized_mass_residual<1e-8
    @test all(0 .<=r.final_concentration .<=1+1e-10)
    @test isapprox(r.outlet_flux_per_area_m_s[end],1.;rtol=1e-8)
end
