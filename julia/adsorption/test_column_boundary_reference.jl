using Test
include("ColumnBoundaryReference.jl")
using .ColumnBoundaryReference
@testset "Independent inlet/outlet reference" begin
    for Pe in (.001,.05,2.,20.,1000.)
        d=dirichlet_neumann_moments(length_m=1.,velocity_m_s=Pe,dispersion_m2_s=1.)
        f=danckwerts_moments(length_m=1.,velocity_m_s=Pe,dispersion_m2_s=1.)
        @test 0<d.mean_s<f.mean_s
        @test d.variance_s2>0
        @test f.mean_s==1/Pe
        @test f.variance_s2>0
    end
    d=dirichlet_neumann_moments(length_m=1.,velocity_m_s=.001,dispersion_m2_s=1.)
    @test isapprox(d.mean_s,.5;rtol=.001)
    @test isapprox(d.variance_s2,1/6;rtol=.001)
    d=dirichlet_neumann_moments(length_m=1.,velocity_m_s=1.,dispersion_m2_s=.00001)
    @test isapprox(d.mean_s,1.;rtol=2e-5)
    @test isapprox(d.variance_s2,2e-5;rtol=1e-4)
    @test_throws ArgumentError danckwerts_moments(length_m=1.,velocity_m_s=0.,dispersion_m2_s=1.)
    @test_throws ArgumentError dirichlet_neumann_moments(length_m=NaN,velocity_m_s=1.,dispersion_m2_s=1.)
end
