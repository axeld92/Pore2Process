using Test
include("AdsorptionColumn.jl")
include("OperatingPointHydraulics.jl")
using .OperatingPointHydraulics
@testset "Native operating-point translation" begin
    p=(native_gradient_Pa_m=100.,native_velocity_m_s=.001,
        native_viscosity_Pa_s=.001,native_density_kg_m3=1000.,
        carrier_viscosity_Pa_s=.0016,carrier_density_kg_m3=800.,
        superficial_velocity_m_s=.002,length_m=.5,area_m2=.01,efficiency_fraction=.75)
    r=translate(;p...)
    @test r.native_to_carrier_velocity_scale==2.
    @test r.native_to_carrier_pressure_scale==3.2
    @test isapprox(r.operating_point_apparent_permeability_m2,1e-8;rtol=1e-14)
    @test isapprox(r.pressure_drop_Pa,160.;rtol=1e-14)
    @test isapprox(r.flow_m3_s,2e-5;rtol=1e-14)
    @test isapprox(r.reference_shaft_power_W,160*2e-5/.75;rtol=1e-14)
    @test r.Reynolds_relative_mismatch==0.
    @test_throws ArgumentError translate(;merge(p,(superficial_velocity_m_s=.003,))...)
    @test_throws ArgumentError translate(;merge(p,(carrier_viscosity_Pa_s=NaN,))...)
    @test_throws ArgumentError translate(;merge(p,(efficiency_fraction=1.1,))...)
end
