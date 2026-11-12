using Test
include("FilmParticle.jl")
include("ExchangeColumn.jl")

@testset "Conservative exchange and film-particle coupling" begin
    feed = (;length_m=.10,area_m2=pi*.011^2/4,porosity=.427,superficial_velocity_m_s=.00768,
            inlet_concentration_kg_m3=.15,dt_s=.05,t_end_s=20.,radial_cells=12,
            diffusivity_m2_s=1.846e-12,film_m_s=4e-5)
    original = FilmParticle.solve_column(;feed...,axial_cells=4)
    coupled = ExchangeColumn.solve_column([1.],[1.],ones(1,1);feed...,periods=2,subcells=2)
    @test coupled.outlet_concentration_kg_m3 ≈ original.outlet_concentration_kg_m3 rtol=1e-8 atol=1e-13
    @test coupled.adsorbed_inventory_kg ≈ original.adsorbed_inventory_kg rtol=1e-8
    @test coupled.liquid_inventory_kg ≈ original.liquid_inventory_kg rtol=1e-8
    @test coupled.normalized_mass_residual < 1e-9
    @test coupled.normalized_interface_flux_mismatch < 1e-9
    phased = ExchangeColumn.solve_column([1.],[1.],ones(1,1);feed...,periods=2,subcells=2,
        startup_s=5.,startup_dt_s=.05)
    @test phased.outlet_concentration_kg_m3 == coupled.outlet_concentration_kg_m3

    # One implicit time step with two nonadsorbing paths and a prescribed
    # exchange matrix has a direct finite-volume solution.
    mix = [.8 .2; .2 .8]
    tracer = ExchangeColumn.solve_column([.5,.5],[.5,1.5],mix;
        feed...,periods=2,subcells=1,dt_s=.5,t_end_s=.5,diffusivity_m2_s=0.)
    tau = tracer.nominal_liquid_transit_s
    factors = .5 ./ (.5 .+ (tau/2).*[.5,1.5])
    expected = sum(factors.*(mix*(.15*factors)))/2
    @test tracer.outlet_concentration_kg_m3[end] ≈ expected rtol=1e-13
    @test tracer.adsorbed_inventory_kg[end]==0
    @test tracer.normalized_mass_residual < 1e-12
    @test_throws ArgumentError ExchangeColumn.solve_column([.5,.5],[.5,1.5],[1. 0.; .2 .8];feed...,periods=2)
end
