using Test
include("../../../julia/adsorption/FilmParticle.jl")
include("../../../julia/adsorption/ExchangeColumn.jl")
include("../../../julia/adsorption/SharedContactColumn.jl")
include("model.jl")
using .MemoryColumn: delays, solve_column

@testset "Pulse moments and conservative reactive column" begin
    T = ones(1,1)
    flow = [1.]
    p = delays(T,fill(1.,1,1),fill(1.16,1,1),flow)
    @test sum(p.share) ≈ 1.
    @test sum(p.share.*p.seconds) ≈ 1.
    @test sum(p.share.*(p.seconds.-1).^2) ≈ .16
    @test minimum(p.seconds)>0
    skew = delays(T,T,2T,flow)
    @test skew.seconds ≈ [.5,3.]
    @test skew.share ≈ [.8,.2]
    @test sum(skew.share.*skew.seconds) ≈ 1.
    @test sum(skew.share.*(skew.seconds.-1).^2) ≈ 1.
    inputs = (blocks=2,liquid_volume_m3=2e-6,dry_mass_kg=.001,flow_m3_s=1e-6,
              inlet_concentration_kg_m3=.15,t_end_s=5.,dt_s=.1)
    run = solve_column([p],flow,T;inputs...,film_m_s=0.)
    # Two independent delays of 0.6 and 1.4 s with equal probability.
    mid = run.time_s[2:end].-.05
    cdf = (mid.>1.2).*.25 + (mid.>2.).*.5 + (mid.>2.8).*.25
    @test run.outlet_concentration_kg_m3[2:end] ./ .15 ≈ cdf atol=1e-12
    @test run.normalized_mass_residual < 1e-12
    @test run.liquid_inventory_kg[end] ≈ 2e-6*.15
    r = solve_column([p],flow,T;inputs...,diffusivity_m2_s=2e-12,film_m_s=5e-5)
    @test r.adsorbed_inventory_kg[end]>0
    @test r.adsorbed_inventory_kg[end] ≈ r.cumulative_film_uptake_kg-r.cumulative_film_release_kg rtol=1e-10
    @test r.normalized_mass_residual < 1e-10
    @test r.normalized_interface_flux_mismatch < 1e-10
    @test r.physical_adsorbent_mass_kg == .001
end
