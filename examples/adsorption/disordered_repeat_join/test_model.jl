using Test
include("../../../julia/adsorption/FilmParticle.jl")
include("../../../julia/adsorption/ExchangeColumn.jl")
include("../../../julia/adsorption/SharedContactColumn.jl")
include("../disordered_memory_column/model.jl")
include("model.jl")
using .MemoryColumn: Delays, solve_column
using .RepeatJoin: path_preserving, solve_preserved, repeat_moments

@testset "Repeat joins conserve flow and change only routing" begin
    flow = [.3,.7]
    first = Delays([1,1],[1,2],flow,[.5,1.5])
    rules = [first]
    area = reshape(flow,2,1)
    p = (blocks=2,liquid_volume_m3=2.4e-6,dry_mass_kg=.001,flow_m3_s=1e-6,
         inlet_concentration_kg_m3=1.,film_m_s=0.,t_end_s=5.,dt_s=.1,seam_history=true)
    moved = path_preserving(rules,2)
    @test moved[1] === first
    @test moved[2].donor == first.receiver
    @test moved[2].seconds == first.seconds
    @test moved[2].share == first.share
    mixed = solve_column(rules,flow,area;p...)
    persisted = solve_preserved(rules,flow,area;p...)
    @test mixed.seam_concentration[1,:] == persisted.seam_concentration[1,:]
    mid = mixed.time_s[2:end].-.05
    @test mixed.outlet_concentration_kg_m3[2:end] ≈
        .09 .* (mid .> 1.) + .42 .* (mid .> 2.) + .49 .* (mid .> 3.) atol=1e-12
    @test persisted.outlet_concentration_kg_m3[2:end] ≈ .3 .* (mid .> 1.) + .7 .* (mid .> 3.) atol=1e-12
    for (keep,var,r) in ((false,.42,mixed),(true,.84,persisted))
        m = repeat_moments(rules,flow,2,2.4,keep)
        @test m.mean ≈ 2.4
        @test m.variance ≈ var
        @test r.physical_adsorbent_mass_kg == p.dry_mass_kg
        @test r.liquid_inventory_kg[end] ≈ p.liquid_volume_m3
        @test r.normalized_mass_residual < 1e-12
    end
    @test_throws ArgumentError path_preserving(rules,0)
end
