using Test
include("../../pentyl_propanoate_transfer/published/model.jl")
include("model.jl")
using .LorenzoPackedBed, .BeadFilm

@testset "Film dimensions and local reaction balance" begin
    u, d, rho, mu, D, porosity = 0.0001, 0.000593, 790.0, 0.001449, 6e-10, 0.64
    wg = film_coefficients(u, d, rho, mu, D, porosity, "wilson_geankoplis")
    @test wg.reynolds ≈ rho * u * d / mu
    @test wg.schmidt ≈ mu / (rho * D)
    @test wg.reynolds * wg.schmidt ≈ u * d / D
    @test wg.sherwood ≈ 1.09 / porosity * cbrt(u * d / D)
    @test wg.coefficient_m_s ≈ wg.sherwood * D / d
    @test film_coefficients(u, d, rho, mu, D, porosity, "isolated_sphere").coefficient_m_s ≈ 2D / d
    @test_throws ArgumentError film_coefficients(u, d, rho, mu, 0, porosity, "wilson_geankoplis")
    @test_throws ArgumentError film_coefficients(u, d, rho, mu, D, porosity, "other")
    rate(x) = 2 * (1 - x)
    @test surface_progress(rate, 0.3, 1.0, 8.0) ≈ (2 + 8 * 0.3) / 10 atol=1e-12
    @test surface_progress(rate, 1.1, 1.0, 8.0) ≈ (2 + 8 * 1.1) / 10 atol=1e-12
end

@testset "Constant catalyst inventory and film limits" begin
    kinetics = PingPong([3520.6, 1746.6, 260.76, 1913.3, 42250, 1.2764e-6, 103740, 776560,
                        848920, 4.4007e-6, 16898, 12813, 6.6093e-6, 12.292, 3.3629, -28091, 313.15])
    temperature, residence, loading, eta, pe = 60.0, 6.0, 333.0, 0.869, 115.0
    base = packed_bed(kinetics, temperature, residence, loading, eta, pe; cells=256)
    fast = packed_bed(FilmKinetics(kinetics, 100_000, eta), temperature, residence, loading, eta, pe; cells=256)
    limited = FilmKinetics(kinetics, 2.0, eta)
    slow = packed_bed(limited, temperature, residence, loading, eta, pe; cells=256)
    @test abs(fast.conversion - base.conversion) < 1e-5
    @test 0 < slow.conversion < fast.conversion < base.conversion
    @test slow.balance_error < 1e-9
    balance = inlet_balance(limited, temperature, loading)
    @test 0 < balance.surface_acid_relative_to_bulk < 1
    @test 0 < balance.rate_relative_to_no_film < 1
    @test balance.film_balance_error_mol_L_min < 1e-12
end
