using Test
include("model.jl")
using .PentylPropanoate

@testset "Apparent reversible batch kinetics" begin
    known = Kinetics(log(0.18), 0.45, log(5.0), 0.2)
    times = [0.0, 5.0, 30.0, 60.0, 150.0, 300.0]
    loading = fill(10.0, length(times))
    temperatures = repeat([40.0, 60.0, 80.0], inner=length(times))
    observed = batch_conversion.(Ref(known), temperatures, repeat(loading, 3), repeat(times, 3))
    pooled = fit_batch(temperatures, repeat(loading, 3), repeat(times, 3), observed)
    @test coefficients(pooled.kinetics) ≈ coefficients(known) rtol=2e-7
    single = fit_at_temperature(60.0, loading, times, observed[7:12], 5.0)
    @test rate_constant(single.kinetics, 60.0) ≈ 0.18 rtol=1e-7
    @test batch_conversion(known, 60, 10, 0) == 0
    @test batch_conversion(known, 60, 10, 100_000) ≈ equilibrium_conversion(known, 60)
    @test_throws ArgumentError fit_at_temperature(60, loading, times, observed, 5.0)
end

@testset "Fixed inventory in plug and dispersed flow" begin
    kinetics = Kinetics(log(0.18), 0, log(5.0), 0)
    mass_kg, acid_mol_s = 0.006, 3.3 / 60_000
    exposure = mass_kg / acid_mol_s
    equivalent_batch_min = exposure / (0.1 * 0.0740785 * 60)
    plug = bed_prediction(kinetics, 60, mass_kg, acid_mol_s).conversion
    @test plug ≈ batch_conversion(kinetics, 60, 10, equivalent_batch_min)

    # Integrate the four-species stoichiometric rate written in mole fractions,
    # independently of the solver's closed-form conversion formula.
    rhs(x) = 0.18 * ((1 - x)/3 * (2 - x)/3 - (x/3)^2/5)
    dx = exposure / 1024
    x = 0.0
    for _ in 1:1024
        a = rhs(x)
        b = rhs(x + dx*a/2)
        c = rhs(x + dx*b/2)
        d = rhs(x + dx*c)
        x += dx * (a + 2b + 2c + d) / 6
    end
    @test plug ≈ x atol=2e-11

    tracer_variance = 2/30 - 2*(1 - exp(-30))/30^2
    @test axial_peclet(tracer_variance) ≈ 30
    coarse = bed_prediction(kinetics, 60, mass_kg, acid_mol_s; peclet=30, cells=2048)
    fine = bed_prediction(kinetics, 60, mass_kg, acid_mol_s; peclet=30, cells=4096)
    @test coarse.conversion ≈ fine.conversion atol=1e-6
    @test 0 < fine.conversion < plug
    @test fine.balance_error < 1e-10
    mixed = bed_prediction(kinetics, 60, mass_kg, acid_mol_s; peclet=0.01, cells=256)
    @test 0 < mixed.conversion < fine.conversion
    @test mixed.balance_error < 1e-10
    @test bed_prediction(kinetics, 60, 0, acid_mol_s; peclet=30).conversion == 0
    @test_throws ArgumentError bed_prediction(kinetics, 60, mass_kg, 0)
end
