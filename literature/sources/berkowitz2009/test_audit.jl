using Test
using TOML
include("audit.jl")
const A = BerkowitzAudit

@testset "Figure axes, interval widths and partial-window moments" begin
    scale = A.axis_map([20.0, 10.0, 0.0], [0, 5, 10])
    @test A.physical(scale, 16.0) ≈ 2.0
    @test scale.max_tick_error < 1e-12
    t = [0.0, 1.0, 2.0, 3.0, 4.0]
    ec = [7.0, 8.0, 9.0, 8.0, 7.0]
    @test A.half_height(t, ec, 7.0) == (left=1.0, right=3.0, width=2.0)
    @test A.half_height(t .+ 110, 3ec .+ 2, 23.0).width ≈ 2.0
    @test A.visible_mean(t, ec, 7.0) ≈ 2.0
    @test A.visible_mean([0.0, 2.0], [0.0, 1.0], 0.0) ≈ 4/3
    @test A.interpolate(t, ec, 2.5) ≈ 8.5
    @test_throws DomainError A.interpolate(t, ec, 8.0)
    @test A.convolve([0.25, 0.75], [0.5, 0.5]) == [0.125, 0.5, 0.375]
end

result = A.evaluate()
source = TOML.parsefile(joinpath(@__DIR__, "evidence.toml"))
saved = TOML.parsefile(joinpath(@__DIR__, "reference.toml"))
@testset "Published inventory and vector readings" begin
    @test result["code_and_input_sha256"] == saved["code_and_input_sha256"]
    @test result["pore_volume_ml"] ≈ 0.39π * 1.35^2 * 40
    @test result["injected_NaCl_mg"] ≈ 0.75
    @test result["column_to_bead_diameter_ratios"] == [27.0, 6.75]
    @test source["pulse"]["flow_ml_min"] .* source["pulse"]["pulse_min"] ≈ fill(1.5, 3)
    @test result["time_axis_tick_error_min"] < 0.12
    @test result["ec_axis_tick_error_uS"] < 0.001
    curves = result["curves"]
    @test length(curves) == 10
    @test sum(s["markers"] for s in curves) == 700
    @test sum(s["sample_times"] for s in curves) == 589
    @test all(s -> (s["times_with_multiple_markers"] > 0) == (s["flow_ml_min"] == 0.4), curves)
    @test all(s -> s["fwhm_reading_bounds_min"][1] <= s["fwhm_min"] <= s["fwhm_reading_bounds_min"][2], curves)
    for (m, archived) in zip(curves, saved["curves"])
        @test m["id"] == archived["id"] && m["fwhm_min"] ≈ archived["fwhm_min"]
    end
end

@testset "Independent series reference" begin
    levels = result["baseline_convergence"]
    for r in levels
        @test abs(r["mass"] - 1) < 1e-10
        @test abs(r["mean_pv"] - r["analytical_mean_pv"]) < 1e-10
        @test abs(r["variance_pv2"] - r["analytical_variance_pv2"]) < 1e-6
        @test r["max_directional_density_difference_per_pv"] < 1e-10
    end
    @test abs(levels[end]["fwhm_pv"] / levels[end-1]["fwhm_pv"] - 1) < 1e-5
    @test abs(levels[end]["peak_pv"] / levels[end-1]["peak_pv"] - 1) < 1e-5
    @test levels[end]["collection_center_fwhm_pv"] > levels[end]["fwhm_pv"]
    @test all(r -> r["pulse_centroid_min"] ≈ (result["pore_volume_ml"] + 0.75)/r["flow_ml_min"], result["rates"])
end

@testset "Published-flow comparison" begin
    slow, middle, fast = result["rates"]
    @test slow["dispersivity_basis"] == "reported single-packing fit"
    @test fast["dispersivity_basis"] == "fixed-alpha transfer scenario"
    @test slow["CF_minus_FC_fwhm_bounds_min"][1] > 25
    @test middle["CF_minus_FC_fwhm_bounds_min"][1] > 13
    @test middle["CF_minus_FC_peak_min"][1] < 0 < middle["CF_minus_FC_peak_min"][2]
    @test 157 < fast["FC_peak_min"][1] < 158
    @test 160 < fast["CF_peak_min"][1] < 164
    @test 84 < fast["baseline_peak_min"] < 86
    @test all(r -> 39 < r["FC_fwhm_pumped_ml"][1] <= r["FC_fwhm_pumped_ml"][2] < 43, result["rates"])
end
