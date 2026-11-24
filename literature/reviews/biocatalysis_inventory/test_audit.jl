using Test
include("audit.jl")
using .InventoryAudit

const A = InventoryAudit

@testset "Axes and geometry units" begin
    @test A.axis_value(70,110,10,0.,100.) == 40
    @test A.axis_value(40,10,110,0.,50.) == 15
    @test A.interp([0.,2.,4.],[0.,10.,14.],3.) == 12
    @test_throws DomainError A.interp([1.,2.],[3.,4.],0.5)
    @test_throws ArgumentError A.interp([1.,1.],[3.,4.],1.)
    # A 10 mL dry layer expands to 20 mL in a 20 mL column.
    @test A.layer_fill(4.,20.,0.4,2.) == (0.5,1.)
    @test_throws ArgumentError A.layer_fill(4.,20.,0.,2.)
end

@testset "Independent analytical activity curves" begin
    t=[0.,1.,2.,4.,10.]
    k=0.41
    a=exp.(-k .* t)
    @test A.early_log_coefficient(t,a) ≈ k rtol=1e-13
    fit=A.fit_activity(t,a,false)
    @test fit.k ≈ k rtol=1e-7
    @test fit.retained == 0
    @test A.activity(fit.k,fit.retained,30.) ≈ exp(-k*30) rtol=1e-6
    f=0.23
    a=f .+ (1-f).*exp.(-k.*t)
    fit=A.fit_activity(t,a,true)
    @test fit.k ≈ k rtol=1e-7
    @test fit.retained ≈ f rtol=1e-7
    @test A.activity(fit.k,fit.retained,30.) ≈ f+(1-f)*exp(-k*30) rtol=1e-7
    @test A.fit_activity(t,ones(length(t)),false).k == 0
    @test_throws ArgumentError A.fit_activity(t,[-1.,1.,1.,1.,1.],true)
end

@testset "Source evidence and independent comparisons" begin
    result=A.run_audit()
    @test result["target_packed_bed_conversion_points_used"] == 0
    wet=result["damstrup"]["wetting"]
    @test wet[1]["free_wet_layer_fraction"] < 0.52
    @test isapprox(wet[2]["free_wet_layer_fraction"],1.;atol=0.02)
    @test wet[3]["free_wet_layer_fraction"] > 1.49
    @test 0.94 < result["damstrup"]["target_analogue_scenario"]["free_wet_layer_fraction"] < 0.97
    hold=result["damstrup"]["figure2_holdouts"]
    @test [q["n"] for q in hold] == [5,5]
    @test all(q->q["rmse_pp"]>10,hold)
    @test all(q->q["outside_reading_envelope_count"]>=3,hold)
    for s in result["rejasse"]["storage_series"]
        @test abs(s["author_half_life_h"]-s["half_life_using_ln10_in_h"]) < 0.001
        @test s["early_log_k_natural_per_min"] ≈ s["author_k_natural_per_min"] rtol=0.05
        @test s["early_log_prediction_10_30_percent"][2] < s["author_30min_remaining_percent"]
        @test all(f->f["grid_refinement_holdout_change_pp"]<1e-5,s["shape_fits"])
    end
    storage=result["rejasse"]["storage_series"]
    retained(s)=only(filter(x->x["shape"]=="fast_plus_retained",s["shape_fits"]))
    @test abs(retained(storage[1])["heldout_error_vs_printed_pp"]) < 2
    @test abs(retained(storage[2])["heldout_error_vs_printed_pp"]) > 5
    lerin=result["lerin"]
    @test length.(getindex.(lerin,"cycle")) == [10,8,10]
    @test lerin[1]["activity_percent"][1] < lerin[1]["activity_percent"][3]
    @test maximum(lerin[3]["activity_percent"]) > 130
    @test lerin[3]["activity_change_from_reference_pp"] < -45
    @test abs(lerin[3]["conversion_change_from_reference_pp"]) < 2
    s=result["strnisa"]
    @test length(s["same_contact_time_holdouts"]) == 5
    @test s["nominal_A_channel_clearance_um"] ≈ 27.9
    for r in s["inventory_and_productivity"]
        @test abs(r["calculated_gamma_U_ul_void"]-r["printed_gamma_U_ul_void"]) < 0.02
        @test !haskey(first(r["cases"]),"BB_production_from_BPN_mmol_h")
        @test haskey(last(r["cases"]),"BB_production_from_BPN_mmol_h")
    end
    a=first(s["inventory_and_productivity"])
    f=last(s["inventory_and_productivity"])
    @test a["inventory_reconstructed_mg"] ≈ 7.104
    @test f["inventory_reconstructed_mg"] ≈ 148.4625
    # Manual unit conversion: 13.56 mmol/g_CaLB * 0.0007104 g / 1.4 min * 60.
    @test last(a["cases"])["BB_production_from_BPN_mmol_h"] ≈ 0.4128438857142857
    @test all(r->last(r["cases"])["BB_to_stoichiometric_BB_ratio"] < 1,
              s["inventory_and_productivity"])
    @test [r["n"] for r in s["same_contact_time_holdouts"]] == [8,10,10,10,9]
end
