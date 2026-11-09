@testset "Measured hydraulic translation and independent capillary oracle" begin
    r,L,mu,A,dp = 50e-6,4e-4,1e-3,4e-8,0.1
    ref=capillary_reference(r,L,mu,dp,A)
    @test ref.flow_m3_s ≈ 6.135923151542565e-13 rtol=1e-14
    @test ref.permeability_cell_m2 ≈ 6.135923151542565e-11 rtol=1e-14
    @test ref.permeability_pore_m2 ≈ 3.125e-10 rtol=1e-14
    @test ref.porosity ≈ pi/16
    @test infer_permeability(L,A,mu,dp,ref.flow_m3_s) ≈ ref.permeability_cell_m2
    @test infer_permeability(L,A,mu,-dp,-ref.flow_m3_s) ≈ ref.permeability_cell_m2
    @test infer_permeability(L,2A,mu,dp,ref.flow_m3_s) ≈ ref.permeability_cell_m2/2
    @test infer_permeability(L,A,mu,dp,0.7ref.flow_m3_s) ≈ 0.7ref.permeability_cell_m2
    @test isnothing(infer_permeability(L,A,mu,0.,0.))
    @test predict_column_flow(ref.permeability_cell_m2,mu,L,A,dp) ≈ ref.flow_m3_s
    @test predict_column_pressure_drop(ref.permeability_cell_m2,mu,L,A,ref.flow_m3_s) ≈ dp
    @test predict_column_pressure_drop(ref.permeability_cell_m2,mu,2L,A,ref.flow_m3_s) ≈ 2dp
    @test predict_column_pressure_drop(ref.permeability_cell_m2,mu,L,A,-ref.flow_m3_s) ≈ -dp
    @test predict_column_pressure_drop(ref.permeability_cell_m2,mu,L,A,0.) == 0.
    @test isnothing(predict_column_pressure_drop(nothing,mu,L,A,1.))
    @test_throws ArgumentError predict_column_pressure_drop(0.,mu,L,A,1.)
    @test_throws ArgumentError predict_column_pressure_drop(ref.permeability_cell_m2,mu,L,A,NaN)
    @test capillary_reference(r,L,mu,-dp,A).flow_m3_s ≈ -ref.flow_m3_s
    @test_throws ArgumentError infer_permeability(L,A,mu,0.,1.)
    @test_throws ArgumentError infer_permeability(L,A,mu,dp,-ref.flow_m3_s)
    @test_throws ArgumentError infer_permeability(0.,A,mu,dp,1.)
    @test_throws ArgumentError infer_permeability(L,A,mu,NaN,1.)
    @test_throws ArgumentError infer_permeability(L,A,mu,dp,Inf)
    @test_throws ArgumentError capillary_reference(r,L,mu,dp,pi*r^2/2)
end
