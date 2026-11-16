using Test
include("model.jl")
using .TravelMemory: predict

@testset "Within-class age retains correlated delay" begin
    # Same local path class for both cohorts, with distinct inlet ages.
    # A tracer from the older cohort acquires a longer delay on every section.
    local_path = 10. .* repeat(1:10;inner=8)
    inlet_age = repeat([1.,1.,1.,1.,3.,3.,3.,3.],10)
    train = hcat(inlet_age,inlet_age+local_path)
    delay = local_path + 3inlet_age
    probe = hcat(train,train[:,2] .+ delay .* transpose(1:6))
    expected = probe[:,3:8]
    average(x) = sum(x)/length(x)
    exact_variance = [average((v.-average(v)).^2) for v in eachcol(expected)]
    for groups in (2,4,8)
        p = predict(train,repeat(probe,2,1),groups).prediction
        @test p[:,1] ≈ ones(6)
        @test p[:,2] ≈ [average(v) for v in eachcol(expected)]
        @test p[:,3] ≈ exact_variance
    end
    collapsed = predict(train,probe,1)
    @test collapsed.prediction ≈ collapsed.baseline_prediction
    @test all(collapsed.prediction[:,3] .< exact_variance)
    @test_throws ArgumentError predict(train,probe[:,1:7],4)
end
