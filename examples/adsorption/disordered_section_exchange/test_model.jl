using Test
include("model.jl")
using .SectionExchange: training, local_rules, propagate

@testset "Section-dependent exchange and tracer moments" begin
    rank = collect(1.:20.)
    train = hcat(rank,3rank)
    fit = training(train)
    # Keep the donor ranks at one seam, then reverse them at the next seam.
    # Delay is deterministic within each class, so the Markov closure is exact.
    donor = repeat(collect(1.:10.),inner=2)
    duration = hcat(donor,2donor,3reverse(donor),4reverse(donor),5donor,6donor,7donor)
    rules = local_rules(repeat(duration,2,1),fit.sizes)
    diag_counts(r) = sum(r.counts[i,i] for i in 1:10)
    @test diag_counts(rules[1]) == 40
    @test diag_counts(rules[2]) == 0
    predicted = propagate(fit.states,rules)
    arrival = train[:,2] .+ cumsum(duration[:,2:end];dims=2)
    avg(values) = sum(values)/length(values)
    @test predicted[:,1] ≈ ones(6)
    @test predicted[:,2] ≈ [avg(a) for a in eachcol(arrival)]
    @test predicted[:,3] ≈ [avg((a.-avg(a)).^2) for a in eachcol(arrival)]
    @test_throws ArgumentError training(hcat(train,rank))
    @test_throws ArgumentError local_rules(fill(NaN,40,7),fit.sizes)
end
