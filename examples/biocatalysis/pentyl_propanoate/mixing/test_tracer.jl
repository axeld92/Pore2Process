using Test

include("tracer.jl")

@testset "Periodic tracer with known transit times" begin
    width = 8
    before = zeros(Int32, 3, width)
    after = similar(before)
    fill!(after, 0)
    before[1, :] .= [width; collect(1:width - 1)]
    after[1, :] .= [collect(2:width); 1]
    speed = zeros(3, width)
    speed[1, :] .= 2.
    flow = Flow(before, after, speed, speed)
    # The synthetic field wraps twice, then exits within an axial voxel.
    dx, length = 0.01, 0.185
    result = sample(flow, collect(0:width - 1), width, 3, dx, length, 0.,
                    0.10, 20., 32, 426620)
    @test all(isapprox.(result[1, :], length / (2dx); atol=2e-12, rtol=0.))
    @test all(iszero, result[2, :])
    @test axial_step((0.9, 0., 0.), 8, (0.1, 0., 0.), 1, collect(0:width - 1), width) ≈ 0.2
    @test axial_step((0.1, 0., 0.), 1, (0.9, 0., 0.), 8, collect(0:width - 1), width) ≈ -0.2
end
