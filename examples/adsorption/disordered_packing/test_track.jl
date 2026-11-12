using Test
include("track.jl")

@testset "Eight unequal fluid-voxel counts" begin
    # One through-flow voxel per cell and successively more closed,
    # stationary padding voxels, to test the bookkeeping for each seam.
    counts = cumsum(1:8)
    owners = [1;counts[1:end-1].+1]
    n = last(counts)
    before = zeros(Int32,3,n)
    after = zeros(Int32,3,n)
    low = zeros(3,n)
    high = zeros(3,n)
    before[1,owners[2:end]] .= owners[1:end-1]
    after[1,owners[1:end-1]] .= owners[2:end]
    low[1,owners] .= 2.
    high[1,owners] .= 2.
    mktempdir() do d
        for (name,array) in (("before.i32",before),("after.i32",after),
                             ("incoming.f64",low),("outgoing.f64",high))
            write(joinpath(d,name),array)
        end
        result = sample(d,1.,0.,.01,8.,10,42,1,counts)
        @test all(result[1:8,:] .≈ collect(.5:.5:4.))
        @test all(isnan,result[9:16,:])
        old = TrackRepeats.sample(Flow(Int32[0 1;0 0;0 0],Int32[2 0;0 0;0 0],
                                  [2. 2.;0. 0.;0. 0.],[2. 2.;0. 0.;0. 0.]),
                                  1.,.04,.01,8.,1000,42,8,1)
        current = sample(d,1.,.04,.01,8.,1000,42,1,counts)
        @test isequal(current,old)
    end
end
