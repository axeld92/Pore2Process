using Test
include("solve.jl")

@testset "Open bounded transport" begin
    # Two well-mixed cells in series, no molecular flux.
    dx = 1.
    dt = .01
    q = 2.
    grid = Parallel.prepare(Int32[1],Int32[2],UInt8[1],[q],2,dx,0.)
    boundary = [q,-q]
    uniform = OpenFlow.evolve_open(grid,ones(2),boundary,dx,dt,10,10)
    @test uniform.concentration==ones(2)
    @test maximum(abs(h[2]+h[6]-h[5]-2) for h in uniform.history)<1e-14
    pulse = OpenFlow.evolve_open(grid,zeros(2),boundary,dx,dt,1,1)
    # After the first Euler stage c¹=(qdt,0). The centered face flux
    # passes the limiter on the second stage, giving c²₂=(qdt)^2/8
    # after the SSP blend. The outlet flux in the first two stages is zero.
    expected = dt*(4*q*(q*dt)^2/8)/6
    @test isapprox(last(pulse.history)[6],expected;rtol=1e-14,atol=0.)
    @test isapprox(last(pulse.history)[8],expected*(dt/2)^2;rtol=1e-14,atol=0.)
    @test isapprox(last(pulse.history)[2]+last(pulse.history)[6],q*dt;rtol=1e-14)
end
