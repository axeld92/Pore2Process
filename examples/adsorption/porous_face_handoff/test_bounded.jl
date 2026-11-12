using Test
include(joinpath(@__DIR__,"../../../julia/adsorption/BoundedFaceTransport.jl"))
using .BoundedFaceTransport

@testset "Limited straight-flow sine and uniform control" begin
    n, side, speed, d = 128, .0013891725543141453, .0128, 1e-9
    dx = side/n
    left, right = collect(1:n), [2:n;1]
    axis, flow = ones(Int,n), fill(speed*dx^2,n)
    x = (left .- .5)*dx
    k = 2π/side
    initial = .5 .+ .25*sin.(k*x)
    times = (.10737255542652152,.21476412830883945)
    a = evolve_bounded(left,right,axis,flow,initial,dx,d,times;time_factor=.5)
    for (t,c) in zip(times,a.snapshots)
        exact = .5 .+ .25*exp(-d*k^2*t)*sin.(k*(x .- speed*t))
        @test sqrt(sum(abs2,c-exact)/n)/.25 < .01
    end
    @test all(h->h[3]>=minimum(initial)-1e-12 && h[4]<=maximum(initial)+1e-12,a.history)
    @test all(h->abs(h[2]-.5)<1e-12,a.history)
    u = evolve_bounded(left,right,axis,flow,fill(.5,n),dx,d,times)
    @test all(h->h[3]==h[4]==.5,u.history)
end
