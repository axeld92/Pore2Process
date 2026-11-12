using Test
include(joinpath(@__DIR__,"../../../julia/adsorption/PorousFaceTransport.jl"))
using .PorousFaceTransport

@testset "Periodic straight-flow translation and diffusion" begin
    n, side, speed, diffusivity = 128, .0013891725543141453, .0128, 1e-9
    dx = side/n
    left = collect(1:n)
    right = [2:n;1]
    flow = fill(speed*dx^2,n)
    matrix = face_operator(left,right,flow,n,dx,diffusivity)
    wave = 2π/side
    x = (left .- .5)*dx
    initial = .5 .+ .25*sin.(wave*x)
    times = (.10737255542652152,.21476412830883945)
    result = evolve(matrix,initial,times,left,right,diffusivity/dx^2)
    independent = exponential_reference(matrix,initial,times)
    lambda = -4diffusivity*sin(wave*dx/2)^2/dx^2
    speed_discrete = speed*sin(wave*dx)/(wave*dx)
    for (t,actual,reference) in zip(times,result.snapshots,independent)
        continuous = .5 .+ .25*exp(-diffusivity*wave^2*t)*sin.(wave*(x .- speed*t))
        semidiscrete = .5 .+ .25*exp(lambda*t)*sin.(wave*(x .- speed_discrete*t))
        @test maximum(abs,actual-continuous)/.25 < .01
        @test maximum(abs,reference-semidiscrete) < 2e-13
        @test maximum(abs,actual-reference)/.25 < 1e-6
    end
    uniform = evolve(matrix,fill(.5,n),times,left,right,diffusivity/dx^2)
    @test maximum(abs,uniform.snapshots[end].-.5) < 1e-13
end
