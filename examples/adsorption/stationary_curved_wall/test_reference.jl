using Test
include("radial_reference.jl")
using .RadialReference: solve,initial,ALPHA

const R=.0004
const D=1e-9
const T=(.10737255542652152,.21476412830883945)

@testset "Independent radial continuum reference" begin
    @testset "Pure diffusion and rigid rotation recover exact modes" begin
        for omega in (0.0,32.34323795541595)
            s,f=solve(256,4096,T; diffusivity=D,radius=R,omega,shear=false)
            for (k,t) in enumerate(T)
                exact=initial.(s).*exp((-D*(ALPHA/R)^2-2im*omega)*t)
                @test maximum(abs.(f[:,k]-exact))<1e-5
            end
        end
    end
    @testset "Diffusion-free shear follows characteristics" begin
        omega=84.0
        s,f=solve(128,8192,T; diffusivity=0.0,radius=R,omega)
        for (k,t) in enumerate(T)
            exact=initial.(s).*exp.(-2im*omega*t.*(1 .- s.^2))
            @test maximum(abs.(f[:,k]-exact))<1e-5
        end
    end
    @testset "Radial diffusion and shear dissipate concentration variance" begin
        s,f=solve(256,2048,T; diffusivity=D,radius=R,omega=84.0)
        energy(a)=sum(s .* abs2.(a))
        @test energy(f[:,2]) < energy(f[:,1]) < energy(initial.(s))
        @test maximum(abs,f)<1
    end
end
