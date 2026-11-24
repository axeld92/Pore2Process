using Test
using Statistics
include("balance.jl")
include("track.jl")

@testset "Axial-only shared-flow projection" begin
    # One periodic ring, no transverse faces: all four rates must
    # converge to the prescribed rate despite nonuniform raw values.
    coords=Int32[0 1 2 3; 0 0 0 0; 0 0 0 0]
    left=Int32[1,2,3,4]; right=Int32[2,3,4,1]; axis=ones(UInt8,4)
    q=project_axial(left,right,axis,[1.1,0.85,0.9,1.0].*1e-12,coords,4e-12)
    @test q ≈ fill(1e-12,4) rtol=1e-12
    @test maximum(abs,divergence(left,right,q,4)) < 1e-23
end

function plug(n,dx,u)
    before=zeros(Int32,3,n);after=zeros(Int32,3,n)
    before[1,2:end]=Int32.(1:n-1);after[1,1:end-1]=Int32.(2:n)
    low=zeros(3,n);high=zeros(3,n)
    low[1,:].=u*dx^2;high[1,:].=u*dx^2
    before,after,low,high
end

@testset "Open tracer flow, zero and nonzero diffusion" begin
    dx=.001;u=.02;n=50;L=n*dx
    b,a,l,h=plug(n,dx,u)
    ts=passage(b,a,l,h,n÷2,dx,0.,.002,5.,32,15)
    @test all(isfinite,ts)
    @test maximum(abs.(ts .- L/u)) < 1e-10
    D=2e-7 # Pe=5000; closed-vessel and high-Pe inverse-Gaussian
             # leading analytical moments differ negligibly.
    ts=passage(b,a,l,h,n÷2,dx,D,.002,5.,4096,21)
    @test all(isfinite,ts)
    @test mean(ts) ≈ L/u rtol=0.01
    @test var(ts;corrected=false) ≈ 2D*L/u^3 rtol=0.15
    @test_throws ErrorException passage(b,a,l,h,n÷2,dx,D,.2,5.,32,15)
end
