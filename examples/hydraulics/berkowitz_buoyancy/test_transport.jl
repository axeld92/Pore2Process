using Test, LinearAlgebra, SparseArrays
include("coupled.jl")
using .SaltBed
@testset "Conservative positive salt density" begin
    h=0.001
    points=[(x,y,z) for x in 0:7 for y in 0:3 for z in 0:3 if !(y==z==1)]
    loc=reduce(hcat,collect.(points))
    n=size(loc,2); k=SaltBed.kernel(loc,8,4,h,.002)
    @test sum(k.volume) ≈ n*h^3
    @test maximum(abs,vec(sum(k.weights;dims=1)).-1)<1e-14
    salt=zeros(n); ids=collect(1:n)
    SaltBed.deposit!(salt,k,ids,ids,trues(n),0.45*h^3)
    @test all(x->isapprox(x,.45,atol=1e-14),salt)
    SaltBed.deposit!(salt,k,ids,fill(1,n),trues(n),.26*h^3)
    @test minimum(salt)>=0 && maximum(salt)<=1+1e-14
    @test sum(salt) ≈ n*.26 rtol=1e-13
    # Moving through the periodic seam or coming up against a closed wall
    # must conserve a compact dose with local demand above capacity.
    mass=[1.6,0.,0.,0.]
    moved=SaltBed.limit_mass!(mass,ones(4),(4,1,1))
    @test moved≈.6
    @test mass≈[1.,.3,0.,.3]
end

@testset "Hydrostatic potential and fixed discharge" begin
    # Periodic axial channel with two layers and no vertical through-flow.
    nx=4; n=2nx; locations=zeros(Int,3,n)
    for j in 0:1, i in 0:nx-1
        locations[:,j*nx+i+1].=(i,0,j)
    end
    l=Int[]; r=Int[]; a=Int[]
    for p in 1:n
        push!(l,p); push!(r,(p-1)÷nx*nx+mod1(p+1,nx)); push!(a,1)
        if p<=nx; push!(l,p); push!(r,p+nx); push!(a,3); end
    end
    base=Float64.(a.==1); faces=findall(==(1),a)
    flow=SaltBed.Flow(zeros(Int32,3,n),zeros(Int32,3,n),zeros(3,n),zeros(3,n))
    graph=(;left=l,right=r,axis=a,base,locations,flow,faces,native_ids=1:n)
    factor=cholesky(Symmetric(SaltBed.graph_laplacian(l,r,n)))
    delta=repeat([.17,0.,-.8],1,n)
    div,plane,excess=SaltBed.project!(graph,factor,delta,1.,2.,nx)
    @test max(div,plane)<1e-13
    @test excess≈.17
    @test all(x->isapprox(x,1.,atol=1e-13),flow.high[1,:])
    @test maximum(abs,flow.high[3,:])<1e-13
end
