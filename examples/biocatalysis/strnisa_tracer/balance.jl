# Reuse the existing multigrid graph projection. Only the axial
# throughput constraint survives with impermeable transverse plates.
include(joinpath(@__DIR__,"../../adsorption/disordered_refinement/balance.jl"))
using .RefinedFaceProjection: hierarchy, project
using .RefinedFaceProjection.PorousFaceTransport: graph_laplacian, divergence
using LinearAlgebra

function project_axial(left,right,axis,raw,coords,target_sum)
    n=size(coords,2)
    scale=maximum(abs,raw)
    scale>0 && target_sum>0 || throw(ArgumentError("Positive flow required"))
    lap=graph_laplacian(left,right,n)
    levels=hierarchy(lap,coords)
    axial=axis .== 1
    q=project(lap,levels,left,right,raw/scale,"rectangular raw")
    basis=project(lap,levels,left,right,Float64.(axial),"axial constraint")
    denom=sum(basis[axial])
    denom>0 || error("No axial periodic flow cycle")
    q .+= ((target_sum/scale-sum(q[axial]))/denom).*basis
    q .*= scale
    maximum(abs,divergence(left,right,q,n)) <= 1e-10maximum(abs,q) ||
        error("The balanced carrier fails local conservation")
    q
end

function main(args)
    length(args)==3 || error("Use faces_directory count target_axial_sum_m3_s")
    dir=args[1]
    n=parse(Int,args[2])
    target=parse(Float64,args[3])
    BLAS.set_num_threads(1)
    load(file,T)=reinterpret(T,read(joinpath(dir,file)))
    left=load("left.i32",Int32); right=load("right.i32",Int32)
    axis=load("axis.u8",UInt8)
    raw=load("unbalanced.f64",Float64)
    coords=reshape(load("coords.i32",Int32),3,n)
    q=project_axial(left,right,axis,raw,coords,target)
    write(joinpath(dir,"balanced.f64"),q)
end
if abspath(PROGRAM_FILE)==@__FILE__
    main(ARGS)
end
