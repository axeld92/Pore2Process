module SaltBed
using LinearAlgebra, SparseArrays, Random, TOML
include(joinpath(@__DIR__,"../berkowitz_homogeneous/first_passage.jl"))
include(joinpath(@__DIR__,"../../../julia/adsorption/PorousFaceTransport.jl"))
using .GlassPassage: Flow, advance, diffuse
using .PorousFaceTransport: graph_laplacian, divergence

loadbin(path,T) = collect(reinterpret(T,read(path)))

"""Tile the saved six-face pore graph, connecting seams to the adjacent copy."""
function poregraph(c)
    folder=c["saved"]; h=c["h"]; nx=c["nx"]; R=c["repeats"]
    rows=reshape(loadbin(joinpath(folder,"field.f64"),Float64),7,:)
    active=Int.(loadbin("active.i32",Int32)); n=length(active)
    ijk=Int.(rows[1:3,active])
    load(name,T)=loadbin(joinpath(folder,"faces",name),T)
    l0=load("left.i32",Int32); r0=load("right.i32",Int32); a0=load("axis.u8",UInt8)
    before=zeros(Int32,3,R*n); after=similar(before); fill!(after,0)
    left=repeat(l0,R); right=repeat(r0,R); axis=repeat(a0,R)
    locations=repeat(ijk,1,R); native_ids=repeat(active,R)
    for rep in 0:R-1
        slice=rep*n+1:(rep+1)*n
        locations[1,slice].+=rep*nx; native_ids[slice].+=rep*size(rows,2)
        for j in eachindex(l0)
            e=rep*length(l0)+j; a=a0[j]
            shift=a==1 && ijk[1,l0[j]]==nx-1
            l=rep*n+l0[j]; r=mod(rep+shift,R)*n+r0[j]
            left[e]=l; right[e]=r
            before[a,r]=l; after[a,l]=r
        end
    end
    base=repeat(load("balanced.f64",Float64),R)./h^3
    base .*= c["q"]/c["bulk_Q_before_projection_m3_s"]
    low=zeros(3,R*n); high=similar(low); fill!(high,0)
    flow=Flow(before,after,low,high)
    for e in eachindex(base)
        low[axis[e],right[e]]=high[axis[e],left[e]]=base[e]
    end
    flow.slope .= high.-low
    inlet=findall(p->locations[1,p]==0 && before[1,p]!=0,axes(low,2))
    faces=findall(==(1),axis)
    (;locations,native_ids,left,right,axis,base,flow,inlet,faces)
end

"""A reciprocal linear kernel spreads marker mass among liquid voxels only.
Each kernel node is normalized by its actual represented pore volume."""
function kernel(loc,nx,nt,h,width)
    nb=(max(1,round(Int,nx*h/width)),max(1,ceil(Int,nt*h/width)),max(1,ceil(Int,nt*h/width)))
    pitches=(nx/nb[1],nt/nb[2],nt/nb[3])
    ids=Matrix{Int32}(undef,8,size(loc,2)); weights=Matrix{Float64}(undef,size(ids))
    volume=zeros(prod(nb))
    for p in axes(loc,2)
        s=ntuple(a->(loc[a,p]+.5)/pitches[a]-.5,3)
        k=floor.(Int,s); f=s.-k
        for b in 0:7
            corner=ntuple(a->(b>>(a-1))&1,3)
            ix=mod(k[1]+corner[1],nb[1])
            iy=clamp(k[2]+corner[2],0,nb[2]-1); iz=clamp(k[3]+corner[3],0,nb[3]-1)
            i=1+ix+nb[1]*(iy+nb[2]*iz)
            w=prod(ntuple(a->corner[a]==1 ? f[a] : 1-f[a],3))
            ids[b+1,p]=i; weights[b+1,p]=w; volume[i]+=h^3*w
        end
    end
    (;ids,weights,volume,nb)
end

"""Redistribute only overfilled kernel nodes to the nearest available liquid
bins. Nonnegative mass, total dose and the upper feed concentration are retained."""
function limit_mass!(mass,volume,shape)
    moved=0.; grid=CartesianIndices(shape); linear=LinearIndices(shape)
    for origin in eachindex(mass)
        excess=mass[origin]-volume[origin]
        excess>0 || continue
        moved+=excess; mass[origin]=volume[origin]
        center=Tuple(grid[origin]); radius=0
        while excess>0
            radius+=1
            radius<=sum(shape) || error("Salt volume exceeds the packing's liquid capacity")
            candidates=Int[]
            for dz in -radius:radius, dy in -(radius-abs(dz)):(radius-abs(dz))
                dx=radius-abs(dz)-abs(dy)
                for sx in (dx==0 ? (0,) : (-dx,dx))
                    ix=mod1(center[1]+sx,shape[1]); iy=center[2]+dy; iz=center[3]+dz
                    1<=iy<=shape[2] && 1<=iz<=shape[3] || continue
                    j=linear[ix,iy,iz]
                    volume[j]>mass[j] && push!(candidates,j)
                end
            end
            unique!(candidates)
            room=sum((volume[j]-mass[j] for j in candidates);init=0.)
            room>0 || continue
            take=min(excess,room)
            for j in candidates
                mass[j]+=(volume[j]-mass[j])*(take/room)
            end
            excess=room>=excess ? 0. : excess-room
        end
    end
    moved
end

function deposit!(salt,k,native_ids,pores,alive,marker_volume,bounded=true)
    mass=zeros(length(k.volume)); fill!(salt,0)
    @inbounds for m in eachindex(pores)
        alive[m] || continue
        p=pores[m]
        for b in 1:8
            mass[k.ids[b,p]]+=marker_volume*k.weights[b,p]
        end
    end
    moved=bounded ? limit_mass!(mass,k.volume,k.nb) : 0.
    for i in eachindex(mass)
        k.volume[i]>0 && (mass[i]/=k.volume[i])
    end
    @inbounds for p in eachindex(native_ids)
        s=0.
        for b in 1:8
            s+=k.weights[b,p]*mass[k.ids[b,p]]
        end
        salt[native_ids[p]]=s
    end
    moved
end

"""Add the solenoidal part of the stock water-velocity change to the saved
balanced carrier. The pump contribution keeps the independent prescribed Q."""
function project!(g,factor,delta,h,qtarget,nx)
    q=similar(g.base)
    @inbounds for e in eachindex(q)
        a=g.axis[e]
        q[e]=(delta[a,g.native_ids[g.left[e]]]+delta[a,g.native_ids[g.right[e]]])/(2h)
    end
    pot=factor\divergence(g.left,g.right,q,size(g.flow.low,2))
    @inbounds for e in eachindex(q)
        q[e]+=pot[g.right[e]]-pot[g.left[e]]
    end
    excess=sum(e->q[e],g.faces)*h^3/(nx*qtarget)
    @. q=q+(1-excess)*g.base
    plane=zeros(nx)
    @inbounds for e in eachindex(q)
        a=g.axis[e]; l=g.left[e]; r=g.right[e]
        g.flow.low[a,r]=g.flow.high[a,l]=q[e]
        a==1 && (plane[1+g.locations[1,l]]+=q[e]*h^3)
    end
    @. g.flow.slope=g.flow.high-g.flow.low
    error=maximum(abs,divergence(g.left,g.right,q,size(g.flow.low,2)))/maximum(abs,q)
    error, maximum(abs,plane./qtarget.-1),excess
end

function await_tag(io,tag,log)
    while !eof(io)
        line=readline(io)
        println(log,line); flush(log)
        startswith(line,tag*" ") && return parse.(Float64,split(line)[2:end])
    end
    error("Native flow solver exited before $tag")
end
function water!(io,log,salt,steps,gravity,pump)
    write("salt.bin",salt)
    println(io,steps," ",gravity," ",pump); flush(io)
    diagnostics=await_tag(io,"SALT_DONE",log)
    u=reshape(loadbin("water.bin",Float64),3,:)
    all(isfinite,u) && all(isfinite,diagnostics) || error("Non-finite native water solution")
    u,diagnostics
end

function markers(c)
    n=c["markers"]
    times=fill(NaN,3,n)
    (;p=zeros(Int,n),xyz=fill((0.,0.,0.),n),ax=zeros(Int,n),next=ones(Int,n),
      entry=(collect(1:n).-0.5).*(c["feed_s"]/n),rng=[Xoshiro(c["seed"]+i) for i in 1:n],
      times,exit_z=fill(NaN,n),low_x=zeros(Int,n))
end

function move!(s,g,c,levels,start,stop)
    inlet=g.inlet
    # The distributor supplies salt proportional to forward water influx.
    cdf=cumsum([max(0.,g.flow.low[1,p]) for p in inlet])
    Threads.@threads :static for m in eachindex(s.p)
        entry=s.entry[m]
        (entry>=stop || s.next[m]>length(levels)) && continue
        t=max(start,entry); p=s.p[m]; xyz=s.xyz[m]; axial=s.ax[m]; next=s.next[m]
        rng=s.rng[m]; times=view(s.times,:,m)
        if p==0
            p=inlet[searchsortedfirst(cdf,rand(rng)*last(cdf))]
            xyz=(0.,rand(rng),rand(rng))
        end
        while t<stop-1e-10
            dt=min(c.dt,stop-t); half=.5*dt
            xyz,p,axial,next=advance(g.flow,xyz,p,axial,t,half,levels,times,next)
            if next<=length(levels)
                sigma=sqrt(2*c.D*dt)/c.h
                for a in (isodd(floor(Int,t/c.dt)) ? (1,2,3) : (3,2,1))
                    xyz,p,axial,next=diffuse(g.flow,xyz,p,axial,a,sigma*randn(rng),t+half,levels,times,next)
                    next>length(levels) && break
                end
            end
            if next<=length(levels)
                xyz,p,axial,next=advance(g.flow,xyz,p,axial,t+half,half,levels,times,next)
            end
            s.low_x[m]=min(s.low_x[m],axial)
            if next>length(levels)
                s.exit_z[m]=(g.locations[3,p]+xyz[3]-c.nt/2)*c.h
                break
            end
            t+=dt
        end
        s.p[m]=p; s.xyz[m]=xyz; s.ax[m]=axial; s.next[m]=next
    end
end

function run()
    c=TOML.parsefile("runtime.toml")
    # Specialize the time-stepping constants instead of indexing a heterogeneous
    # dictionary inside each marker substep.
    g=poregraph(c); h=c["h"]; nx=c["nx"]*c["repeats"]; n=c["markers"]
    cfg=(;dt=Float64(c["dt"]),D=Float64(c["diffusivity"]),h=Float64(h),nt=Int(c["nt"]))
    levels=round.(Int,[.05,.1,c["exit_m"]]./h)
    all(abs.(levels.*h.-[.05,.1,c["exit_m"]]).<1e-10) || error("Observation plane falls between faces")
    active=c["gravity_accel"]!=0
    factor=cholesky(Symmetric(graph_laplacian(g.left,g.right,size(g.flow.low,2))))
    k=active ? kernel(g.locations,nx,c["nt"],h,c["kernel_m"]) : nothing
    salt=zeros(c["native_rows"]*c["repeats"]); s=markers(c)
    io=open(Cmd(readlines("native_args.txt")),"r+")
    log=open("native.log","w"); stats=open("history.csv","w")
    println(stats,"time_s,inside,recovered,centroid_z_m,fraction_below_midplane,feed_fraction_max,salt_volume_error_m3,native_Q_m3_s,max_Mach,max_abs_density_departure,face_div_relative,plane_Q_relative,pump_correction_fraction,redistributed_fraction_of_dose,reverse_entry_Q_fraction")
    try
        await_tag(io,"SALT_READY",log)
        base,native=water!(io,log,salt,c["warm"],0.,1.)
        old=reshape(loadbin(joinpath(c["saved"],"field.f64"),Float64),7,:)
        archive=repeat(old[4:6,:]./c["native_speed_scale"],1,c["repeats"])
        println("Warm water velocity change relative to saved field ",norm(base.-archive)/norm(archive))
        length(base)==3length(salt) || error("Unexpected water array")
        maxdiv,maxplane,pumpcorr=project!(g,factor,base.-archive,h,c["q"],nx)
        masserr=0.; maxC=0.; moved=0.
        t=0.; interval=c["coupling_s"]
        while t<c["stop_s"]
            mid=t+interval/2
            move!(s,g,cfg,levels,t,mid)
            alive=(s.p.>0).&(s.next.<=length(levels))
            if active
                moved=deposit!(salt,k,g.native_ids,s.p,alive,c["source_m3"]/n,c["bounded"]==1)/c["source_m3"]
                masserr=sum(salt)*h^3-count(alive)*c["source_m3"]/n
                maxC=maximum(salt)
                c["bounded"]==0 || maxC<=1+1e-12 || error("Feed concentration exceeded")
                u,native=water!(io,log,salt,c["sweeps"],c["gravity_accel"],1.)
                if mid-interval/2 < 600 <= mid+interval/2 || mid-interval/2 < 3600 <= mid+interval/2
                    label=mid<1000 ? "early" : "hour"
                    write("salt_"*label*".bin",salt)
                    write("water_"*label*".bin",u)
                end
                u.-=archive
                maxdiv,maxplane,pumpcorr=project!(g,factor,u,h,c["q"],nx)
                maxdiv<1e-9 && maxplane<1e-9 || error("Water face balance failed")
                abs(masserr)<1e-10*c["source_m3"] || error("Density deposition lost salt")
            end
            inside=count(alive)
            if inside>0
                height=[(g.locations[3,s.p[m]]+s.xyz[m][3]-c["nt"]/2)*h for m in findall(alive)]
                z=sum(height)/inside; below=count(<(0),height)/inside
                reverse=sum(p->max(0.,-g.flow.low[1,p]),g.inlet)*h^3/c["q"]
                println(stats,join((mid,inside,count(>(length(levels)),s.next),z,below,maxC,masserr,native...,maxdiv,maxplane,pumpcorr,moved,reverse),","))
                flush(stats)
            end
            move!(s,g,cfg,levels,mid,t+interval)
            t+=interval
            all(>(length(levels)),s.next) && break
        end
        all(isfinite,s.times) || error("Incomplete first-passage recovery at observation cutoff")
        guard=c["repeats"]*c["nx"]-last(levels)
        minimum(s.low_x)>-guard || error("Upstream pulse reached the downstream image")
        write("arrivals.f64",s.times); write("outlet_z.f64",s.exit_z); write("entries.f64",s.entry)
        println("All ",n," salt markers recovered. Furthest upstream voxel ",minimum(s.low_x),"; guard ",guard)
    finally
        close(stats); close(log)
        if process_running(io)
            println(io,"-1 0 1"); flush(io)
        end
        close(io)
    end
end
if abspath(PROGRAM_FILE)==@__FILE__
    run()
end
end
