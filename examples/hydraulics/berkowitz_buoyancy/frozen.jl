# Hold the density field fixed and relax the stock water independently.
include(joinpath(@__DIR__,"coupled.jl"))
using .SaltBed, SparseArrays, LinearAlgebra, TOML
function audit(parent,label)
    c=TOML.parsefile(joinpath(parent,"runtime.toml"))
    g=SaltBed.poregraph(c); h=c["h"]; nx=c["nx"]*c["repeats"]
    factor=cholesky(Symmetric(SaltBed.graph_laplacian(g.left,g.right,size(g.flow.low,2))))
    old=reshape(SaltBed.loadbin(joinpath(c["saved"],"field.f64"),Float64),7,:)
    archive=repeat(old[4:6,:]./c["native_speed_scale"],1,c["repeats"])
    salt=SaltBed.loadbin(joinpath(parent,"salt_"*label*".bin"),Float64)
    online=reshape(SaltBed.loadbin(joinpath(parent,"water_"*label*".bin"),Float64),3,:)
    w=salt[g.native_ids]; w./=sum(w)
    function transport_velocity(u)
        SaltBed.project!(g,factor,u.-archive,h,c["q"],nx)
        (g.flow.low.+g.flow.high).*(h/2)
    end
    live=transport_velocity(online)
    io=open(Cmd(readlines(joinpath(parent,"native_args.txt"))),"r+")
    log=open("native.log","w")
    try
        SaltBed.await_tag(io,"SALT_READY",log)
        SaltBed.water!(io,log,zero(salt),c["warm"],0.,1.)
        output=open("relaxation.csv","w")
        println(output,"sweeps,mean_z_m_s,live_mean_z_m_s,weighted_rms_live_error_m_s,weighted_rms_previous_change_m_s,native_Q_m3_s,max_Mach,max_abs_density_departure")
        previous=live; steps=0
        for step in (96,288,768,1536,3072,6144)
            u,raw=SaltBed.water!(io,log,salt,step-steps,c["gravity_accel"],1.)
            v=transport_velocity(u)
            rms(x)=sqrt(sum(w.*vec(sum(abs2,x;dims=1))))
            dz=dot(w,view(v,3,:))
            change=rms(v.-previous)
            println(output,join((step,dz,dot(w,view(live,3,:)),rms(v.-live),change,raw...),","))
            flush(output)
            step>=1536 && change<1e-7 && break
            previous=v; steps=step
        end
        close(output)
    finally
        close(log)
        if process_running(io); println(io,"-1 0 1"); flush(io); end
        close(io)
    end
end
audit(abspath(ARGS[1]),ARGS[2])
