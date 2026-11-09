"""Independent finite-difference reference for fully developed duct flow.

This is a verification oracle, never an OpenLB initializer or the process
hydraulic input. Solve -laplacian(w)=1 on a polygon with zero wall values.
Shortley-Weller cut distances locate the declared polygon walls. The returned
integral C has units m^4, and Q=C*G/mu for an axial pressure gradient G.
"""
function duct_stokes_reference(vertices::AbstractMatrix, resolution::Integer)
    size(vertices,2)==2 && size(vertices,1)>=3 || throw(ArgumentError("Require a polygon (N,2)."))
    16<=resolution<=512 || throw(ArgumentError("Reference resolution must be 16..512."))
    v=Float64.(vertices)
    all(isfinite,v) || throw(ArgumentError("Polygon must be finite."))
    span=maximum(v[:,1])-minimum(v[:,1]); span>0 || throw(ArgumentError("Degenerate polygon."))
    # Work in dimensionless coordinates to condition the sparse operator.
    v=v/span; h=1/resolution
    xmin,xmax=extrema(v[:,1]); ymin,ymax=extrema(v[:,2])
    xs=collect(xmin+h/2:h:xmax); ys=collect(ymin+h/2:h:ymax)
    function inside(x,y)
        flag=false
        for i in axes(v,1)
            j=mod1(i+1,size(v,1)); a,b=v[i,:],v[j,:]
            if (a[2]>y)!=(b[2]>y) && x<(b[1]-a[1])*(y-a[2])/(b[2]-a[2])+a[1]
                flag=!flag
            end
        end
        flag
    end
    ids=zeros(Int,length(xs),length(ys)); count=0
    for j in eachindex(ys), i in eachindex(xs)
        if inside(xs[i],ys[j]); count+=1; ids[i,j]=count; end
    end
    count>0 || throw(ArgumentError("No reference interior points."))
    function wall_distance(x,y,axis,sgn)
        distances=Float64[]
        for k in axes(v,1)
            a,b=v[k,:],v[mod1(k+1,size(v,1)),:]
            other=3-axis; query=other==1 ? x : y; origin=axis==1 ? x : y
            if min(a[other],b[other])<=query<=max(a[other],b[other]) && a[other]!=b[other]
                t=(query-a[other])/(b[other]-a[other]); d=sgn*(a[axis]+t*(b[axis]-a[axis])-origin)
                d>1e-12 && push!(distances,d)
            end
        end
        isempty(distances) && throw(ArgumentError("Unclosed polygon ray."))
        minimum(distances)
    end
    rows=Int[]; cols=Int[]; values=Float64[]
    for j in eachindex(ys), i in eachindex(xs)
        id=ids[i,j]; id==0 && continue
        diagonal=0.
        for axis in 1:2
            neighbors=Int[]; distances=Float64[]
            for sign in (-1,1)
                ni,nj=axis==1 ? (i+sign,j) : (i,j+sign)
                nid=checkbounds(Bool,ids,ni,nj) ? ids[ni,nj] : 0
                push!(neighbors,nid)
                push!(distances,nid>0 ? h : wall_distance(xs[i],ys[j],axis,sign))
            end
            a,b=distances
            diagonal+=2/(a*b)
            for (nid,d) in zip(neighbors,distances)
                if nid>0; push!(rows,id); push!(cols,nid); push!(values,-2/(d*(a+b))); end
            end
        end
        push!(rows,id);push!(cols,id);push!(values,diagonal)
    end
    matrix=sparse(rows,cols,values,count,count); w=matrix\ones(count)
    residual=norm(matrix*w.-1,Inf)
    (; conductance_geometry_m4=sum(w)*h^2*span^4, residual, nodes=count,
       area_quadrature_m2=count*h^2*span^2)
end

"""Axially distributed, leaky laminar channel with constant membrane impedance.

Input C is measured from OpenLB: mu_native*Q_native/G_native (m^4).
All lengths, flows, pressures, viscosity, areas and impedances use SI units.
No entrance, turbulent, osmotic or physical fouling-time model is implied.
The exact hyperbolic solution obeys p'=-mu*Q/C and Q'=-P*p/Z.
"""
function star_membrane_load(C_m4, area_channel_m2, perimeter_m, length_m,
        viscosity_Pa_s,density_kg_m3,Re_in,flux_m_s,membrane_area_m2,Z_support,Z_fixed)
    C=positive("native geometry coefficient_m4",C_m4)
    A=positive("channel open area_m2",area_channel_m2); P=positive("wetted perimeter_m",perimeter_m)
    L=positive("active length_m",length_m); mu=positive("viscosity_Pa_s",viscosity_Pa_s)
    rho=positive("density_kg_m3",density_kg_m3); re=positive("inlet Reynolds number",Re_in)
    J=positive("flux_m_s",flux_m_s); Am=positive("membrane area_m2",membrane_area_m2)
    zs=positive("support impedance_Pa_s_m",Z_support); zf=positive("fixed impedance_Pa_s_m",Z_fixed)
    Z=zs+zf; Dh=4A/P; u=re*mu/(rho*Dh); qi=u*A; qp=J*P*L
    qi>qp || throw(ArgumentError("Inlet channel flow must exceed permeate extraction."))
    resistance=mu/C; beta=sqrt(resistance*P/Z); t=beta*L
    t<20 || throw(ArgumentError("Axial model outside supported hyperbolic range."))
    pin=(qi*2sinh(t/2)^2+qp)*beta*Z/(P*sinh(t))
    pout=pin*cosh(t)-resistance*qi/beta*sinh(t)
    pout>0 || throw(ArgumentError("Negative outlet TMP: lower inlet speed, shorten channel, or increase membrane impedance/flux."))
    n=101; x=collect(range(0.,L,length=n))
    pressure=pin*cosh.(beta*x).-resistance*qi/beta*sinh.(beta*x)
    flow=qi*cosh.(beta*x).-P/(Z*beta)*pin*sinh.(beta*x)
    parallel=Am/(P*L); Qp=J*Am; Qm=qi*parallel; drop=pin-pout
    # High-circulation DWSIM lumping dissipates Qin*delta_p. Exact distributed
    # channel dissipation is integral(mu/C*Q^2 dx); trapezoidal error is O(h^2).
    friction=parallel*resistance*sum((flow[1:end-1].^2+flow[2:end].^2)/2)*(L/(n-1))
    membrane_loss=parallel*P/Z*sum((pressure[1:end-1].^2+pressure[2:end].^2)/2)*(L/(n-1))
    lumped=Qm*drop
    (; geometry_coefficient_m4=C, hydraulic_diameter_m=Dh,inlet_velocity_m_s=u,
       inlet_Re=re,outlet_Re=re*(qi-qp)/qi,viscosity_Pa_s=mu,density_kg_m3=rho,
       total_membrane_impedance_Pa_s_m=Z,support_impedance_Pa_s_m=zs,
       fixed_impedance_Pa_s_m=zf,selective_and_other_fixed_fraction=zf/Z,
       mean_TMP_Pa=Z*J,inlet_TMP_Pa=pin,outlet_TMP_Pa=pout,
       channel_drop_Pa=drop,axial_loss_over_mean_TMP=drop/(Z*J),
       channel_resistance_Pa_s_m3=drop/Qm,module_flow_m3_s=Qm,
       permeate_flow_m3_s=Qp,module_to_permeate_ratio=Qm/Qp,
       equivalent_parallel_channels=parallel,channel_inlet_flow_m3_s=qi,
       channel_permeate_flow_m3_s=qp,channel_outlet_flow_m3_s=qi-qp,
       slender_channel_inertia_screen=re*Dh/L*(qp/qi),
       wall_suction_velocity_ratio=J/u,
       distributed_channel_dissipation_W=friction,lumped_channel_dissipation_W=lumped,
       distributed_membrane_dissipation_W=membrane_loss,
       distributed_port_energy_residual_W=friction+membrane_loss-(Qm*drop+Qp*pout),
       lumped_channel_dissipation_relative_excess=lumped/friction-1,
       x_m=x,pressure_Pa=pressure,local_flux_m_s=pressure/Z,channel_flow_m3_s=flow)
end

"""Conservative dilute COD tracer on the distributed channel/recycle solution.

Constant local rejection is an explicitly transferred scenario parameter.
The supplied concentration is the module inlet, as in the source feed table;
fresh concentration is calculated, not silently held at the same value.
"""
function star_feed_balance(module_concentration_kg_m3,rejection,load,recovery)
    cm=positive("module feed concentration_kg_m3",module_concentration_kg_m3)
    r=finite("local rejection",rejection); 0<=r<=1 || throw(ArgumentError("Rejection must be in [0,1]."))
    y=positive("recovery",recovery); y<1 || throw(ArgumentError("Recovery must be below one."))
    qi,qo,qp=load.channel_inlet_flow_m3_s,load.channel_outlet_flow_m3_s,load.channel_permeate_flow_m3_s
    cr=cm*(qi/qo)^r; cp=(qi*cm-qo*cr)/qp
    Qm,Qp=load.module_flow_m3_s,load.permeate_flow_m3_s
    Qf=Qp/y; Qreturn=Qm-Qf; Qpurge=Qf-Qp
    Qreturn>0 || throw(ArgumentError("Positive recycle return required."))
    cf=(Qm*cm-Qreturn*cr)/Qf
    cf>=0 || throw(ArgumentError("Invalid fresh-feed concentration."))
    (; module_inlet_COD_kg_m3=cm,assumed_local_rejection=r,
       fresh_COD_kg_m3=cf,permeate_COD_kg_m3=cp,retentate_COD_kg_m3=cr,
       apparent_module_rejection=1-cp/cm,
       overall_COD_removal_fraction=1-y*cp/cf,
       module_residual_kg_s=Qm*cm-(Qm-Qp)*cr-Qp*cp,
       mixer_residual_kg_s=Qf*cf+Qreturn*cr-Qm*cm,
       overall_residual_kg_s=Qf*cf-Qp*cp-Qpurge*cr,
       axial_COD_kg_m3=cm*(qi./load.channel_flow_m3_s).^r)
end
