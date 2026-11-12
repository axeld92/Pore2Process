include(joinpath(@__DIR__, "../../../julia/adsorption/DiskFaceTransport.jl"))
using .DiskFaceTransport
using LinearAlgebra

function main(args)
    length(args) == 10 || error("Use output_dir nr na radius depth omega diffusivity dt time_mid time_final")
    directory = args[1]
    nr, na = parse.(Int, args[2:3])
    r, depth, omega, d, dt, t_mid, t_final = parse.(Float64, args[4:10])
    Base.ENDIAN_BOM == 0x04030201 || error("The saved-field contract requires a little-endian host")
    BLAS.set_num_threads(1)
    grid = sector_grid(nr, na; radius_m=r, depth_m=depth, omega_rad_s=omega, diffusivity_m2_s=d)
    initial = Matrix{Float64}(undef, nr, na)
    open(joinpath(directory, "initial.f64"), "r") do io
        read!(io, initial)
        eof(io) || error("Unexpected bytes after initial field")
    end
    result = transport(grid, initial, (t_mid,t_final); dt_max_s=dt)
    for (name, field) in zip(("middle","final"), result.fields)
        write(joinpath(directory, name*".f64"), field)
    end
    open(joinpath(directory, "history.csv"), "w") do io
        println(io, "time_s,liquid_inventory_kg,min_concentration_kg_m3,max_concentration_kg_m3")
        for row in result.history
            println(io, join(row, ','))
        end
    end
    open(joinpath(directory, "rings.csv"), "w") do io
        println(io, "ring,r_inner_m,r_outer_m,r_center_m,cell_volume_m3,theta_flow_m3_s,theta_diffusion_m3_s,outer_radial_diffusion_m3_s,outer_radial_flow_m3_s")
        for i in 1:nr
            println(io, join((i,grid.edges[i],grid.edges[i+1],grid.centers[i],grid.volume[i],
                              grid.q[i],grid.g_angular[i],grid.g_radial[i+1],0.0), ','))
        end
    end
end

main(ARGS)
