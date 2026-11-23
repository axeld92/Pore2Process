module BeadMoisture

using LinearAlgebra: Tridiagonal

export relax_water

"""
Transient water storage in a spherical bead with reversible local sorption.

c in kg/m³ of pore liquid; storage in kg/m³ of bead. Both free and
sorbed water are mobile in the balance only through the pore-liquid
flux D*eps/tau times the concentration gradient. qmax is kg/kg catalyst.
The film coefficient (m/s) may be Inf for a fixed surface concentration.
The shell output describes water storage, not an enzyme activity.
"""
function relax_water(radius, eps, tau, rho_cat, d_bulk, kfilm,
                     qmax, K, c_start, c_bath;
                     cells=64, dt=2.0, duration=10800.0, sample_step=60.0,
                     shell_depth=80e-6)
    all(isfinite, (radius, eps, tau, rho_cat, d_bulk, qmax, K, c_start, c_bath,
                   dt, duration, sample_step, shell_depth)) &&
        radius > 0 && 0 < eps <= 1 && tau >= 1 && rho_cat > 0 &&
        d_bulk > 0 && kfilm > 0 && qmax >= 0 && K > 0 &&
        c_start >= 0 && c_bath >= 0 && c_start != c_bath &&
        cells >= 2 && dt > 0 && duration > 0 && sample_step >= dt &&
        0 < shell_depth <= radius ||
        throw(ArgumentError("Invalid sphere, storage, film or discretisation inputs."))
    nsteps, skip = round(Int, duration/dt), round(Int, sample_step/dt)
    nsteps * dt ≈ duration && skip * dt ≈ sample_step && nsteps % skip == 0 ||
        throw(ArgumentError("Use time steps that divide both the duration and sampling interval."))

    # The common 4π in integrated fluxes and volumes cancels.
    edges = range(0.0, radius; length=cells+1)
    dr = radius/cells
    vol = diff(edges.^3) ./ 3
    shell_vol = diff(max.(edges, radius-shell_depth).^3) ./ 3
    conductance = d_bulk*eps/tau .* edges[2:end-1].^2 ./ dr
    g_surface = radius^2 / (dr/(2*d_bulk*eps/tau) + 1/kfilm)
    incoming = vcat(0.0, conductance)
    outgoing = vcat(conductance, g_surface)
    # L has units 1/s, acts on kg/m³ liquid and gives the rate of
    # change of the kg/m³ bead inventory.
    L = Tridiagonal(incoming[2:end]./vol[2:end],
                    -(incoming+outgoing)./vol,
                    outgoing[1:end-1]./vol[1:end-1])
    boundary = zeros(cells)
    boundary[end] = g_surface*c_bath/vol[end]
    stored(c) = eps*c + rho_cat*qmax*c/(K+c)
    derivative(c) = eps + rho_cat*qmax*K/(K+c)^2
    c = fill(Float64(c_start), cells)
    initial, equilibrium = stored(c_start), stored(c_bath)
    inventory_scale = abs(equilibrium-initial)*sum(vol)
    tolerance = 5e-12 * max(initial, equilibrium, 1.0)
    lo, hi = minmax(c_start, c_bath)
    remaining, remaining_shell, pore_remaining, times = Float64[], Float64[], Float64[], Float64[]
    function record!(t)
        s = stored.(c) .- equilibrium
        push!(times, t)
        push!(remaining, sum(s.*vol)/(sum(vol)*(initial-equilibrium)))
        push!(remaining_shell, sum(s.*shell_vol)/(sum(shell_vol)*(initial-equilibrium)))
        push!(pore_remaining, sum((c.-c_bath).*vol)/(sum(vol)*(c_start-c_bath)))
    end
    record!(0.0)
    net_exchange = 0.0
    max_balance = 0.0
    max_newton = 0
    for step in 1:nsteps
        old_storage = stored.(c)
        residual(v) = stored.(v) .- old_storage .- dt.*(L*v .+ boundary)
        for i in 1:40
            r = residual(c)
            err = maximum(abs, r)
            max_newton = max(max_newton, i)
            err <= tolerance && break
            J = Tridiagonal(-dt .* L.dl, derivative.(c) .- dt .* L.d, -dt .* L.du)
            update = J \ r
            damping = 1.0
            accepted = false
            while damping >= 2.0^-30
                trial = clamp.(c .- damping.*update, lo, hi)
                if maximum(abs, residual(trial)) < err
                    c .= trial
                    accepted = true
                    break
                end
                damping /= 2
            end
            accepted || error("Could not reduce a moisture-balance residual.")
        end
        maximum(abs, residual(c)) <= tolerance ||
            error("Moisture storage step did not converge.")
        net_exchange += dt*g_surface*(c_bath-c[end])
        imbalance = sum((stored.(c).-initial).*vol) - net_exchange
        max_balance = max(max_balance, abs(imbalance)/inventory_scale)
        step % skip == 0 && record!(step*dt)
    end
    (; times_s=times, remaining, remaining_shell, pore_remaining,
       max_relative_balance_error=max_balance, max_newton_iterations=max_newton)
end

end
