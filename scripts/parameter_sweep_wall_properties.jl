"""
Script to perform parameter sweep for wall properties (WALL_STIFFNESS and WALL_DAMPING).

This script runs multiple simulations with different wall stiffness and damping values
to analyze the sensitivity of the system to these parameters.

Usage: Run this script to generate parameter sweep results.
"""

using Pkg
Pkg.activate(".")

using GLMakie
using Printf
using DifferentialEquations
using LinearAlgebra
using Dates

# Paths
src_dir = joinpath(@__DIR__, "..", "src")
figures_dir = joinpath(@__DIR__, "..", "papers", "figures")
data_dir = joinpath(@__DIR__, "..", "data")
mkpath(figures_dir)
mkpath(data_dir)

println("Loading 11×11 simulation (geometric NN springs + backplate)...")
include(joinpath(src_dir, "lattice_simulation_11x11.jl"))

# Base values
const BASE_WALL_STIFFNESS = 10000.0
const BASE_WALL_DAMPING = 10.0

# Multiplier values to test
const MULTIPLIERS = [0.5, 1.0, 2.0, 5.0, 10.0]

function create_modified_rhs(wall_stiffness, wall_damping)
    """
    Create a modified ODE RHS matching `lattice_2d_rhs_nn_backplate!` but with custom wall (k, c).
    """
    function modified_rhs!(du, u, p, t)
        pos = reshape(view(u, 1:TOTAL_DOF), 2, TOTAL_MASSES)
        vel = reshape(view(u, TOTAL_DOF+1:2*TOTAL_DOF), 2, TOTAL_MASSES)
        dpos = reshape(view(du, 1:TOTAL_DOF), 2, TOTAL_MASSES)
        dvel = reshape(view(du, TOTAL_DOF+1:2*TOTAL_DOF), 2, TOTAL_MASSES)
        dpos .= vel
        dvel .= 0.0
        material_order = collect(1:N)
        materials = DEFAULT_MATERIALS
        Δ = GRID_SPACING
        for k in 1:TOTAL_MASSES
            i, j = lattice_i(k), lattice_j(k)
            pos_k = pos[:, k]
            vel_k = vel[:, k]
            function add_spring_force!(neighbor_i, neighbor_j, R0, k_spring, alpha_spring, add_damping=false, damping_coeff=C_DAMPING)
                if 1 <= neighbor_i <= N && 1 <= neighbor_j <= N
                    neighbor_idx = lattice_idx(neighbor_i, neighbor_j)
                    u_neighbor = pos[:, neighbor_idx]
                    vel_neighbor = vel[:, neighbor_idx]
                    force = spring_force_2d_geometric(pos_k, u_neighbor, R0, k_spring, alpha_spring)
                    dvel[:, k] += force
                    if add_damping
                        dvel[:, k] += damping_force_2d(vel_k, vel_neighbor, damping_coeff)
                    end
                end
            end
            add_spring_force!(i, j-1, [-Δ, 0.0], K_COUPLING, ALPHA_COUPLING, true, C_DAMPING)
            add_spring_force!(i, j+1, [Δ, 0.0], K_COUPLING, ALPHA_COUPLING, true, C_DAMPING)
            k_col = get_column_k_coupling(j, material_order, materials)
            c_col = get_column_c_damping(j, material_order, materials)
            alpha_col = get_column_alpha_coupling(j, material_order, materials)
            add_spring_force!(i-1, j, [0.0, -Δ], k_col, alpha_col, true, c_col)
            add_spring_force!(i+1, j, [0.0, Δ], k_col, alpha_col, true, c_col)
            if j == N
                equilibrium_x = (N - 1) * GRID_SPACING
                backplate_x = equilibrium_x + BACKPLATE_DISTANCE
                current_x = equilibrium_x + pos_k[1]
                if current_x > backplate_x
                    penetration = current_x - backplate_x
                    dvel[1, k] += -wall_stiffness * penetration
                    if vel_k[1] > 0
                        dvel[1, k] -= wall_damping * vel_k[1]
                    end
                end
            end
        end
        if t <= F_ACTIVE_TIME
            fx_unit, fy_unit = calculate_force_components(1.0, FORCE_ANGLE_DEGREES)
            for row in 1:N
                force_mag = get_distributed_force_magnitude(row)
                target_idx = lattice_idx(row, 1)
                dvel[1, target_idx] += force_mag * fx_unit
                dvel[2, target_idx] += force_mag * fy_unit
            end
        end
        dvel ./= MASS
        return nothing
    end
    return modified_rhs!
end

function potential_energy_with_wall_params(pos_matrix, wall_stiffness)
    """
    Geometric NN spring PE (same as `potential_energy_2d_nn_backplate`) plus wall penalty with given stiffness.
    """
    pe = 0.0
    material_order = collect(1:N)
    materials = DEFAULT_MATERIALS
    Δ = GRID_SPACING
    function add_pe(uA, uB, R0, k_spring, alpha_spring)
        R = R0 .+ (uB .- uA)
        L = norm(R)
        L0 = norm(R0)
        s = L - L0
        return (k_spring / alpha_spring) * (exp(alpha_spring * s) - alpha_spring * s - 1.0)
    end
    for i in 1:N, j in 1:N-1
        k1 = lattice_idx(i, j)
        k2 = lattice_idx(i, j+1)
        pe += add_pe(pos_matrix[:, k1], pos_matrix[:, k2], [Δ, 0.0], K_COUPLING, ALPHA_COUPLING)
    end
    for i in 1:N-1, j in 1:N
        k1 = lattice_idx(i, j)
        k2 = lattice_idx(i+1, j)
        k_j = get_column_k_coupling(j, material_order, materials)
        a_j = get_column_alpha_coupling(j, material_order, materials)
        pe += add_pe(pos_matrix[:, k1], pos_matrix[:, k2], [0.0, Δ], k_j, a_j)
    end
    # Wall potential energy (penalty for penetration)
    equilibrium_x = (N - 1) * GRID_SPACING
    backplate_x = equilibrium_x + BACKPLATE_DISTANCE
    for i in 1:N
        k = lattice_idx(i, N)  # Rightmost column
        current_x = equilibrium_x + pos_matrix[1, k]  # absolute x position
        if current_x > backplate_x
            penetration = current_x - backplate_x
            # Penalty potential energy (like a very stiff spring)
            pe += 0.5 * wall_stiffness * penetration^2
        end
    end
    
    return pe
end

function run_simulation_with_wall_params(wall_stiffness, wall_damping)
    """
    Run simulation with specified wall parameters.
    Returns: (solution, total_work, final_energy, max_displacement)
    """
    # Create modified RHS function with specified wall parameters
    modified_rhs! = create_modified_rhs(wall_stiffness, wall_damping)
    
    # Initial state: all masses at rest at equilibrium positions
    u0 = zeros(2 * TOTAL_DOF)
    tspan = (0.0, T_END)
    
    # Create and solve ODE problem
    prob = ODEProblem(modified_rhs!, u0, tspan)
    sol = solve(prob, Vern9();
                reltol = REL_TOL, 
                abstol = ABS_TOL,
                saveat = OUTPUT_INTERVAL,
                dense = false)
    
    if sol.retcode != :Success
        error("Solver failed: $(sol.retcode)")
    end
    
    total_work = work_done_2d_distributed(sol)
    
    # Calculate final energy
    pos_final = reshape(view(sol.u[end], 1:TOTAL_DOF), 2, TOTAL_MASSES)
    vel_final = reshape(view(sol.u[end], TOTAL_DOF+1:2*TOTAL_DOF), 2, TOTAL_MASSES)
    final_energy = kinetic_energy_2d(vel_final) + potential_energy_with_wall_params(pos_final, wall_stiffness)
    
    # Calculate maximum displacement of rightmost column
    equilibrium_grid = create_equilibrium_grid()
    max_x_displacement = 0.0
    
    for i in 1:length(sol.t)
        pos = reshape(view(sol.u[i], 1:TOTAL_DOF), 2, TOTAL_MASSES)
        current_positions = pos .+ equilibrium_grid
        
        # Find maximum x-displacement of rightmost column (j = N)
        for row in 1:N
            k = lattice_idx(row, N)
            x_pos = current_positions[1, k]
            equilibrium_x = (N - 1) * GRID_SPACING
            displacement = x_pos - equilibrium_x
            max_x_displacement = max(max_x_displacement, displacement)
        end
    end
    
    return sol, total_work, final_energy, max_x_displacement
end

function extract_energy_data(sol, total_work, wall_stiffness)
    """Extract energy statistics from solution."""
    time_values = sol.t
    energy_values = Float64[]
    
    for i in 1:length(sol.t)
        pos = reshape(view(sol.u[i], 1:TOTAL_DOF), 2, TOTAL_MASSES)
        vel = reshape(view(sol.u[i], TOTAL_DOF+1:2*TOTAL_DOF), 2, TOTAL_MASSES)
        
        ke = kinetic_energy_2d(vel)
        pe = potential_energy_with_wall_params(pos, wall_stiffness)
        push!(energy_values, ke + pe)
    end
    
    final_energy = energy_values[end]
    energy_dissipated = total_work - final_energy
    dissipation_percent = total_work > 0 ? (energy_dissipated / total_work * 100) : 0.0
    
    return Dict(
        "time" => time_values,
        "energy" => energy_values,
        "total_work" => total_work,
        "final_energy" => final_energy,
        "energy_dissipated" => energy_dissipated,
        "dissipation_percent" => dissipation_percent
    )
end

# Main parameter sweep
println("="^70)
println("Parameter Sweep: Wall Properties")
println("="^70)
println()

results = []

for mult_stiff in MULTIPLIERS
    for mult_damp in MULTIPLIERS
        wall_stiff = BASE_WALL_STIFFNESS * mult_stiff
        wall_damp = BASE_WALL_DAMPING * mult_damp
        
        println("Running simulation: WALL_STIFFNESS = $(wall_stiff) N/m ($(mult_stiff)×), WALL_DAMPING = $(wall_damp) N·s/m ($(mult_damp)×)")
        
        try
            sol, total_work, final_energy, max_disp = run_simulation_with_wall_params(wall_stiff, wall_damp)
            energy_data = extract_energy_data(sol, total_work, wall_stiff)
            
            push!(results, Dict(
                "wall_stiffness" => wall_stiff,
                "wall_stiffness_mult" => mult_stiff,
                "wall_damping" => wall_damp,
                "wall_damping_mult" => mult_damp,
                "total_work" => total_work,
                "final_energy" => final_energy,
                "energy_dissipated" => energy_data["energy_dissipated"],
                "dissipation_percent" => energy_data["dissipation_percent"],
                "max_displacement" => max_disp,
                "solution" => sol,
                "energy_data" => energy_data
            ))
            
            println("  ✓ Completed: Dissipation = $(round(energy_data["dissipation_percent"], digits=2))%, Max displacement = $(round(max_disp, digits=6)) m")
        catch e
            println("  ✗ Failed: $e")
        end
        
        println()
    end
end

# Generate comparison plots
println("Generating comparison plots...")

# Plot 1: Energy dissipation vs wall stiffness (for different damping values)
fig1 = Figure(size = (1200, 800), fontsize = 14)

# Subplot 1: Dissipation percentage vs wall stiffness
ax1 = Axis(fig1[1, 1],
           title = "Energy Dissipation vs Wall Stiffness",
           xlabel = "Wall Stiffness Multiplier",
           ylabel = "Dissipation Percentage (%)")

for mult_damp in MULTIPLIERS
    stiffness_vals = Float64[]
    dissipation_vals = Float64[]
    
    for r in results
        if r["wall_damping_mult"] == mult_damp
            push!(stiffness_vals, r["wall_stiffness_mult"])
            push!(dissipation_vals, r["dissipation_percent"])
        end
    end
    
    if !isempty(stiffness_vals)
        # Sort by stiffness
        sorted_indices = sortperm(stiffness_vals)
        stiffness_vals = stiffness_vals[sorted_indices]
        dissipation_vals = dissipation_vals[sorted_indices]
        
        scatterlines!(ax1, stiffness_vals, dissipation_vals, 
                      label = "Damping = $(mult_damp)×", 
                      linewidth = 2,
                      markersize = 8)
    end
end

axislegend(ax1, position = :rt)

# Subplot 2: Max displacement vs wall stiffness
ax2 = Axis(fig1[1, 2],
           title = "Maximum Displacement vs Wall Stiffness",
           xlabel = "Wall Stiffness Multiplier",
           ylabel = "Max Displacement (m)")

for mult_damp in MULTIPLIERS
    stiffness_vals = Float64[]
    disp_vals = Float64[]
    
    for r in results
        if r["wall_damping_mult"] == mult_damp
            push!(stiffness_vals, r["wall_stiffness_mult"])
            push!(disp_vals, r["max_displacement"])
        end
    end
    
    if !isempty(stiffness_vals)
        sorted_indices = sortperm(stiffness_vals)
        stiffness_vals = stiffness_vals[sorted_indices]
        disp_vals = disp_vals[sorted_indices]
        
        scatterlines!(ax2, stiffness_vals, disp_vals,
                      label = "Damping = $(mult_damp)×",
                      linewidth = 2,
                      markersize = 8)
    end
end

axislegend(ax2, position = :rt)

# Subplot 3: Dissipation percentage vs wall damping
ax3 = Axis(fig1[2, 1],
           title = "Energy Dissipation vs Wall Damping",
           xlabel = "Wall Damping Multiplier",
           ylabel = "Dissipation Percentage (%)")

for mult_stiff in MULTIPLIERS
    damping_vals = Float64[]
    dissipation_vals = Float64[]
    
    for r in results
        if r["wall_stiffness_mult"] == mult_stiff
            push!(damping_vals, r["wall_damping_mult"])
            push!(dissipation_vals, r["dissipation_percent"])
        end
    end
    
    if !isempty(damping_vals)
        sorted_indices = sortperm(damping_vals)
        damping_vals = damping_vals[sorted_indices]
        dissipation_vals = dissipation_vals[sorted_indices]
        
        scatterlines!(ax3, damping_vals, dissipation_vals,
                      label = "Stiffness = $(mult_stiff)×",
                      linewidth = 2,
                      markersize = 8)
    end
end

axislegend(ax3, position = :rt)

# Subplot 4: Max displacement vs wall damping
ax4 = Axis(fig1[2, 2],
           title = "Maximum Displacement vs Wall Damping",
           xlabel = "Wall Damping Multiplier",
           ylabel = "Max Displacement (m)")

for mult_stiff in MULTIPLIERS
    damping_vals = Float64[]
    disp_vals = Float64[]
    
    for r in results
        if r["wall_stiffness_mult"] == mult_stiff
            push!(damping_vals, r["wall_damping_mult"])
            push!(disp_vals, r["max_displacement"])
        end
    end
    
    if !isempty(damping_vals)
        sorted_indices = sortperm(damping_vals)
        damping_vals = damping_vals[sorted_indices]
        disp_vals = disp_vals[sorted_indices]
        
        scatterlines!(ax4, damping_vals, disp_vals,
                      label = "Stiffness = $(mult_stiff)×",
                      linewidth = 2,
                      markersize = 8)
    end
end

axislegend(ax4, position = :rt)

save(joinpath(figures_dir, "wall_parameter_sweep.png"), fig1)
println("Saved parameter sweep plots to: $(joinpath(figures_dir, "wall_parameter_sweep.png"))")

# Print LaTeX table
println("\n" * "="^70)
println("LaTeX TABLE FOR PARAMETER SWEEP RESULTS")
println("="^70)
println()
println("\\begin{table}[h]")
println("\\centering")
println("\\begin{tabular}{|c|c|c|c|c|c|}")
println("\\hline")
println("\\textbf{Stiffness} & \\textbf{Damping} & \\textbf{Work (J)} & \\textbf{Final Energy (J)} & \\textbf{Dissipation (\\%)} & \\textbf{Max Disp. (m)} \\\\")
println("\\hline")

for r in results
    @printf("%.1f× & %.1f× & %.6f & %.6f & %.2f & %.6f \\\\\n",
            r["wall_stiffness_mult"], r["wall_damping_mult"],
            r["total_work"], r["final_energy"],
            r["dissipation_percent"], r["max_displacement"])
end

println("\\hline")
println("\\end{tabular}")
println("\\caption{Parameter sweep results for wall stiffness and damping.}")
println("\\label{tab:wall_parameter_sweep}")
println("\\end{table}")

# Save results to file
results_file = joinpath(data_dir, "wall_parameter_sweep_results.jl")
open(results_file, "w") do f
    println(f, "# Parameter sweep results for wall properties")
    println(f, "# Generated on: $(Dates.now())")
    println(f)
    println(f, "results = [")
    for r in results
        println(f, "    Dict(")
        println(f, "        \"wall_stiffness\" => $(r["wall_stiffness"]),")
        println(f, "        \"wall_stiffness_mult\" => $(r["wall_stiffness_mult"]),")
        println(f, "        \"wall_damping\" => $(r["wall_damping"]),")
        println(f, "        \"wall_damping_mult\" => $(r["wall_damping_mult"]),")
        println(f, "        \"total_work\" => $(r["total_work"]),")
        println(f, "        \"final_energy\" => $(r["final_energy"]),")
        println(f, "        \"energy_dissipated\" => $(r["energy_dissipated"]),")
        println(f, "        \"dissipation_percent\" => $(r["dissipation_percent"]),")
        println(f, "        \"max_displacement\" => $(r["max_displacement"]),")
        println(f, "    ),")
    end
    println(f, "]")
end

println("\nSaved results to: $results_file")
println("="^70)
println("Parameter sweep completed!")
println("="^70)

