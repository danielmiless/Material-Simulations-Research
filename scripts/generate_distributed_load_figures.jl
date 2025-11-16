"""
Script to generate comparison figures for the distributed load update paper.

This script:
1. Generates a snapshot of the 10×10 system with point load (for comparison)
2. Generates a snapshot of the 11×11 system with distributed load
3. Creates a force distribution visualization
4. Generates comparison snapshots at key time points

Usage: Run this script to generate all figures needed for the distributed load update paper.
"""

using Pkg
Pkg.activate(".")

using CairoMakie
using Printf
using DifferentialEquations
using LinearAlgebra

# Paths
src_dir = joinpath(@__DIR__, "..", "src")
deprecated_dir = joinpath(@__DIR__, "..", "Deprecated Scripts")
figures_dir = joinpath(@__DIR__, "..", "papers", "figures")
mkpath(figures_dir)

println("="^70)
println("Generating Figures for Distributed Load Update Paper")
println("="^70)
println()

# Helper function to load and run simulation
function run_simulation_and_get_snapshot(file_path, title, output_path, snapshot_time=0.2)
    """Load simulation file, run it, and generate a snapshot at specified time."""
    println("Loading and running: $file_path")
    
    # Include the simulation file
    include(file_path)
    
    # Run simulation (suppress output)
    println("  Running simulation...")
    u0 = zeros(2 * TOTAL_DOF)
    tspan = (0.0, T_END)
    prob = ODEProblem(lattice_2d_rhs_with_diagonals_and_backplate!, u0, tspan)
    sol = solve(prob, Vern9();
                reltol = REL_TOL, 
                abstol = ABS_TOL,
                saveat = OUTPUT_INTERVAL,
                dense = false)
    
    if sol.retcode != :Success
        error("Solver failed: $(sol.retcode)")
    end
    
    # Find frame closest to snapshot_time
    frame_idx = argmin(abs.(sol.t .- snapshot_time))
    frame_idx = max(1, min(frame_idx, length(sol.t)))
    actual_time = sol.t[frame_idx]
    
    println("  Generating snapshot at t=$(round(actual_time, digits=3))s (frame $frame_idx)")
    
    # Get equilibrium grid and backplate positions
    equilibrium_grid = Base.invokelatest(create_equilibrium_grid)
    backplate_positions = Base.invokelatest(create_backplate_positions)
    
    # Determine system size
    total_masses = size(equilibrium_grid, 2)
    total_dof = total_masses * 2
    n_size = Int(sqrt(total_masses))
    
    # Create figure
    fig = Figure(size = (1000, 1000), fontsize = 14)
    ax = Axis(fig[1, 1],
              title = title,
              xlabel = "X Position (m)",
              ylabel = "Y Position (m)",
              aspect = DataAspect())
    
    # Get positions at this frame
    pos = reshape(view(sol.u[frame_idx], 1:total_dof), 2, total_masses)
    current_positions = pos .+ equilibrium_grid
    
    # Generate spring connections
    function idx(i, j, n)
        return (i - 1) * n + j
    end
    
    nearest_neighbor_connections = Tuple{Int, Int}[]
    diagonal_connections = Tuple{Int, Int}[]
    
    for i in 1:n_size, j in 1:n_size
        k = idx(i, j, n_size)
        
        # Horizontal springs
        if j < n_size
            push!(nearest_neighbor_connections, (k, idx(i, j+1, n_size)))
        end
        
        # Vertical springs
        if i < n_size
            push!(nearest_neighbor_connections, (k, idx(i+1, j, n_size)))
        end
        
        # Diagonal springs
        if i < n_size && j < n_size
            push!(diagonal_connections, (k, idx(i+1, j+1, n_size)))
        end
        if i < n_size && j > 1
            push!(diagonal_connections, (k, idx(i+1, j-1, n_size)))
        end
    end
    
    # Plot springs
    for (k1, k2) in nearest_neighbor_connections
        spring_points = Point2f[current_positions[:, k1], current_positions[:, k2]]
        lines!(ax, spring_points, color = :blue, linewidth = 2)
    end
    
    for (k1, k2) in diagonal_connections
        spring_points = Point2f[current_positions[:, k1], current_positions[:, k2]]
        lines!(ax, spring_points, color = :red, linewidth = 1.5)
    end
    
    # Plot backplate
    backplate_line_points = Point2f[backplate_positions[:, i] for i in 1:n_size]
    lines!(ax, backplate_line_points, color = :gray, linewidth = 8)
    
    # Plot masses with column-based coloring
    mass_points = Point2f[current_positions[:, k] for k in 1:total_masses]
    colors = Vector{Symbol}(undef, total_masses)
    column_colors = [:lightblue, :lightcyan, :cyan, :blue, :mediumblue, :darkblue, :navy, :midnightblue, :darkslateblue, :darkblue, :black]
    
    for k in 1:total_masses
        i = 1 + div(k-1, n_size)
        j = 1 + mod(k-1, n_size)
        color_idx = min(j, length(column_colors))
        colors[k] = column_colors[color_idx]
    end
    
    scatter!(ax, mass_points, markersize = 12, color = colors, strokewidth = 1.5, strokecolor = :white)
    
    # Add force visualization
    if actual_time <= F_ACTIVE_TIME
        # Check if this is distributed load (11x11) or point load (10x10)
        if n_size == 11
            # Distributed load - show all forces
            fx_unit, fy_unit = Base.invokelatest(calculate_force_components, 1.0, FORCE_ANGLE_DEGREES)
            scale_factor = 0.3
            
            for i in 1:n_size
                force_mag = Base.invokelatest(get_distributed_force_magnitude, i)
                target_idx = idx(i, 1, n_size)
                target_pos = equilibrium_grid[:, target_idx]
                arrow_length = scale_factor * force_mag / F_MAG
                force_end = target_pos + arrow_length * [fx_unit, fy_unit]
                
                lines!(ax, [Point2f(target_pos), Point2f(force_end)], color = :green, linewidth = 3)
                scatter!(ax, [Point2f(force_end)], markersize = 8, color = :green, marker = :circle)
            end
        else
            # Point load - show single force
            if isdefined(Main, :FORCE_TARGET_ROW) && isdefined(Main, :FORCE_TARGET_COL)
                target_idx = idx(FORCE_TARGET_ROW, FORCE_TARGET_COL, n_size)
                fx, fy = Base.invokelatest(calculate_force_components, F_MAG, FORCE_ANGLE_DEGREES)
                scale_factor = 0.3
                target_pos = equilibrium_grid[:, target_idx]
                force_end = target_pos + scale_factor * [fx, fy] / F_MAG
                
                lines!(ax, [Point2f(target_pos), Point2f(force_end)], color = :green, linewidth = 4)
                scatter!(ax, [Point2f(force_end)], markersize = 12, color = :green, marker = :circle)
            end
        end
    end
    
    # Set axis limits
    grid_min = minimum(equilibrium_grid) - 0.5
    grid_max_x = maximum(backplate_positions[1, :]) + 0.5
    grid_max_y = maximum(backplate_positions[2, :]) + 0.5
    limits!(ax, grid_min, grid_max_x, grid_min, grid_max_y)
    
    # Save figure
    save(output_path, fig, px_per_unit = 2)
    println("  Saved to: $output_path")
    println()
    
    return sol, equilibrium_grid, backplate_positions
end

# Generate force distribution diagram
function generate_force_distribution_diagram(output_path, f_mag=1200.0)
    """Generate a diagram showing the force distribution pattern."""
    println("Generating force distribution diagram...")
    
    fig = Figure(size = (800, 600), fontsize = 14)
    ax = Axis(fig[1, 1],
              title = "Distributed Load Pattern (Trapezoidal Distribution)",
              xlabel = "Row Index",
              ylabel = "Force Magnitude (N)",
              xticks = 1:11)
    
    # Force values
    rows = 1:11
    forces = [i == 1 || i == 11 ? f_mag/20.0 : f_mag/10.0 for i in rows]
    
    # Plot as bar chart
    barplot!(ax, rows, forces, color = :blue, strokewidth = 2, strokecolor = :black)
    
    # Add labels
    for (i, f) in enumerate(forces)
        text!(ax, i, f + f_mag/50, text = "F/$(i == 1 || i == 11 ? "20" : "10")", 
              fontsize = 10, align = (:center, :bottom))
    end
    
    # Add horizontal line for reference
    hlines!(ax, [f_mag/20.0, f_mag/10.0], color = :gray, linestyle = :dash, linewidth = 1)
    
    save(output_path, fig, px_per_unit = 2)
    println("  Saved to: $output_path")
    println()
end

# Main execution
println("Step 1: Generating 10×10 point load snapshot (for comparison)...")
point_load_path = joinpath(figures_dir, "point_load_10x10_snapshot.png")
try
    sol_10x10, eq_grid_10, bp_10 = run_simulation_and_get_snapshot(
        joinpath(deprecated_dir, "lattice_simulation_10x10.jl"),
        "10×10 Lattice with Point Load (Comparison)",
        point_load_path,
        0.2
    )
catch e
    println("  Warning: Could not generate 10×10 snapshot: $e")
    println("  (This is okay if you only want 11×11 figures)")
end

println("Step 2: Generating 11×11 distributed load snapshot...")
distributed_load_path = joinpath(figures_dir, "distributed_load_11x11_snapshot.png")
sol_11x11, eq_grid_11, bp_11 = run_simulation_and_get_snapshot(
    joinpath(src_dir, "lattice_simulation_11x11.jl"),
    "11×11 Lattice with Distributed Load",
    distributed_load_path,
    0.2
)

println("Step 3: Generating force distribution diagram...")
force_dist_path = joinpath(figures_dir, "force_distribution_pattern.png")
generate_force_distribution_diagram(force_dist_path, 1200.0)  # Use default F_MAG value

println("="^70)
println("Figure generation complete!")
println("="^70)
println()
println("Generated figures:")
println("  - $point_load_path")
println("  - $distributed_load_path")
println("  - $force_dist_path")
println()

