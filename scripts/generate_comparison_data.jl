"""
Script to generate comparison data and figures for the LaTeX update document.

Usage: Run this script after ensuring both simulation files are available.
It will:
1. Run both simulations (with and without backplate)
2. Extract energy statistics
3. Generate comparison figures
4. Output data in LaTeX table format

Note: This script temporarily modifies the simulation files to prevent
automatic execution, then restores them.
"""

using Pkg
Pkg.activate(".")

using GLMakie
using Printf
using DifferentialEquations
using LinearAlgebra

# Paths
src_dir = joinpath(@__DIR__, "..", "src")
figures_dir = joinpath(@__DIR__, "..", "papers", "figures")
mkpath(figures_dir)

# Read simulation files and extract function definitions
function load_simulation_code()
    """Load simulation code without executing it."""
    sim_file_path = joinpath(src_dir, "lattice_simulation.jl")
    bp_file_path = joinpath(src_dir, "lattice_simulation_with_backplate.jl")
    
    # Read files
    sim_code = read(sim_file_path, String)
    bp_code = read(bp_file_path, String)
    
    # Split at execution marker to get only function definitions
    sim_functions = split(sim_code, "#  RUN SIMULATION")[1]
    bp_functions = split(bp_code, "#  RUN SIMULATION")[1]
    
    # Write temporary files with just the functions
    temp_sim = tempname() * ".jl"
    temp_bp = tempname() * ".jl"
    write(temp_sim, sim_functions)
    write(temp_bp, bp_functions)
    
    try
        # Include the function definitions
        include(temp_sim)
        include(temp_bp)
    finally
        # Clean up temp files
        rm(temp_sim, force=true)
        rm(temp_bp, force=true)
    end
end

println("Loading simulation functions...")
load_simulation_code()

function extract_energy_data(sol, total_work, has_backplate = false)
    """Extract energy statistics from solution."""
    time_values = sol.t
    energy_values = Float64[]
    
    for i in 1:length(sol.t)
        pos = reshape(view(sol.u[i], 1:TOTAL_DOF), 2, TOTAL_MASSES)
        vel = reshape(view(sol.u[i], TOTAL_DOF+1:2*TOTAL_DOF), 2, TOTAL_MASSES)
        
        ke = kinetic_energy_2d(vel)
        # Check which potential energy function is available
        if has_backplate && isdefined(Main, :potential_energy_2d_with_diagonals_and_backplate)
            pe = potential_energy_2d_with_diagonals_and_backplate(pos)
        else
            pe = potential_energy_2d_with_diagonals(pos)
        end
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

function extract_displacement_data(sol)
    """Extract maximum displacement of rightmost column."""
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
    
    return max_x_displacement
end

function generate_energy_comparison_figure(data_no_backplate, data_with_backplate, output_path)
    """Generate energy evolution comparison figure."""
    fig = Figure(size = (1000, 600), fontsize = 14)
    ax = Axis(fig[1, 1],
              title = "Energy Evolution: Comparison with and without Backplate",
              xlabel = "Time (s)",
              ylabel = "Energy (J)")
    
    # Plot energy evolution
    lines!(ax, data_no_backplate["time"], data_no_backplate["energy"],
           color = :blue, linewidth = 2, label = "Without Backplate")
    lines!(ax, data_with_backplate["time"], data_with_backplate["energy"],
           color = :red, linewidth = 2, label = "With Backplate")
    
    # Plot work input lines
    lines!(ax, data_no_backplate["time"], fill(data_no_backplate["total_work"], length(data_no_backplate["time"])),
           color = :blue, linestyle = :dash, linewidth = 1.5, label = "Work Input (No Backplate)")
    lines!(ax, data_with_backplate["time"], fill(data_with_backplate["total_work"], length(data_with_backplate["time"])),
           color = :red, linestyle = :dash, linewidth = 1.5, label = "Work Input (With Backplate)")
    
    axislegend(ax, position = :rt)
    
    save(output_path, fig)
    println("Saved energy comparison figure to: $output_path")
end

function generate_lattice_snapshot(sol, equilibrium_grid, backplate_positions, title, output_path, has_backplate = false, frame_idx = nothing)
    """Generate a snapshot of the lattice at a representative time."""
    # Use a time point near the middle of the simulation (or first frame if only one)
    if frame_idx === nothing
        frame_idx = length(sol.t) > 1 ? div(length(sol.t), 2) : 1
    end
    # Ensure frame_idx is valid (1-indexed)
    frame_idx = max(1, min(frame_idx, length(sol.t)))
    
    fig = Figure(size = (800, 800), fontsize = 14)
    ax = Axis(fig[1, 1],
              title = title,
              xlabel = "X Position (m)",
              ylabel = "Y Position (m)",
              aspect = DataAspect())
    
    # Get positions at this frame
    pos = reshape(view(sol.u[frame_idx], 1:TOTAL_DOF), 2, TOTAL_MASSES)
    current_positions = pos .+ equilibrium_grid
    
    # Get spring connections
    nearest_neighbor_connections, diagonal_connections = create_spring_connections_with_diagonals()
    
    # Plot springs
    for (k1, k2) in nearest_neighbor_connections
        spring_points = Point2f[current_positions[:, k1], current_positions[:, k2]]
        lines!(ax, spring_points, color = :blue, linewidth = 2)
    end
    
    for (k1, k2) in diagonal_connections
        spring_points = Point2f[current_positions[:, k1], current_positions[:, k2]]
        lines!(ax, spring_points, color = :red, linewidth = 1.5)
    end
    
    # Plot backplate if present
    if has_backplate && backplate_positions !== nothing
        backplate_line_points = Point2f[backplate_positions[:, i] for i in 1:N]
        lines!(ax, backplate_line_points, color = :gray, linewidth = 8)
    end
    
    # Plot masses
    mass_points = Point2f[current_positions[:, k] for k in 1:TOTAL_MASSES]
    colors = fill(:black, TOTAL_MASSES)
    target_idx = lattice_idx(FORCE_TARGET_ROW, FORCE_TARGET_COL)
    colors[target_idx] = :orange
    scatter!(ax, mass_points, markersize = 15, color = colors, strokewidth = 2, strokecolor = :white)
    
    # Set axis limits
    grid_min = minimum(equilibrium_grid) - 0.5
    if has_backplate && backplate_positions !== nothing
        grid_max_x = maximum(backplate_positions[1, :]) + 0.5
        grid_max_y = maximum(backplate_positions[2, :]) + 0.5
    else
        grid_max_x = maximum(equilibrium_grid[1, :]) + 0.5
        grid_max_y = maximum(equilibrium_grid[2, :]) + 0.5
    end
    limits!(ax, grid_min, grid_max_x, grid_min, grid_max_y)
    
    save(output_path, fig)
    println("Saved snapshot to: $output_path")
end

function print_latex_table_data(data_no_backplate, data_with_backplate, max_disp_no_backplate, max_disp_with_backplate)
    """Print data formatted for LaTeX table."""
    println("\n" * "="^70)
    println("DATA FOR LaTeX TABLE")
    println("="^70)
    println()
    println("Copy this into the table in your LaTeX document:")
    println()
    println("\\begin{tabular}{|l|c|c|}")
    println("\\hline")
    println("\\textbf{Metric} & \\textbf{Without Backplate} & \\textbf{With Backplate} \\\\")
    println("\\hline")
    @printf("Total Work Input (J) & %.6f & %.6f \\\\\n", 
            data_no_backplate["total_work"], data_with_backplate["total_work"])
    @printf("Final Total Energy (J) & %.6f & %.6f \\\\\n", 
            data_no_backplate["final_energy"], data_with_backplate["final_energy"])
    @printf("Energy Dissipated (J) & %.6f & %.6f \\\\\n", 
            data_no_backplate["energy_dissipated"], data_with_backplate["energy_dissipated"])
    @printf("Dissipation Percentage (\\%%) & %.2f & %.2f \\\\\n", 
            data_no_backplate["dissipation_percent"], data_with_backplate["dissipation_percent"])
    println("\\hline")
    println("\\end{tabular}")
    println()
    println("Maximum x-displacement of rightmost column:")
    @printf("  Without backplate: %.6f m\n", max_disp_no_backplate)
    @printf("  With backplate: %.6f m\n", max_disp_with_backplate)
    println()
end

# Main execution
println("="^70)
println("Generating Comparison Data and Figures")
println("="^70)
println()

println("Running simulation WITHOUT backplate...")
println("(This will display the animation window - you can close it)")
fig_no_bp, sol_no_bp, work_no_bp, energy_no_bp = run_2d_simulation_with_configurable_force();
data_no_backplate = extract_energy_data(sol_no_bp, work_no_bp, false)
max_disp_no_backplate = extract_displacement_data(sol_no_bp)
equilibrium_grid = create_equilibrium_grid()

println("\nRunning simulation WITH backplate...")
println("(This will display the animation window - you can close it)")
fig_with_bp, sol_with_bp, work_with_bp, energy_with_bp = run_2d_simulation_with_configurable_force_and_backplate();
data_with_backplate = extract_energy_data(sol_with_bp, work_with_bp, true)
max_disp_with_backplate = extract_displacement_data(sol_with_bp)
backplate_positions = create_backplate_positions()

println("\nGenerating figures...")

# Generate energy comparison figure
energy_fig_path = joinpath(figures_dir, "energy_comparison_damping.png")
generate_energy_comparison_figure(data_no_backplate, data_with_backplate, energy_fig_path)

# Generate initial setup figure (with backplate at t=0)
println("\nGenerating initial setup figure...")
# Use frame 1 (initial state) from the backplate simulation
initial_fig_path = joinpath(figures_dir, "backplate_initial_setup.png")
generate_lattice_snapshot(sol_with_bp, equilibrium_grid, backplate_positions, 
                         "Initial Configuration with Backplate", initial_fig_path, true, 1)

# Generate snapshot without backplate
snapshot_no_bp_path = joinpath(figures_dir, "no_backplate_snapshot.png")
mid_time = round(sol_no_bp.t[div(length(sol_no_bp.t), 2)], digits=2)
generate_lattice_snapshot(sol_no_bp, equilibrium_grid, nothing, 
                         "Lattice Configuration Without Backplate (t = $mid_time s)", 
                         snapshot_no_bp_path, false)

# Generate snapshot with backplate
snapshot_with_bp_path = joinpath(figures_dir, "with_backplate_snapshot.png")
mid_time_bp = round(sol_with_bp.t[div(length(sol_with_bp.t), 2)], digits=2)
generate_lattice_snapshot(sol_with_bp, equilibrium_grid, backplate_positions, 
                         "Lattice Configuration With Backplate (t = $mid_time_bp s)", 
                         snapshot_with_bp_path, true)

# Print data for LaTeX table
print_latex_table_data(data_no_backplate, data_with_backplate, max_disp_no_backplate, max_disp_with_backplate)

println("="^70)
println("Done! All figures saved to: $figures_dir")
println("="^70)
