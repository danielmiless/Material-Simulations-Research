"""
Script to generate all materials (figures, data, tables) for the comprehensive update paper.

This script:
1. Runs 5×5 simulations with 1.5× material multiplier
2. Runs 5×5 simulations with 2/3× material multiplier  
3. Runs 10×10 simulation with material scaling
4. Generates comparison figures
5. Extracts energy statistics and displacement data
6. Creates LaTeX tables
7. Generates snapshot figures at key time points

Usage: Run this script to generate all materials needed for the update paper.
"""

using Pkg
Pkg.activate(".")

# Use CairoMakie for static figure generation (non-interactive, better for saving files)
using CairoMakie
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

println("="^70)
println("Generating Materials for Update Paper")
println("="^70)
println()

# Helper function to load and modify simulation code
function load_and_modify_simulation(file_path, material_multiplier=nothing, n_size=nothing)
    """Load simulation code and optionally modify constants."""
    code = read(file_path, String)
    
    # Split at execution marker
    functions = split(code, "#  RUN SIMULATION")[1]
    
    # Modify constants if needed
    if material_multiplier !== nothing
        # Replace MATERIAL_MULTIPLIER constant
        functions = replace(functions, 
            r"const MATERIAL_MULTIPLIER = [0-9.]+" => 
            "const MATERIAL_MULTIPLIER = $(material_multiplier)")
    end
    
    if n_size !== nothing
        # Replace N constant
        functions = replace(functions,
            r"const N = [0-9]+" =>
            "const N = $(n_size)")
        # Update TOTAL_MASSES and TOTAL_DOF
        total_masses = n_size * n_size
        total_dof = total_masses * 2
        functions = replace(functions,
            r"const TOTAL_MASSES = N \* N[^\n]*" =>
            "const TOTAL_MASSES = $(total_masses)")
        functions = replace(functions,
            r"const TOTAL_DOF = TOTAL_MASSES \* DOF_PER_MASS[^\n]*" =>
            "const TOTAL_DOF = $(total_dof)")
    end
    
    # Write temporary file
    temp_file = tempname() * ".jl"
    write(temp_file, functions)
    
    try
        include(temp_file)
    finally
        rm(temp_file, force=true)
    end
end

# Function to run simulation and extract data
function run_simulation_and_extract_data(sim_type, material_multiplier=nothing, n_size=5)
    """
    Run simulation and extract data.
    sim_type: "backplate" or "10x10"
    """
    println("Running simulation: $(sim_type), multiplier=$(material_multiplier), N=$(n_size)")
    
    if sim_type == "backplate"
        file_path = joinpath(src_dir, "lattice_simulation_with_backplate.jl")
        load_and_modify_simulation(file_path, material_multiplier, n_size)
        
        # Run simulation
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
        
        # Use invokelatest to handle world age issues with dynamically loaded functions
        total_work = Base.invokelatest(work_done_2d_configurable, sol)
        equilibrium_grid = Base.invokelatest(create_equilibrium_grid)
        
        # Extract energy data
        time_values = sol.t
        energy_values = Float64[]
        for i in 1:length(sol.t)
            pos = reshape(view(sol.u[i], 1:TOTAL_DOF), 2, TOTAL_MASSES)
            vel = reshape(view(sol.u[i], TOTAL_DOF+1:2*TOTAL_DOF), 2, TOTAL_MASSES)
            ke = Base.invokelatest(kinetic_energy_2d, vel)
            pe = Base.invokelatest(potential_energy_2d_with_diagonals_and_backplate, pos)
            push!(energy_values, ke + pe)
        end
        
        final_energy = energy_values[end]
        energy_dissipated = total_work - final_energy
        dissipation_percent = total_work > 0 ? (energy_dissipated / total_work * 100) : 0.0
        
        # Extract max displacement
        max_x_displacement = 0.0
        for i in 1:length(sol.t)
            pos = reshape(view(sol.u[i], 1:TOTAL_DOF), 2, TOTAL_MASSES)
            current_positions = pos .+ equilibrium_grid
            for row in 1:N
                k = Base.invokelatest(lattice_idx, row, N)
                x_pos = current_positions[1, k]
                equilibrium_x = (N - 1) * GRID_SPACING
                displacement = x_pos - equilibrium_x
                max_x_displacement = max(max_x_displacement, displacement)
            end
        end
        
        return Dict(
            "solution" => sol,
            "total_work" => total_work,
            "final_energy" => final_energy,
            "energy_dissipated" => energy_dissipated,
            "dissipation_percent" => dissipation_percent,
            "max_displacement" => max_x_displacement,
            "time" => time_values,
            "energy" => energy_values,
            "equilibrium_grid" => equilibrium_grid
        )
        
    elseif sim_type == "10x10"
        file_path = joinpath(src_dir, "lattice_simulation_10x10.jl")
        load_and_modify_simulation(file_path, material_multiplier, 10)
        
        # Run simulation
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
        
        # Use invokelatest to handle world age issues with dynamically loaded functions
        total_work = Base.invokelatest(work_done_2d_configurable, sol)
        equilibrium_grid = Base.invokelatest(create_equilibrium_grid)
        
        # Extract energy data
        time_values = sol.t
        energy_values = Float64[]
        for i in 1:length(sol.t)
            pos = reshape(view(sol.u[i], 1:TOTAL_DOF), 2, TOTAL_MASSES)
            vel = reshape(view(sol.u[i], TOTAL_DOF+1:2*TOTAL_DOF), 2, TOTAL_MASSES)
            ke = Base.invokelatest(kinetic_energy_2d, vel)
            pe = Base.invokelatest(potential_energy_2d_with_diagonals_and_backplate, pos)
            push!(energy_values, ke + pe)
        end
        
        final_energy = energy_values[end]
        energy_dissipated = total_work - final_energy
        dissipation_percent = total_work > 0 ? (energy_dissipated / total_work * 100) : 0.0
        
        return Dict(
            "solution" => sol,
            "total_work" => total_work,
            "final_energy" => final_energy,
            "energy_dissipated" => energy_dissipated,
            "dissipation_percent" => dissipation_percent,
            "time" => time_values,
            "energy" => energy_values,
            "equilibrium_grid" => equilibrium_grid
        )
    end
end

# Function to generate snapshot figure
function generate_snapshot(sol, equilibrium_grid, backplate_positions, title, output_path, frame_idx=nothing)
    """Generate a snapshot of the lattice at a representative time."""
    if frame_idx === nothing
        frame_idx = length(sol.t) > 1 ? div(length(sol.t), 2) : 1
    end
    frame_idx = max(1, min(frame_idx, length(sol.t)))
    
    # Determine system size from equilibrium_grid (2 × TOTAL_MASSES)
    total_masses = size(equilibrium_grid, 2)
    total_dof = total_masses * 2  # 2 DOF per mass (x, y)
    n_size = Int(sqrt(total_masses))  # Assuming square lattice
    
    fig = Figure(size = (800, 800), fontsize = 14)
    ax = Axis(fig[1, 1],
              title = title,
              xlabel = "X Position (m)",
              ylabel = "Y Position (m)",
              aspect = DataAspect())
    
    # Get positions at this frame
    pos = reshape(view(sol.u[frame_idx], 1:total_dof), 2, total_masses)
    current_positions = pos .+ equilibrium_grid
    
    # Generate spring connections based on actual system size (not global constants)
    # Helper function to convert (i, j) to linear index
    function idx(i, j, n)
        return (i - 1) * n + j
    end
    
    nearest_neighbor_connections = Tuple{Int, Int}[]
    diagonal_connections = Tuple{Int, Int}[]
    
    # Generate connections for N×N lattice
    for i in 1:n_size, j in 1:n_size
        k = idx(i, j, n_size)
        
        # Nearest neighbors (horizontal and vertical)
        if j < n_size  # right neighbor
            push!(nearest_neighbor_connections, (k, idx(i, j+1, n_size)))
        end
        if i < n_size  # down neighbor
            push!(nearest_neighbor_connections, (k, idx(i+1, j, n_size)))
        end
        
        # Diagonal neighbors (upper-left to lower-right)
        if i < n_size && j < n_size
            push!(diagonal_connections, (k, idx(i+1, j+1, n_size)))
        end
        # Diagonal neighbors (upper-right to lower-left)
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
    
    # Plot backplate if present
    if backplate_positions !== nothing
        backplate_n = size(backplate_positions, 2)
        backplate_line_points = Point2f[backplate_positions[:, i] for i in 1:backplate_n]
        lines!(ax, backplate_line_points, color = :gray, linewidth = 8)
    end
    
    # Plot masses
    mass_points = Point2f[current_positions[:, k] for k in 1:total_masses]
    colors = fill(:black, total_masses)
    # Try to highlight target mass if we can determine it
    # Note: We skip this if constants aren't available (e.g., when generating snapshots
    # for data from a different simulation size)
    try
        # Constants should be in scope from the last include, but may not match this solution
        if isdefined(Main, :FORCE_TARGET_ROW) && isdefined(Main, :FORCE_TARGET_COL) && isdefined(Main, :lattice_idx)
            target_row = Main.FORCE_TARGET_ROW
            target_col = Main.FORCE_TARGET_COL
            target_idx = Base.invokelatest(lattice_idx, target_row, target_col)
            if 1 <= target_idx <= total_masses
                colors[target_idx] = :orange
            end
        end
    catch
        # If we can't get target index, just use all black
    end
    scatter!(ax, mass_points, markersize = 15, color = colors, strokewidth = 2, strokecolor = :white)
    
    # Set axis limits
    grid_min = minimum(equilibrium_grid) - 0.5
    if backplate_positions !== nothing
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

# Main execution
println("Step 1: Running 5×5 simulation with 1.5× material multiplier...")
data_1_5 = run_simulation_and_extract_data("backplate", 1.5, 5)
backplate_positions_1_5 = Base.invokelatest(create_backplate_positions)
println("  ✓ Completed: Dissipation = $(round(data_1_5["dissipation_percent"], digits=2))%")
println()

println("Step 2: Running 5×5 simulation with 2/3× material multiplier...")
data_2_3 = run_simulation_and_extract_data("backplate", 2/3, 5)
backplate_positions_2_3 = Base.invokelatest(create_backplate_positions)
println("  ✓ Completed: Dissipation = $(round(data_2_3["dissipation_percent"], digits=2))%")
println()

println("Step 3: Running 10×10 simulation with 1.5× material multiplier...")
data_10x10 = run_simulation_and_extract_data("10x10", 1.5, 10)
backplate_positions_10x10 = Base.invokelatest(create_backplate_positions)
println("  ✓ Completed: Dissipation = $(round(data_10x10["dissipation_percent"], digits=2))%")
println()

# Generate comparison figure for material scaling
println("Step 4: Generating material scaling comparison figure...")
fig = Figure(size = (1200, 600), fontsize = 14)
ax = Axis(fig[1, 1],
          title = "Energy Evolution: Material Scaling Comparison",
          xlabel = "Time (s)",
          ylabel = "Energy (J)")

# Use lines with different styles and occasional markers to make differences more visible
# Plot lines first
lines!(ax, data_1_5["time"], data_1_5["energy"],
       color = :blue, linewidth = 3, label = "1.5× Multiplier", linestyle = :solid)
lines!(ax, data_2_3["time"], data_2_3["energy"],
       color = :red, linewidth = 3, label = "2/3× Multiplier", linestyle = :dash)

# Add markers at regular intervals (every 20th point) to distinguish lines without overwhelming the plot
marker_interval = max(1, div(length(data_1_5["time"]), 30))  # ~30 markers total
scatter!(ax, data_1_5["time"][1:marker_interval:end], data_1_5["energy"][1:marker_interval:end],
         color = :blue, markersize = 10, marker = :circle, 
         strokewidth = 1, strokecolor = :white)
scatter!(ax, data_2_3["time"][1:marker_interval:end], data_2_3["energy"][1:marker_interval:end],
         color = :red, markersize = 10, marker = :rect,
         strokewidth = 1, strokecolor = :white)

work_line_1_5 = fill(data_1_5["total_work"], length(data_1_5["time"]))
work_line_2_3 = fill(data_2_3["total_work"], length(data_2_3["time"]))
lines!(ax, data_1_5["time"], work_line_1_5,
       color = :blue, linestyle = :dash, linewidth = 1.5, label = "Work Input (1.5×)")
lines!(ax, data_2_3["time"], work_line_2_3,
       color = :red, linestyle = :dash, linewidth = 1.5, label = "Work Input (2/3×)")

axislegend(ax, position = :rt)
save(joinpath(figures_dir, "material_scaling_comparison.png"), fig)
println("  ✓ Saved: material_scaling_comparison.png")
println()

# Generate snapshot figures
println("Step 5: Generating snapshot figures...")
generate_snapshot(data_1_5["solution"], data_1_5["equilibrium_grid"], backplate_positions_1_5,
                 "5×5 Lattice with 1.5× Material Scaling (t = $(round(data_1_5["time"][div(length(data_1_5["time"]), 2)], digits=2)) s)",
                 joinpath(figures_dir, "material_scaling_1_5_snapshot.png"))

generate_snapshot(data_2_3["solution"], data_2_3["equilibrium_grid"], backplate_positions_2_3,
                 "5×5 Lattice with 2/3× Material Scaling (t = $(round(data_2_3["time"][div(length(data_2_3["time"]), 2)], digits=2)) s)",
                 joinpath(figures_dir, "material_scaling_2_3_snapshot.png"))

# For 10×10, take snapshot shortly after force is applied (around 0.2-0.25s)
# Find the frame index closest to 0.2s after force application
target_time_10x10 = 0.2  # Shortly after F_ACTIVE_TIME (0.15s)
frame_idx_10x10 = 1
for i in 1:length(data_10x10["time"])
    if data_10x10["time"][i] >= target_time_10x10
        frame_idx_10x10 = i
        break
    end
end
generate_snapshot(data_10x10["solution"], data_10x10["equilibrium_grid"], backplate_positions_10x10,
                 "10×10 Lattice with 1.5× Material Scaling (t = $(round(data_10x10["time"][frame_idx_10x10], digits=2)) s)",
                 joinpath(figures_dir, "material_scaling_10x10_snapshot.png"), frame_idx_10x10)
println()

# Print LaTeX table
println("="^70)
println("LaTeX TABLE FOR MATERIAL SCALING RESULTS")
println("="^70)
println()
println("\\begin{table}[h]")
println("\\centering")
println("\\begin{tabular}{|l|c|c|c|}")
println("\\hline")
println("\\textbf{Configuration} & \\textbf{Work (J)} & \\textbf{Final Energy (J)} & \\textbf{Dissipation (\\%)} \\\\")
println("\\hline")
@printf("5×5, 1.5× multiplier & %.6f & %.6f & %.2f \\\\\n",
        data_1_5["total_work"], data_1_5["final_energy"], data_1_5["dissipation_percent"])
@printf("5×5, 2/3× multiplier & %.6f & %.6f & %.2f \\\\\n",
        data_2_3["total_work"], data_2_3["final_energy"], data_2_3["dissipation_percent"])
@printf("10×10, 1.5× multiplier & %.6f & %.6f & %.2f \\\\\n",
        data_10x10["total_work"], data_10x10["final_energy"], data_10x10["dissipation_percent"])
println("\\hline")
println("\\end{tabular}")
println("\\caption{Material scaling comparison results.}")
println("\\label{tab:material_scaling}")
println("\\end{table}")
println()

println("="^70)
println("All materials generated successfully!")
println("="^70)

