"""
Script to generate analysis plots for the research paper:
- Position vs time for representative mass points
- Velocity vs time for representative mass points
- Kinetic and potential energy vs time
- Total energy vs time
- Work done by external forces vs time
- Backplate reaction force vs time

Usage: Run this script to generate all required analysis figures.
"""

using Pkg
Pkg.activate(".")

using GLMakie
# Use GLMakie for plotting (same as other scripts)
using Printf
using DifferentialEquations
using LinearAlgebra

# Paths
src_dir = joinpath(@__DIR__, "..", "src")
figures_dir = joinpath(@__DIR__, "..", "papers", "figures")
mkpath(figures_dir)

# Include the simulation code
include(joinpath(src_dir, "lattice_simulation_11x11.jl"))

println("="^70)
println("Generating Analysis Plots")
println("="^70)
println()

# Run simulation with default material ordering
println("Running simulation...")
material_order = collect(1:N)  # Sequential ordering
u0 = zeros(2 * TOTAL_DOF)
tspan = (0.0, T_END)
p = (material_order=material_order, materials=DEFAULT_MATERIALS)
prob = ODEProblem(lattice_2d_rhs_with_diagonals_and_backplate!, u0, tspan, p)
sol = solve(prob, Vern9();
            reltol = REL_TOL, 
            abstol = ABS_TOL,
            saveat = OUTPUT_INTERVAL,
            dense = false)

if sol.retcode != :Success
    error("Solver failed: $(sol.retcode)")
end

println("Simulation completed successfully!")
println("Time steps: $(length(sol.t))")
println()

# Create equilibrium grid
equilibrium_grid = create_equilibrium_grid()

# Helper function to get mass index from (i, j)
function get_mass_idx(i, j)
    return lattice_idx(i, j)
end

# Representative mass points
mass_points = [
    (1, 1, "Corner (1,1)"),
    (6, 1, "Left Edge (6,1)"),
    (6, 6, "Center (6,6)"),
    (6, 11, "Right Edge (6,11)")
]

# Extract position and velocity data
println("Extracting position and velocity data...")
positions_x = Dict()
positions_y = Dict()
velocities_x = Dict()
velocities_y = Dict()
velocity_magnitudes = Dict()

for (i, j, label) in mass_points
    k = get_mass_idx(i, j)
    pos_x = Float64[]
    pos_y = Float64[]
    vel_x = Float64[]
    vel_y = Float64[]
    vel_mag = Float64[]
    
    for n in 1:length(sol.t)
        pos = reshape(view(sol.u[n], 1:TOTAL_DOF), 2, TOTAL_MASSES)
        vel = reshape(view(sol.u[n], TOTAL_DOF+1:2*TOTAL_DOF), 2, TOTAL_MASSES)
        
        # Get equilibrium position
        eq_pos = equilibrium_grid[:, k]
        current_pos = pos[:, k] .+ eq_pos
        current_vel = vel[:, k]
        
        push!(pos_x, current_pos[1])
        push!(pos_y, current_pos[2])
        push!(vel_x, current_vel[1])
        push!(vel_y, current_vel[2])
        push!(vel_mag, norm(current_vel))
    end
    
    positions_x[label] = pos_x
    positions_y[label] = pos_y
    velocities_x[label] = vel_x
    velocities_y[label] = vel_y
    velocity_magnitudes[label] = vel_mag
end

# Generate position vs time plot
println("Generating position vs time plot...")
fig_pos = Figure(size = (1200, 800), fontsize = 14)
ax_x = Axis(fig_pos[1, 1], title = "X Position vs Time", xlabel = "Time (s)", ylabel = "X Position (m)")
ax_y = Axis(fig_pos[2, 1], title = "Y Position vs Time", xlabel = "Time (s)", ylabel = "Y Position (m)")

colors = [:blue, :red, :green, :orange]
for (idx, (i, j, label)) in enumerate(mass_points)
    lines!(ax_x, sol.t, positions_x[label], label = label, color = colors[idx], linewidth = 2)
    lines!(ax_y, sol.t, positions_y[label], label = label, color = colors[idx], linewidth = 2)
end

axislegend(ax_x, position = :rt)
axislegend(ax_y, position = :rt)

save(joinpath(figures_dir, "position_vs_time.png"), fig_pos)
println("  Saved: position_vs_time.png")

# Generate velocity vs time plot
println("Generating velocity vs time plot...")
fig_vel = Figure(size = (1200, 800), fontsize = 14)
ax_vx = Axis(fig_vel[1, 1], title = "X Velocity vs Time", xlabel = "Time (s)", ylabel = "X Velocity (m/s)")
ax_vy = Axis(fig_vel[2, 1], title = "Y Velocity vs Time", xlabel = "Time (s)", ylabel = "Y Velocity (m/s)")

for (idx, (i, j, label)) in enumerate(mass_points)
    lines!(ax_vx, sol.t, velocities_x[label], label = label, color = colors[idx], linewidth = 2)
    lines!(ax_vy, sol.t, velocities_y[label], label = label, color = colors[idx], linewidth = 2)
end

axislegend(ax_vx, position = :rt)
axislegend(ax_vy, position = :rt)

save(joinpath(figures_dir, "velocity_vs_time.png"), fig_vel)
println("  Saved: velocity_vs_time.png")

# Extract energy data
println("Extracting energy data...")
kinetic_energies = Float64[]
potential_energies = Float64[]
total_energies = Float64[]

for n in 1:length(sol.t)
    pos = reshape(view(sol.u[n], 1:TOTAL_DOF), 2, TOTAL_MASSES)
    vel = reshape(view(sol.u[n], TOTAL_DOF+1:2*TOTAL_DOF), 2, TOTAL_MASSES)
    
    ke = kinetic_energy_2d(vel)
    pe = potential_energy_2d_with_diagonals_and_backplate(pos)
    
    push!(kinetic_energies, ke)
    push!(potential_energies, pe)
    push!(total_energies, ke + pe)
end

# Calculate work done
println("Calculating work done by external forces...")
total_work = work_done_2d_distributed(sol)
work_values = Float64[]

# Calculate cumulative work at each time step
fx_unit, fy_unit = calculate_force_components(1.0, FORCE_ANGLE_DEGREES)
let cumulative_work = 0.0
    for n in 1:length(sol.t)
    if n == 1
        push!(work_values, 0.0)
        continue
    end
    
    if sol.t[n] <= F_ACTIVE_TIME
        # Calculate work increment
        dt = sol.t[n] - sol.t[n-1]
        for i in 1:N
            force_mag = get_distributed_force_magnitude(i)
            target_idx = lattice_idx(i, 1)
            
            # Get velocities at midpoint
            vel_prev = reshape(view(sol.u[n-1], TOTAL_DOF+1:2*TOTAL_DOF), 2, TOTAL_MASSES)[:, target_idx]
            vel_curr = reshape(view(sol.u[n], TOTAL_DOF+1:2*TOTAL_DOF), 2, TOTAL_MASSES)[:, target_idx]
            vel_avg = 0.5 * (vel_prev + vel_curr)
            
            # Power = F · v
            power = force_mag * (fx_unit * vel_avg[1] + fy_unit * vel_avg[2])
            cumulative_work += power * dt
        end
    end
    
        push!(work_values, cumulative_work)
    end
end

# Generate energy plots
println("Generating energy plots...")
fig_energy = Figure(size = (1200, 600), fontsize = 14)
ax_energy = Axis(fig_energy[1, 1], title = "Energy vs Time", xlabel = "Time (s)", ylabel = "Energy (J)")

lines!(ax_energy, sol.t, kinetic_energies, label = "Kinetic Energy", color = :blue, linewidth = 2)
lines!(ax_energy, sol.t, potential_energies, label = "Potential Energy", color = :red, linewidth = 2)
lines!(ax_energy, sol.t, total_energies, label = "Total Energy", color = :green, linewidth = 2)
lines!(ax_energy, sol.t, work_values, label = "Work Input", color = :orange, linestyle = :dash, linewidth = 2)

axislegend(ax_energy, position = :rt)

save(joinpath(figures_dir, "energy_vs_time.png"), fig_energy)
println("  Saved: energy_vs_time.png")

# Calculate backplate force
println("Calculating backplate reaction force...")
backplate_forces = Float64[]

equilibrium_x = (N - 1) * GRID_SPACING
backplate_x = equilibrium_x + BACKPLATE_DISTANCE

for n in 1:length(sol.t)
    pos = reshape(view(sol.u[n], 1:TOTAL_DOF), 2, TOTAL_MASSES)
    vel = reshape(view(sol.u[n], TOTAL_DOF+1:2*TOTAL_DOF), 2, TOTAL_MASSES)
    
    total_force = 0.0
    
    for i in 1:N
        k = lattice_idx(i, N)  # Rightmost column
        pos_k = pos[:, k]
        vel_k = vel[:, k]
        
        current_x = equilibrium_x + pos_k[1]
        
        if current_x > backplate_x
            penetration = current_x - backplate_x
            wall_force = WALL_STIFFNESS * penetration
            
            if vel_k[1] > 0
                wall_force += WALL_DAMPING * vel_k[1]
            end
            
            total_force += abs(wall_force)
        end
    end
    
    push!(backplate_forces, total_force)
end

# Generate backplate force plot
println("Generating backplate force plot...")
fig_force = Figure(size = (1000, 600), fontsize = 14)
ax_force = Axis(fig_force[1, 1], title = "Backplate Reaction Force vs Time", 
                xlabel = "Time (s)", ylabel = "Force Magnitude (N)")

lines!(ax_force, sol.t, backplate_forces, color = :purple, linewidth = 2)

save(joinpath(figures_dir, "backplate_force_vs_time.png"), fig_force)
println("  Saved: backplate_force_vs_time.png")

# Generate work plot separately
println("Generating work vs time plot...")
fig_work = Figure(size = (1000, 600), fontsize = 14)
ax_work = Axis(fig_work[1, 1], title = "Cumulative Work Done by External Forces vs Time",
               xlabel = "Time (s)", ylabel = "Work (J)")

lines!(ax_work, sol.t, work_values, color = :orange, linewidth = 2)

save(joinpath(figures_dir, "work_vs_time.png"), fig_work)
println("  Saved: work_vs_time.png")

println()
println("="^70)
println("All plots generated successfully!")
println("="^70)
println("Figures saved to: $figures_dir")
