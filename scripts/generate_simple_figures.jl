using Pkg
Pkg.activate(".")

using CairoMakie
CairoMakie.activate!()

# Paths
figures_dir = joinpath(@__DIR__, "..", "papers", "figures")
mkpath(figures_dir)

# Figure 1: Two-Mass Spring System
function generate_two_mass_figure()
    fig = Figure(size = (800, 600), fontsize = 14)
    ax = Axis(fig[1, 1],
              title = "Two-Mass Spring System",
              xlabel = "X Position (m)",
              ylabel = "Y Position (m)",
              aspect = DataAspect(),
              limits = (-0.5, 2.5, -0.5, 0.5))
    
    # Mass positions
    m1_pos = Point2f(0.0, 0.0)
    m2_pos = Point2f(1.0, 0.0)
    
    # Draw spring (wavy line)
    spring_points = Point2f[]
    n_spring_points = 20
    for i in 0:n_spring_points
        x = i / n_spring_points
        y = 0.05 * sin(2π * 3 * x)  # Wavy spring
        push!(spring_points, Point2f(x, y))
    end
    lines!(ax, spring_points, color = :blue, linewidth = 3, label = "Spring")
    
    # Draw masses
    scatter!(ax, [m1_pos], markersize = 30, color = :red, label = "Mass 1 (m₁)")
    scatter!(ax, [m2_pos], markersize = 30, color = :red, label = "Mass 2 (m₂)")
    
    # Add labels
    text!(ax, "m₁", position = Point2f(-0.2, 0.0), align = (:center, :center), fontsize = 16)
    text!(ax, "m₂", position = Point2f(1.2, 0.0), align = (:center, :center), fontsize = 16)
    text!(ax, "k, α", position = Point2f(0.5, 0.15), align = (:center, :center), fontsize = 12, color = :blue)
    
    # Add coordinate axes (simple lines with arrowheads)
    lines!(ax, Point2f[Point2f(-0.3, 0.0), Point2f(0.2, 0.0)], color = :black, linewidth = 1)
    lines!(ax, Point2f[Point2f(0.0, -0.3), Point2f(0.0, 0.2)], color = :black, linewidth = 1)
    scatter!(ax, [Point2f(0.2, 0.0)], markersize = 8, color = :black, marker = :utriangle)
    scatter!(ax, [Point2f(0.0, 0.2)], markersize = 8, color = :black, marker = :rtriangle)
    text!(ax, "x", position = Point2f(0.25, -0.15), align = (:center, :center), fontsize = 12)
    text!(ax, "y", position = Point2f(-0.15, 0.25), align = (:center, :center), fontsize = 12)
    
    axislegend(ax, position = :rt)
    
    save(joinpath(figures_dir, "two_mass_spring_system.png"), fig)
    println("Saved: two_mass_spring_system.png")
end

# Figure 2: 5x5 Lattice System
function generate_5x5_lattice_figure()
    fig = Figure(size = (800, 800), fontsize = 14)
    ax = Axis(fig[1, 1],
              title = "5×5 Mass-Spring Lattice System",
              xlabel = "X Position (m)",
              ylabel = "Y Position (m)",
              aspect = DataAspect())
    
    N = 5
    grid_spacing = 1.0
    
    # Create grid positions
    positions = Point2f[]
    for i in 1:N
        for j in 1:N
            push!(positions, Point2f((j-1)*grid_spacing, (i-1)*grid_spacing))
        end
    end
    
    # Function to convert (i,j) to index
    idx(i, j) = (i-1)*N + j
    
    # Draw nearest neighbor springs (horizontal and vertical)
    for i in 1:N
        for j in 1:N
            # Horizontal connections
            if j < N
                k1 = idx(i, j)
                k2 = idx(i, j+1)
                lines!(ax, Point2f[positions[k1], positions[k2]], color = :blue, linewidth = 2)
            end
            # Vertical connections
            if i < N
                k1 = idx(i, j)
                k2 = idx(i+1, j)
                lines!(ax, Point2f[positions[k1], positions[k2]], color = :blue, linewidth = 2)
            end
        end
    end
    
    # Draw diagonal springs
    for i in 1:N
        for j in 1:N
            # Diagonal down-right
            if i < N && j < N
                k1 = idx(i, j)
                k2 = idx(i+1, j+1)
                lines!(ax, Point2f[positions[k1], positions[k2]], color = :red, linewidth = 1.5, linestyle = :dash)
            end
            # Diagonal down-left
            if i < N && j > 1
                k1 = idx(i, j)
                k2 = idx(i+1, j-1)
                lines!(ax, Point2f[positions[k1], positions[k2]], color = :red, linewidth = 1.5, linestyle = :dash)
            end
        end
    end
    
    # Draw masses
    scatter!(ax, positions, markersize = 20, color = :black, strokewidth = 2, strokecolor = :white)
    
    # Add legend
    lines!(ax, Point2f[Point2f(-0.5, -0.5), Point2f(-0.3, -0.5)], color = :blue, linewidth = 2, label = "Nearest Neighbor")
    lines!(ax, Point2f[Point2f(-0.5, -0.6), Point2f(-0.3, -0.6)], color = :red, linewidth = 1.5, linestyle = :dash, label = "Diagonal")
    axislegend(ax, position = :rt)
    
    save(joinpath(figures_dir, "5x5_lattice_system.png"), fig)
    println("Saved: 5x5_lattice_system.png")
end

# Generate both figures
println("Generating simple system figures...")
generate_two_mass_figure()
generate_5x5_lattice_figure()
println("Done!")
