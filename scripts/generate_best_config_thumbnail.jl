#!/usr/bin/env julia
# Save papers/figures/best_optimized_configuration_thumbnail.png (default: overnight Optim best ordering).
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using GLMakie

src_dir = joinpath(@__DIR__, "..", "src")
figures_dir = joinpath(@__DIR__, "..", "papers", "figures")
mkpath(figures_dir)
include(joinpath(src_dir, "lattice_simulation_11x11.jl"))

material_order = [11, 10, 9, 8, 7, 6, 4, 3, 2, 1, 5]
equilibrium_grid = create_equilibrium_grid()
backplate_positions = create_backplate_positions()
current_positions = equilibrium_grid

nn_connections = Tuple{Int, Int}[]
for i in 1:N, j in 1:N - 1
    push!(nn_connections, (lattice_idx(i, j), lattice_idx(i, j + 1)))
end
for i in 1:N - 1, j in 1:N
    push!(nn_connections, (lattice_idx(i, j), lattice_idx(i + 1, j)))
end

fig = Figure(size = (800, 800), fontsize = 14)
ax = Axis(fig[1, 1];
          title = "Best-found material ordering (initial geometry, t = 0)",
          xlabel = "X Position (m)",
          ylabel = "Y Position (m)",
          aspect = DataAspect())

for (k1, k2) in nn_connections
    pts = Point2f[current_positions[:, k1], current_positions[:, k2]]
    lines!(ax, pts; color = :steelblue, linewidth = 2)
end

bp_line = Point2f[backplate_positions[:, i] for i in 1:N]
lines!(ax, bp_line; color = :gray, linewidth = 8)

mass_points = Point2f[current_positions[:, k] for k in 1:TOTAL_MASSES]
mat_ids = Float64[material_order[lattice_j(k)] for k in 1:TOTAL_MASSES]
scatter!(ax, mass_points; markersize = 12, color = mat_ids, colormap = :tab10,
         colorrange = (1, 11))

# Highlight a representative loaded node (column 1, row 6), consistent with analysis plots.
load_idx = lattice_idx(6, 1)
scatter!(ax, [Point2f(current_positions[:, load_idx])]; markersize = 18,
         color = (:orange, 0.9), strokewidth = 2, strokecolor = :black)

grid_min = minimum(equilibrium_grid) - 0.5
grid_max_x = maximum(backplate_positions[1, :]) + 0.5
grid_max_y = maximum(backplate_positions[2, :]) + 0.5
limits!(ax, grid_min, grid_max_x, grid_min, grid_max_y)

out = joinpath(figures_dir, "best_optimized_configuration_thumbnail.png")
save(out, fig)
println("Saved: ", out)
