#!/usr/bin/env julia
"""
Compute baseline peak force for sequential material ordering [1, 2, 3, ..., 11]
This provides a comparison point for the optimization results.
"""

using Printf
using Dates

# Load simulation
include(joinpath(@__DIR__, "..", "src", "lattice_simulation_11x11.jl"))

println("="^80)
println("Baseline Sequential Material Ordering Evaluation")
println("="^80)
println()

# Sequential ordering: [1, 2, 3, ..., 11]
baseline_order = collect(1:11)
println("Material Ordering: $baseline_order")
println("(Sequential: stiffest to softest)")
println()

# Run simulation
println("Running simulation...")
start_time = time()

peak_force = run_simulation_with_material_ordering(baseline_order)

elapsed_time = time() - start_time

println()
println("="^80)
println("Baseline Results")
println("="^80)
@printf("Material Ordering: %s\n", baseline_order)
@printf("Peak Backplate Force: %.2f N\n", peak_force)
@printf("Simulation Time: %.2f seconds (%.2f minutes)\n", elapsed_time, elapsed_time/60)
println("="^80)

# Save to file
output_file = joinpath(@__DIR__, "optimization", "comparison_output", "baseline_result.txt")
open(output_file, "w") do f
    println(f, "# Baseline Sequential Material Ordering Result")
    println(f, "# Date: $(now())")
    println(f, "")
    println(f, "Material Ordering: $baseline_order")
    println(f, "Peak Backplate Force (N): $peak_force")
    println(f, "Simulation Time (s): $elapsed_time")
end

@printf("\nResults saved to: %s\n", output_file)

