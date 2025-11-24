"""
Basic test script to verify the optimization setup works.

This script tests the basic functionality without requiring all packages.
"""

println("="^80)
println("Basic Optimization Setup Test")
println("="^80)

# Test 1: Load simulation
println("\n1. Testing simulation loading...")
try
    include(joinpath(@__DIR__, "..", "..", "src", "lattice_simulation_11x11.jl"))
    println("   ✓ Simulation loaded successfully")
catch e
    println("   ✗ Failed to load simulation: $e")
    exit(1)
end

# Test 2: Test material ordering function
println("\n2. Testing material ordering function...")
try
    test_order = collect(1:11)
    force = run_simulation_with_material_ordering(test_order)
    println("   ✓ Force calculation works: $force")
catch e
    println("   ✗ Force calculation failed: $e")
    exit(1)
end

# Test 3: Test optimization utilities
println("\n3. Testing optimization utilities...")
try
    include(joinpath(@__DIR__, "optimization_utils.jl"))
    println("   ✓ Optimization utilities loaded")
catch e
    println("   ✗ Failed to load utilities: $e")
    exit(1)
end

# Test 4: Check available optimizers
println("\n4. Checking available optimizers...")
available = []

# Check Optim
try
    using Optim
    push!(available, "Optim.jl")
    println("   ✓ Optim.jl available")
catch
    println("   ✗ Optim.jl not available")
end

# Check Evolutionary
try
    using Evolutionary
    push!(available, "Evolutionary.jl")
    println("   ✓ Evolutionary.jl available")
catch
    println("   ✗ Evolutionary.jl not available")
end

# Check Metaheuristics
try
    using Metaheuristics
    push!(available, "Metaheuristics.jl")
    println("   ✓ Metaheuristics.jl available")
catch
    println("   ✗ Metaheuristics.jl not available")
end

# Check BlackBoxOptim
try
    using BlackBoxOptim
    push!(available, "BlackBoxOptim.jl")
    println("   ✓ BlackBoxOptim.jl available")
catch
    println("   ✗ BlackBoxOptim.jl not available (optional)")
end

println("\n" * "="^80)
println("Test Summary")
println("="^80)
println("Available optimizers: $(length(available))")
for opt in available
    println("  - $opt")
end

if length(available) == 0
    println("\n⚠️  No optimizers available. Please install at least one:")
    println("   using Pkg; Pkg.add(\"Optim\")")
    exit(1)
else
    println("\n✓ Basic setup test passed!")
    println("You can now run optimization scripts with available optimizers.")
end

