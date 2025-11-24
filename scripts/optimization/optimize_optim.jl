"""
Optim.jl Optimization Script for Material Ordering

This script uses Optim.jl to find the optimal material ordering
that minimizes peak force transferred to the backplate.

For armor applications (e.g., Kevlar vest simulation), peak force is the
critical metric as it determines injury thresholds and penetration risk.
"""

using Optim
using Printf

# Load simulation and utilities
include(joinpath(@__DIR__, "..", "..", "src", "lattice_simulation_11x11.jl"))
include(joinpath(@__DIR__, "optimization_utils.jl"))


########################################################################
#  OPTIMIZATION CONFIGURATION
########################################################################
# Use regular variables to avoid conflicts when multiple scripts are included
if !@isdefined(MAX_ITERATIONS)
    const MAX_ITERATIONS = 100  # Maximum number of iterations
end
if !@isdefined(VERBOSE)
    const VERBOSE = true
end
if !@isdefined(TRACK_HISTORY)
    const TRACK_HISTORY = true
end


########################################################################
#  MAIN OPTIMIZATION FUNCTION
########################################################################
function optimize_material_ordering_optim(; 
    materials=DEFAULT_MATERIALS,
    max_iterations=MAX_ITERATIONS,
    verbose=VERBOSE,
    track_history=TRACK_HISTORY
)
    """
    Optimize material ordering using Optim.jl.
    
    Note: Optim.jl is primarily for continuous optimization, so we use
    SimulatedAnnealing which can handle discrete problems better.
    
    Args:
        materials: Array of material property tuples (default: DEFAULT_MATERIALS)
        max_iterations: Maximum number of iterations
        verbose: If true, print progress
        track_history: If true, track evaluation history
    
    Returns:
        Tuple (best_permutation, best_force, state, result)
    """
    println("="^80)
    println("Optim.jl Material Ordering Optimization")
    println("="^80)
    println("Max Iterations: $max_iterations")
    println("Materials: $(length(materials)) distinct materials")
    println()
    
    # Create objective function
    obj_func, state = create_objective_function(
        run_simulation_with_material_ordering,
        N,
        DEFAULT_MATERIALS;
        materials=materials,
        verbose=verbose,
        track_history=track_history
    )
    
    # For Optim.jl, we'll use a continuous encoding and convert to permutation
    # Start with random permutation converted to continuous
    initial_perm = random_permutation(N)
    initial_x = permutation_to_continuous(initial_perm, N)
    
    # Define continuous objective function
    function continuous_objective(x)
        # Convert continuous to permutation
        perm = continuous_to_permutation(x, N)
        return obj_func(perm)
    end
    
    # Use SimulatedAnnealing for better handling of discrete-like problems
    # or use a gradient-free method
    println("Starting optimization with SimulatedAnnealing...")
    start_time = time()
    
    # Try SimulatedAnnealing (if available) or use Nelder-Mead
    try
        result = Optim.optimize(
            continuous_objective,
            initial_x,
            Optim.SimulatedAnnealing(),
            Optim.Options(
                iterations = max_iterations,
                show_trace = verbose,
                store_trace = track_history
            )
        )
    catch
        # Fallback to Nelder-Mead if SimulatedAnnealing not available
        println("SimulatedAnnealing not available, using Nelder-Mead...")
        result = Optim.optimize(
            continuous_objective,
            initial_x,
            Optim.NelderMead(),
            Optim.Options(
                iterations = max_iterations,
                show_trace = verbose,
                store_trace = track_history
            )
        )
    end
    
    elapsed_time = time() - start_time
    
    # Extract best solution
    best_continuous = Optim.minimizer(result)
    best_permutation = continuous_to_permutation(best_continuous, N)
    best_force = Optim.minimum(result)
    
    # Update state with final results
    if best_force < state.best_force
        state.best_force = best_force
        state.best_permutation = best_permutation
    end
    
    # Print summary
    println()
    print_optimization_summary(state, "Optim.jl")
    println("Optimization Time: $(round(elapsed_time, digits=2)) seconds")
    
    return best_permutation, best_force, state, result
end


########################################################################
#  RUN OPTIMIZATION
########################################################################
if abspath(PROGRAM_FILE) == @__FILE__
    # Run optimization
    best_perm, best_f, state, result = optimize_material_ordering_optim()
    
    # Save results
    config = Dict(
        "optimizer" => "Optim",
        "method" => "SimulatedAnnealing/NelderMead",
        "max_iterations" => MAX_ITERATIONS
    )
    
    result_file = joinpath(@__DIR__, "results_optim.txt")
    save_optimization_result(result_file, state, "Optim.jl", config)
    
    println("\nOptimization complete! Results saved to $result_file")
end

