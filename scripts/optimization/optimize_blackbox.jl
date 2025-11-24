"""
BlackBoxOptim.jl Optimization Script for Material Ordering

This script uses BlackBoxOptim.jl to find the optimal material ordering
that minimizes peak force transferred to the backplate.

For armor applications (e.g., Kevlar vest simulation), peak force is the
critical metric as it determines injury thresholds and penetration risk.
"""

# BlackBoxOptim is required
using BlackBoxOptim

using Printf

# Load simulation and utilities
include(joinpath(@__DIR__, "..", "..", "src", "lattice_simulation_11x11.jl"))
include(joinpath(@__DIR__, "optimization_utils.jl"))


########################################################################
#  OPTIMIZATION CONFIGURATION
########################################################################
# Use regular variables to avoid conflicts when multiple scripts are included
if !@isdefined(MAX_EVALUATIONS)
    const MAX_EVALUATIONS = 100  # Maximum number of function evaluations
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
function optimize_material_ordering_blackbox(; 
    materials=DEFAULT_MATERIALS,
    max_evaluations=MAX_EVALUATIONS,
    verbose=VERBOSE,
    track_history=TRACK_HISTORY
)
    """
    Optimize material ordering using BlackBoxOptim.jl.
    
    Args:
        materials: Array of material property tuples (default: DEFAULT_MATERIALS)
        max_evaluations: Maximum number of function evaluations
        verbose: If true, print progress
        track_history: If true, track evaluation history
    
    Returns:
        Tuple (best_permutation, best_force, state, result)
    """
    println("="^80)
    println("BlackBoxOptim.jl Material Ordering Optimization")
    println("="^80)
    println("Max Evaluations: $max_evaluations")
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
    
    # BlackBoxOptim works with continuous variables, so we need to encode permutations
    # We'll use a continuous encoding where we optimize N continuous values
    # and convert to permutation using argsort
    
    # Search range: [0, 1] for each dimension (will be converted to permutation)
    search_range = [(0.0, 1.0) for _ in 1:N]
    
    # Define continuous objective function
    function continuous_objective(x)
        # Convert continuous to permutation
        perm = sortperm(x)
        return obj_func(perm)
    end
    
    # Run optimization
    println("Starting optimization...")
    start_time = time()
    
    result = bboptimize(
        continuous_objective;
        SearchRange = search_range,
        MaxFuncEvals = max_evaluations,
        Method = :adaptive_de_rand_1_bin_radiuslimited,  # Good for mixed problems
        TraceMode = verbose ? :verbose : :silent
    )
    
    elapsed_time = time() - start_time
    
    # Extract best solution
    best_continuous = best_candidate(result)
    best_permutation = sortperm(best_continuous)
    best_force = best_fitness(result)
    
    # Update state with final results
    if best_force < state.best_force
        state.best_force = best_force
        state.best_permutation = best_permutation
    end
    
    # Print summary
    println()
    print_optimization_summary(state, "BlackBoxOptim")
    println("Best Continuous Encoding: $best_continuous")
    println("Optimization Time: $(round(elapsed_time, digits=2)) seconds")
    
    return best_permutation, best_force, state, result
end


########################################################################
#  RUN OPTIMIZATION
########################################################################
if abspath(PROGRAM_FILE) == @__FILE__
    # Run optimization
    best_perm, best_f, state, result = optimize_material_ordering_blackbox()
    
    # Save results
    config = Dict(
        "optimizer" => "BlackBoxOptim",
        "method" => "adaptive_de_rand_1_bin_radiuslimited",
        "max_evaluations" => MAX_EVALUATIONS
    )
    
    result_file = joinpath(@__DIR__, "results_blackbox.txt")
    save_optimization_result(result_file, state, "BlackBoxOptim", config)
    
    println("\nOptimization complete! Results saved to $result_file")
end

