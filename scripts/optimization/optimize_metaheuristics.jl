"""
Metaheuristics.jl Optimization Script for Material Ordering

This script uses Metaheuristics.jl to find the optimal material ordering
that minimizes peak force transferred to the backplate.

For armor applications (e.g., Kevlar vest simulation), peak force is the
critical metric as it determines injury thresholds and penetration risk.
"""

using Metaheuristics
using Metaheuristics: boxconstraints  # Explicitly import boxconstraints
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
function optimize_material_ordering_metaheuristics(; 
    materials=DEFAULT_MATERIALS,
    max_evaluations=MAX_EVALUATIONS,
    verbose=VERBOSE,
    track_history=TRACK_HISTORY,
    algorithm=:PSO  # Options: :PSO, :DE, :ECA, :ES
)
    """
    Optimize material ordering using Metaheuristics.jl.
    
    Args:
        materials: Array of material property tuples (default: DEFAULT_MATERIALS)
        max_evaluations: Maximum number of function evaluations
        verbose: If true, print progress
        track_history: If true, track evaluation history
        algorithm: Algorithm to use (:PSO, :DE, :ECA, :ES)
    
    Returns:
        Tuple (best_permutation, best_force, state, result)
    """
    println("="^80)
    println("Metaheuristics.jl Material Ordering Optimization")
    println("="^80)
    println("Algorithm: $algorithm")
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
    
    # Metaheuristics.jl works with continuous variables
    # We'll use continuous encoding and convert to permutation
    
    # Bounds: [0, 1] for each dimension
    # Metaheuristics.jl expects bounds as lower and upper bounds
    lb = zeros(N)  # Lower bounds
    ub = ones(N)   # Upper bounds
    
    # Define continuous objective function with hard evaluation limit
    # Metaheuristics.jl's f_calls_limit might not be strictly enforced,
    # so we enforce it ourselves in the objective function
    evaluation_count = Ref(0)
    function continuous_objective(x)
        # Enforce hard evaluation limit
        evaluation_count[] += 1
        if evaluation_count[] > max_evaluations
            # Return a large penalty value to signal we've exceeded the limit
            # This should cause the optimizer to stop
            return 1e10
        end
        
        # Convert continuous to permutation
        perm = continuous_to_permutation(x, N)
        return obj_func(perm)
    end
    
    # Select algorithm
    println("Starting optimization with $algorithm...")
    start_time = time()
    
    # Metaheuristics.jl API: optimize(f, bounds, method)
    # Create bounds using boxconstraints
    bounds = boxconstraints(lb = lb, ub = ub)
    
    # Configure algorithms with appropriate population sizes
    # Population size should be reasonable for the problem dimension
    # But must not exceed max_evaluations (population-based algorithms need at least pop_size evaluations)
    base_pop_size = min(50, max(10, 4 * N))  # Adaptive population size
    
    # For very low max_evaluations, we need special handling
    # DE and ECA require initialization (pop_size) + at least 1 iteration
    # PSO can work with just pop_size evaluations
    # Minimum requirements:
    # - PSO: pop_size (typically 2-4) evaluations
    # - DE/ECA: pop_size + 1 (need initialization + at least one iteration)
    # - ES: pop_size (typically 2-4) evaluations
    if algorithm in [:DE, :ECA]
        # DE and ECA need pop_size for initialization + at least 1 more for iteration
        min_pop_size = 4  # Minimum viable population size
        min_evaluations_needed = min_pop_size + 1  # Need initialization + at least 1 iteration
    else
        # PSO and ES can work with just pop_size evaluations
        min_pop_size = 2
        min_evaluations_needed = min_pop_size
    end
    
    # Check if we can actually run this algorithm with the given constraints
    if max_evaluations < min_evaluations_needed
        error_msg = "Algorithm $algorithm requires at least $min_evaluations_needed evaluations (minimum population size $min_pop_size + iterations), but max_evaluations=$max_evaluations. Increase max_evaluations to at least $min_evaluations_needed or use a different algorithm."
        if verbose
            println("ERROR: $error_msg")
        end
        throw(ErrorException(error_msg))
    end
    
    # Calculate actual population size
    # For DE/ECA, we need to leave room for at least 1 iteration after initialization
    # So pop_size should be at most (max_evaluations - 1) for DE/ECA
    if algorithm in [:DE, :ECA]
        # Reserve at least 1 evaluation for iteration, rest for population
        max_pop_for_evaluations = max(min_pop_size, max_evaluations - 1)
        pop_size = min(base_pop_size, max_pop_for_evaluations)
    else
        # PSO and ES can use all evaluations for population if needed
        pop_size = max(min_pop_size, min(base_pop_size, max_evaluations))
    end
    
    if base_pop_size > max_evaluations && verbose
        @warn "Requested population size ($base_pop_size) exceeds max_evaluations ($max_evaluations). Reducing to $pop_size."
    end
    
    # Create options with max function evaluations
    # Options are passed to the algorithm constructor
    options = Metaheuristics.Options(; f_calls_limit = max_evaluations, seed = 1)
    
    if algorithm == :PSO
        # PSO needs population size N parameter and options
        method = Metaheuristics.PSO(N = pop_size, options = options)
        result = Metaheuristics.optimize(continuous_objective, bounds, method)
    elseif algorithm == :DE
        # DE needs population size N parameter and options
        method = Metaheuristics.DE(N = pop_size, options = options)
        result = Metaheuristics.optimize(continuous_objective, bounds, method)
    elseif algorithm == :ECA
        # ECA needs population size N parameter and options
        method = Metaheuristics.ECA(N = pop_size, options = options)
        result = Metaheuristics.optimize(continuous_objective, bounds, method)
    elseif algorithm == :ES
        # Try CMA_ES which is the Evolution Strategy in Metaheuristics.jl
        # Check if CMA_ES is available
        if isdefined(Metaheuristics, :CMA_ES)
            method = Metaheuristics.CMA_ES(N = pop_size, options = options)
            result = Metaheuristics.optimize(continuous_objective, bounds, method)
        else
            error("ES/CMA_ES algorithm not available in this version of Metaheuristics.jl. Use :PSO, :DE, or :ECA instead.")
        end
    else
        error("Unknown algorithm: $algorithm. Use :PSO, :DE, :ECA, or :ES")
    end
    
    elapsed_time = time() - start_time
    
    # Extract best solution
    best_continuous = Metaheuristics.minimizer(result)
    best_permutation = continuous_to_permutation(best_continuous, N)
    best_force = Metaheuristics.minimum(result)
    
    # Update state with final results
    if best_force < state.best_force
        state.best_force = best_force
        state.best_permutation = best_permutation
    end
    
    # Print summary
    println()
    print_optimization_summary(state, "Metaheuristics.jl ($algorithm)")
    println("Optimization Time: $(round(elapsed_time, digits=2)) seconds")
    
    return best_permutation, best_force, state, result
end


########################################################################
#  RUN OPTIMIZATION
########################################################################
if abspath(PROGRAM_FILE) == @__FILE__
    # Test multiple algorithms (ES/CMA_ES might not be available)
    algorithms = [:PSO, :DE, :ECA, :ES]
    results = Dict()
    
    for alg in algorithms
        println("\n" * "="^80)
        println("Testing algorithm: $alg")
        println("="^80)
        
        try
            best_perm, best_f, state, result = optimize_material_ordering_metaheuristics(algorithm=alg)
            results[alg] = (best_perm, best_f, state, result)
            
            # Save results
            config = Dict(
                "optimizer" => "Metaheuristics",
                "method" => string(alg),
                "max_evaluations" => MAX_EVALUATIONS
            )
            
            result_file = joinpath(@__DIR__, "results_metaheuristics_$(alg).txt")
            save_optimization_result(result_file, state, "Metaheuristics.jl ($alg)", config)
            
        catch e
            println("Error with algorithm $alg: $e")
            results[alg] = nothing
        end
    end
    
    # Print comparison
    println("\n" * "="^80)
    println("Algorithm Comparison")
    println("="^80)
    for (alg, res) in results
        if res !== nothing
            _, best_f, _, _ = res
            println("$alg: Force = $best_f")
        else
            println("$alg: Failed")
        end
    end
    
    println("\nOptimization complete!")
end

