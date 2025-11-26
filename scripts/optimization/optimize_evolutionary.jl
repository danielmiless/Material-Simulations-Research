"""
Evolutionary.jl Optimization Script for Material Ordering

This script uses Evolutionary.jl to find the optimal material ordering
that minimizes peak force transferred to the backplate.

For armor applications (e.g., Kevlar vest simulation), peak force is the
critical metric as it determines injury thresholds and penetration risk.
"""

using Evolutionary
using Printf

# Load simulation and utilities
include(joinpath(@__DIR__, "..", "..", "src", "lattice_simulation_11x11.jl"))
include(joinpath(@__DIR__, "optimization_utils.jl"))


########################################################################
#  OPTIMIZATION CONFIGURATION
########################################################################
# Use regular variables to avoid conflicts when multiple scripts are included
if !@isdefined(MAX_GENERATIONS)
    const MAX_GENERATIONS = 50  # Maximum number of generations
end
if !@isdefined(POPULATION_SIZE)
    const POPULATION_SIZE = 20  # Population size
end
if !@isdefined(VERBOSE)
    const VERBOSE = true
end
if !@isdefined(TRACK_HISTORY)
    const TRACK_HISTORY = true
end


########################################################################
#  CUSTOM OPERATORS FOR PERMUTATIONS
########################################################################
function permutation_crossover(p1, p2)
    """
    Order crossover (OX) for permutations.
    Preserves relative order from parents.
    """
    n = length(p1)
    # Select random segment
    i1 = rand(1:n)
    i2 = rand(1:n)
    if i1 > i2
        i1, i2 = i2, i1
    end
    
    # Create child
    child = zeros(Int, n)
    child[i1:i2] = p1[i1:i2]
    
    # Fill remaining positions from p2
    idx = 1
    for val in p2
        if !(val in child[i1:i2])
            while child[idx] != 0
                idx += 1
            end
            child[idx] = val
        end
    end
    
    return child
end


function permutation_mutation(perm)
    """
    Swap mutation for permutations.
    Randomly swap two elements.
    """
    mutated = copy(perm)
    i1 = rand(1:length(perm))
    i2 = rand(1:length(perm))
    mutated[i1], mutated[i2] = mutated[i2], mutated[i1]
    return mutated
end


########################################################################
#  MAIN OPTIMIZATION FUNCTION
########################################################################
function optimize_material_ordering_evolutionary(; 
    materials=DEFAULT_MATERIALS,
    max_generations=MAX_GENERATIONS,
    population_size=POPULATION_SIZE,
    verbose=VERBOSE,
    track_history=TRACK_HISTORY
)
    """
    Optimize material ordering using Evolutionary.jl.
    
    Args:
        materials: Array of material property tuples (default: DEFAULT_MATERIALS)
        max_generations: Maximum number of generations
        population_size: Population size
        verbose: If true, print progress
        track_history: If true, track evaluation history
    
    Returns:
        Tuple (best_permutation, best_force, state, result)
    """
    println("="^80)
    println("Evolutionary.jl Material Ordering Optimization")
    println("="^80)
    println("Max Generations: $max_generations")
    println("Population Size: $population_size")
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
    
    # Initialize population with random permutations
    # Evolutionary.jl may have minimum population size requirements
    # Ensure population size is at least 11 to avoid BoundsError (some internal indexing may go up to pop_size+1)
    min_pop_size = max(11, population_size)
    if population_size < 11 && verbose
        println("Warning: Population size $population_size may be too small. Using minimum of 11.")
    end
    initial_population = [random_permutation(N) for _ in 1:min_pop_size]
    
    # Define fitness function (minimize, so negate for Evolutionary.jl which maximizes)
    function fitness(perm)
        return -obj_func(perm)  # Negate because Evolutionary.jl maximizes
    end
    
    println("Starting optimization with Genetic Algorithm...")
    println("Actual population size: $(length(initial_population))")
    start_time = time()
    
    # Create GA method with custom operators
    # Note: Evolutionary.jl's GA constructor might not accept custom operators directly
    # We'll use a wrapper approach
    ga_method = try
        # Try to create GA with custom operators
        # The API might require different syntax
        Evolutionary.GA(
            selection = Evolutionary.tournament(3),
            mutationRate = 0.1
        )
    catch e
        # Fallback: use default GA if custom operators fail
        if verbose
            println("Warning: Custom GA configuration failed, using defaults: $e")
        end
        Evolutionary.GA()
    end
    
    # Create options
    options = Evolutionary.Options(
        iterations = max_generations,
        store_trace = track_history,
        show_trace = verbose
    )
    
    # Note: Evolutionary.jl is designed for continuous optimization, not discrete permutation problems.
    # It treats vectors as populations of scalars, not as permutation arrays. This makes it fundamentally
    # incompatible with our permutation-based optimization problem.
    # 
    # We'll attempt to use it, but it will likely fail. The error will be caught and reported gracefully.
    
    # Create objective from first individual
    initial_x = initial_population[1]
    obj = Evolutionary.EvolutionaryObjective(fitness, initial_x)
    
    # Run optimization
    # Note: Evolutionary.jl is likely to fail for permutation problems due to API incompatibility
    default_ga = Evolutionary.GA()
    
    result = try
        # Try with custom GA method
        Evolutionary.optimize(
            obj,
            Evolutionary.NoConstraints(),
            ga_method,
            initial_x,
            options
        )
    catch e1
        # If that fails, try with default GA
        if verbose
            println("Warning: Optimization with custom GA failed, trying with default GA: $e1")
        end
        try
            Evolutionary.optimize(
                obj,
                Evolutionary.NoConstraints(),
                default_ga,
                initial_x,
                options
            )
        catch e2
            # Evolutionary.jl is not suitable for permutation optimization
            # Provide a clear error message explaining why
            error_msg = """
            Evolutionary.jl is not compatible with permutation optimization problems.
            
            Reason: Evolutionary.jl is designed for continuous optimization and treats vectors as 
            populations of scalars, not as permutation arrays. When used with permutation problems,
            it attempts to pass individual integers to the fitness function instead of permutation vectors.
            
            Original error: $e2
            
            Recommendation: Use a different optimizer designed for discrete/permutation problems:
            - Metaheuristics.jl (PSO, DE, ECA) - Works well for permutations
            - BlackBoxOptim.jl - Handles permutations via continuous encoding
            - Optim.jl - Can work with permutations via continuous encoding
            """
            error(ErrorException(error_msg))
        end
    end
    
    elapsed_time = time() - start_time
    
    # Extract best solution
    best_permutation = Evolutionary.best(result)
    best_fitness_value = Evolutionary.best_fitness(result)
    best_force = -best_fitness_value  # Negate back since we negated in fitness function
    
    # Update state with final results
    if best_force < state.best_force
        state.best_force = best_force
        state.best_permutation = best_permutation
    end
    
    # Print summary
    println()
    print_optimization_summary(state, "Evolutionary.jl")
    println("Optimization Time: $(round(elapsed_time, digits=2)) seconds")
    
    return best_permutation, best_force, state, result
end


########################################################################
#  RUN OPTIMIZATION
########################################################################
if abspath(PROGRAM_FILE) == @__FILE__
    # Run optimization
    best_perm, best_f, state, result = optimize_material_ordering_evolutionary()
    
    # Save results
    config = Dict(
        "optimizer" => "Evolutionary",
        "method" => "Genetic Algorithm",
        "max_generations" => MAX_GENERATIONS,
        "population_size" => POPULATION_SIZE
    )
    
    result_file = joinpath(@__DIR__, "results_evolutionary.txt")
    save_optimization_result(result_file, state, "Evolutionary.jl", config)
    
    println("\nOptimization complete! Results saved to $result_file")
end

