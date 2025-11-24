"""
Optimization Utilities for Material Ordering

This module provides shared utilities for optimizing material ordering to minimize
backplate force. It includes objective function wrappers, permutation validation,
and result tracking utilities.
"""

using Printf
using Random
using Dates
using Base: time

# Import the simulation functions
# Note: This assumes the simulation file is already loaded or will be loaded
# The optimization scripts should include both this file and the simulation file

# These will be set when the simulation is loaded
const N_OPT = 11  # Will be overridden by actual N from simulation
const DEFAULT_MATERIALS_OPT = nothing  # Will be set by simulation


########################################################################
#  PERMUTATION VALIDATION AND UTILITIES
########################################################################
function is_valid_permutation(perm, n_val=N_OPT)
    """
    Check if an array is a valid permutation of [1, 2, ..., N].
    
    Args:
        perm: Array to check
        n_val: Value of N (default: N_OPT)
    
    Returns:
        true if valid permutation, false otherwise
    """
    if length(perm) != n_val
        return false
    end
    
    # Check that all values are in range [1, N]
    if !all(1 .<= perm .<= n_val)
        return false
    end
    
    # Check that all values are unique
    if length(Set(perm)) != n_val
        return false
    end
    
    return true
end


function validate_permutation(perm, n_val=N_OPT)
    """
    Validate permutation and throw error if invalid.
    
    Args:
        perm: Permutation array to validate
        n_val: Value of N (default: N_OPT)
    
    Throws:
        ArgumentError if permutation is invalid
    """
    if !is_valid_permutation(perm, n_val)
        error("Invalid permutation: $perm. Must be a permutation of [1, 2, ..., $n_val]")
    end
end


function random_permutation(n_val=N_OPT)
    """
    Generate a random permutation of [1, 2, ..., N].
    
    Args:
        n_val: Value of N (default: N_OPT)
    
    Returns:
        Random permutation array
    """
    return shuffle(collect(1:n_val))
end


########################################################################
#  OBJECTIVE FUNCTION WRAPPER
########################################################################
mutable struct OptimizationState
    """Track optimization progress"""
    iteration::Int
    best_force::Float64
    best_permutation::Vector{Int}
    evaluation_count::Int
    start_time::Float64
    force_history::Vector{Float64}
    permutation_history::Vector{Vector{Int}}
    
    function OptimizationState()
        new(0, Inf, Int[], 0, time(), Float64[], Vector{Int}[])
    end
end


function create_objective_function(run_sim_func, n_val, default_materials; materials=nothing, verbose=false, track_history=false)
    """
    Create an objective function for optimization.
    
    Args:
        run_sim_func: Function to run simulation (e.g., run_simulation_with_material_ordering)
        n_val: Value of N (number of columns)
        default_materials: Default materials array
        materials: Array of material property tuples (default: default_materials)
        verbose: If true, print progress (default: false)
        track_history: If true, track evaluation history (default: false)
    
    Returns:
        Tuple (objective_function, state) where state tracks optimization progress
    """
    state = OptimizationState()
    mat_props = materials === nothing ? default_materials : materials
    
    function objective_function(perm)
        """
        Objective function: minimize peak backplate force.
        
        Args:
            perm: Permutation array [1, 2, ..., N] or continuous encoding
        
        Returns:
            Peak backplate force (maximum force magnitude, to be minimized)
        """
        state.evaluation_count += 1
        
        # Handle continuous encodings (some optimizers use continuous variables)
        # If perm is continuous, convert to permutation
        if eltype(perm) <: AbstractFloat
            perm = continuous_to_permutation(perm, n_val)
        end
        
        # Validate permutation
        try
            validate_permutation(perm, n_val)
        catch e
            # Return large penalty for invalid permutations
            return 1e10
        end
        
        # Run simulation
        try
            force = run_sim_func(perm; materials=mat_props, return_solution=false)
            
            # Update best if better
            if force < state.best_force
                state.best_force = force
                state.best_permutation = copy(perm)
                state.iteration = state.evaluation_count
                
                if verbose
                    elapsed = time() - state.start_time
                    @printf("Iteration %d: Force = %.6f (best so far)\n", 
                           state.evaluation_count, force)
                    @printf("  Permutation: %s\n", perm)
                    @printf("  Time: %.2f s\n", elapsed)
                end
            end
            
            # Track history if requested
            if track_history
                push!(state.force_history, force)
                push!(state.permutation_history, copy(perm))
            end
            
            return force
            
        catch e
            # Return large penalty for failed simulations
            if verbose
                @printf("Warning: Simulation failed for permutation %s: %s\n", perm, e)
            end
            return 1e10
        end
    end
    
    return objective_function, state
end


########################################################################
#  CONTINUOUS ENCODING FOR PERMUTATIONS
########################################################################
function continuous_to_permutation(continuous, n_val=N_OPT)
    """
    Convert continuous encoding to permutation using argsort.
    
    Some optimizers work with continuous variables. We can encode a permutation
    as a continuous array and convert it using argsort.
    
    Args:
        continuous: Array of continuous values
        n_val: Value of N (default: N_OPT)
    
    Returns:
        Permutation array
    """
    if length(continuous) != n_val
        error("Continuous encoding must have length $n_val, got $(length(continuous))")
    end
    
    # Use argsort to convert continuous values to permutation
    # This gives indices that would sort the array
    perm = sortperm(continuous)
    
    return perm
end


function permutation_to_continuous(perm, n_val=N_OPT)
    """
    Convert permutation to continuous encoding.
    
    Args:
        perm: Permutation array
        n_val: Value of N (default: N_OPT)
    
    Returns:
        Continuous encoding array
    """
    validate_permutation(perm, n_val)
    
    # Create continuous encoding: use permutation indices as values
    continuous = zeros(Float64, n_val)
    for (i, p) in enumerate(perm)
        continuous[p] = Float64(i)
    end
    
    return continuous
end


########################################################################
#  RESULT SAVING AND LOADING
########################################################################
function save_optimization_result(filename, state, optimizer_name, config)
    """
    Save optimization results to a file.
    
    Args:
        filename: Output filename
        state: OptimizationState object
        optimizer_name: Name of optimizer used
        config: Configuration dictionary
    """
    open(filename, "w") do f
        println(f, "# Optimization Results")
        println(f, "# Optimizer: $optimizer_name")
        println(f, "# Date: $(now())")
        println(f, "# Objective: Minimize Peak Force")
        println(f, "")
        println(f, "Best Peak Force: $(state.best_force)")
        println(f, "Best Permutation: $(state.best_permutation)")
        println(f, "Iterations: $(state.evaluation_count)")
        println(f, "Elapsed Time: $(time() - state.start_time) seconds")
        println(f, "")
        println(f, "# Configuration:")
        for (key, value) in config
            println(f, "$key: $value")
        end
        println(f, "")
        println(f, "# Peak Force History:")
        for (i, force) in enumerate(state.force_history)
            println(f, "$i\t$force")
        end
    end
    @printf("Results saved to %s\n", filename)
end


########################################################################
#  HELPER FUNCTIONS
########################################################################
function print_optimization_summary(state, optimizer_name)
    """
    Print a summary of optimization results.
    
    Args:
        state: OptimizationState object
        optimizer_name: Name of optimizer used
    """
    println("\n" * "="^80)
    println("Optimization Summary: $optimizer_name")
    println("="^80)
    println("Best Peak Force: $(state.best_force)")
    println("Best Permutation: $(state.best_permutation)")
    println("Total Evaluations: $(state.evaluation_count)")
    println("Elapsed Time: $(round(time() - state.start_time, digits=2)) seconds")
    println("="^80)
end

