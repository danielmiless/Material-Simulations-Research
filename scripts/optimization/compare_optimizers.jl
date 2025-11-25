"""
Compare All Optimizers for Material Ordering

This script runs all optimization algorithms and compares their performance
in terms of convergence, final objective value, and computation time.

Generates:
- Convergence plots (force vs evaluations)
- LaTeX comparison tables
- CSV data files
- Summary statistics
"""

using Printf
using DelimitedFiles
using CairoMakie
using Dates
using Statistics
using Base.Threads

# Load simulation and utilities
include(joinpath(@__DIR__, "..", "..", "src", "lattice_simulation_11x11.jl"))
include(joinpath(@__DIR__, "optimization_utils.jl"))

# Load optimizer scripts
include(joinpath(@__DIR__, "optimize_blackbox.jl"))
include(joinpath(@__DIR__, "optimize_optim.jl"))
include(joinpath(@__DIR__, "optimize_evolutionary.jl"))
include(joinpath(@__DIR__, "optimize_metaheuristics.jl"))


########################################################################
#  COMPARISON CONFIGURATION
########################################################################
# Can be overridden via command line or environment variable
const DEFAULT_MAX_EVALUATIONS = 20  # Reduced for faster testing (increase to 50-100 for final runs)
const VERBOSE = true
const TRACK_HISTORY = true

# Check for command line argument, environment variable, or Julia variable
function get_max_evaluations()
    # First check if MAX_EVALUATIONS is already defined as a Julia variable in Main
    # (e.g., set via -e 'MAX_EVALUATIONS=50; include(...)')
    if isdefined(Main, :MAX_EVALUATIONS)
        return Main.MAX_EVALUATIONS
    end
    # Check command line args
    if length(ARGS) > 0
        try
            return parse(Int, ARGS[1])
        catch
            @warn "Could not parse max_evaluations from command line, using default"
        end
    end
    # Check environment variable
    if haskey(ENV, "MAX_EVALUATIONS")
        try
            return parse(Int, ENV["MAX_EVALUATIONS"])
        catch
            @warn "Could not parse MAX_EVALUATIONS from environment, using default"
        end
    end
    return DEFAULT_MAX_EVALUATIONS
end


########################################################################
#  MAIN COMPARISON FUNCTION
########################################################################
function compare_all_optimizers(;
    materials=DEFAULT_MATERIALS,
    max_evaluations=nothing,
    verbose=VERBOSE,
    track_history=TRACK_HISTORY
)
    """
    Run all optimizers and compare results.
    
    Args:
        materials: Array of material property tuples
        max_evaluations: Maximum number of function evaluations (default: from get_max_evaluations())
        verbose: If true, print progress
        track_history: If true, track evaluation history
    
    Returns:
        Dictionary of results for each optimizer
    """
    # Use provided value or get from config
    if max_evaluations === nothing
        max_evaluations = get_max_evaluations()
    end
    
    println("="^80)
    println("Optimizer Comparison: Material Ordering Optimization")
    println("="^80)
    println("Max Evaluations: $max_evaluations")
    println("Materials: $(length(materials)) distinct materials")
    println()
    
    results = Dict()
    overall_start_time = time()
    total_optimizers = 1 + 2 + 4  # BlackBoxOptim + Optim + Evolutionary + 4 Metaheuristics
    
    println("Running all optimizers in parallel for faster execution...")
    println("  All optimizers will run concurrently")
    println("  Total time ≈ longest optimizer (not sum of all)")
    println()
    
    # Thread-safe progress tracking
    completed_optimizers = Threads.Atomic{Int}(0)
    started_optimizers = Threads.Atomic{Int}(0)
    progress_lock = ReentrantLock()
    running_optimizers = Set{String}()  # Track which optimizers are currently running
    
    # Progress tracking function (thread-safe)
    function update_progress(optimizer_name, status, elapsed=nothing, increment_completed=false, increment_started=false)
        lock(progress_lock)
        try
            if increment_started
                Threads.atomic_add!(started_optimizers, 1)
                push!(running_optimizers, optimizer_name)
            end
            if increment_completed
                Threads.atomic_add!(completed_optimizers, 1)
                delete!(running_optimizers, optimizer_name)
            end
            
            started = started_optimizers[]
            completed = completed_optimizers[]
            progress_pct = round(100 * completed / total_optimizers, digits=1)
            timestamp = Dates.format(now(), "HH:MM:SS")
            
            # Show running optimizers if any
            running_str = length(running_optimizers) > 0 ? " | Running: $(join(sort(collect(running_optimizers)), ", "))" : ""
            
            if elapsed !== nothing
                println("[$timestamp] Progress: $completed/$total_optimizers ($progress_pct%) | $optimizer_name: $status (Time: $(round(elapsed, digits=1))s)$running_str")
            else
                println("[$timestamp] Progress: $completed/$total_optimizers ($progress_pct%) | $optimizer_name: $status$running_str")
            end
            flush(stdout)
        finally
            unlock(progress_lock)
        end
    end
    
    # PARALLEL EXECUTION - All optimizers run concurrently
    println("Launching all optimizers in parallel...")
    println()
    
    # Create tasks for all optimizers
    tasks = Task[]
    
    # Task 1: BlackBoxOptim
    push!(tasks, @async begin
        opt_start = time()
        update_progress("BlackBoxOptim", "Starting...", nothing, false, true)
        try
            best_perm, best_f, state, opt_result = optimize_material_ordering_blackbox(
                materials=materials,
                max_evaluations=max_evaluations,
                verbose=false,
                track_history=track_history
            )
            opt_elapsed = time() - opt_start
            update_progress("BlackBoxOptim", "✓ Completed", opt_elapsed, true, false)
            return (:BlackBoxOptim, (permutation=best_perm, force=best_f, state=state, result=opt_result, success=true))
        catch e
            opt_elapsed = time() - opt_start
            error_msg = string(e)
            if occursin('\n', error_msg)
                error_msg = split(error_msg, '\n')[1]
            end
            update_progress("BlackBoxOptim", "✗ Failed: $error_msg", opt_elapsed, true, false)
            return (:BlackBoxOptim, (success=false, error=e))
        end
    end)
    
    # Task 2: Optim.jl
    push!(tasks, @async begin
        opt_start = time()
        update_progress("Optim.jl", "Starting...", nothing, false, true)
        try
            best_perm, best_f, state, opt_result = optimize_material_ordering_optim(
                materials=materials,
                max_iterations=max_evaluations,
                verbose=false,
                track_history=track_history
            )
            opt_elapsed = time() - opt_start
            update_progress("Optim.jl", "✓ Completed", opt_elapsed, true, false)
            return (:Optim, (permutation=best_perm, force=best_f, state=state, result=opt_result, success=true))
        catch e
            opt_elapsed = time() - opt_start
            update_progress("Optim.jl", "✗ Failed", opt_elapsed, true, false)
            return (:Optim, (success=false, error=e))
        end
    end)
    
    # Task 3: Evolutionary.jl
    push!(tasks, @async begin
        opt_start = time()
        update_progress("Evolutionary.jl", "Starting...", nothing, false, true)
        try
            best_perm, best_f, state, opt_result = optimize_material_ordering_evolutionary(
                materials=materials,
                max_generations=div(max_evaluations, 2),
                population_size=10,
                verbose=false,
                track_history=track_history
            )
            opt_elapsed = time() - opt_start
            update_progress("Evolutionary.jl", "✓ Completed", opt_elapsed, true, false)
            return (:Evolutionary, (permutation=best_perm, force=best_f, state=state, result=opt_result, success=true))
        catch e
            opt_elapsed = time() - opt_start
            update_progress("Evolutionary.jl", "✗ Failed", opt_elapsed, true, false)
            return (:Evolutionary, (success=false, error=e))
        end
    end)
    
    # Tasks 4-7: Metaheuristics algorithms
    metaheuristics_algorithms = [:PSO, :DE, :ECA, :ES]
    for alg in metaheuristics_algorithms
        push!(tasks, @async begin
            alg_name = "Metaheuristics_$alg"
            opt_start = time()
            update_progress(alg_name, "Starting...", nothing, false, true)
            try
                best_perm, best_f, state, opt_result = optimize_material_ordering_metaheuristics(
                    materials=materials,
                    max_evaluations=max_evaluations,
                    verbose=false,
                    track_history=track_history,
                    algorithm=alg
                )
                opt_elapsed = time() - opt_start
                update_progress(alg_name, "✓ Completed", opt_elapsed, true, false)
                return (Symbol(alg_name), (permutation=best_perm, force=best_f, state=state, result=opt_result, success=true))
            catch e
                opt_elapsed = time() - opt_start
                error_msg = string(e)
                if occursin('\n', error_msg)
                    error_msg = split(error_msg, '\n')[1]
                end
                # Truncate very long error messages
                if length(error_msg) > 100
                    error_msg = error_msg[1:97] * "..."
                end
                
                if alg == :ES && occursin("not available", string(e))
                    update_progress(alg_name, "⚠ Skipped: $error_msg", opt_elapsed, true, false)
                    return (Symbol(alg_name), (success=false, error=e, skipped=true))
                else
                    update_progress(alg_name, "✗ Failed: $error_msg", opt_elapsed, true, false)
                    return (Symbol(alg_name), (success=false, error=e))
                end
            end
        end)
    end
    
    # Wait for all tasks to complete and collect results
    println("\nWaiting for all optimizers to complete...")
    # Use a channel or periodic updates to show progress while waiting
    for (i, task) in enumerate(tasks)
        try
            name, result = fetch(task)
            results[name] = result
            # Print a brief update when each task completes
            if haskey(result, :success) && result.success
                println("  ✓ $(name) completed")
            elseif haskey(result, :success) && !result.success
                println("  ✗ $(name) failed")
            end
        catch e
            println("  ✗ Task $i failed with error: $e")
        end
    end
    
    # Final summary
    total_elapsed = time() - overall_start_time
    println("\n" * "="^80)
    println("All optimizers completed!")
    println("Total time: $(round(total_elapsed, digits=1)) seconds ($(round(total_elapsed/60, digits=1)) minutes)")
    println("="^80)
    
    return results
end


########################################################################
#  PRINT COMPARISON SUMMARY
########################################################################
function print_comparison_summary(results)
    """
    Print a summary table comparing all optimizers.
    """
    println("\n" * "="^80)
    println("COMPARISON SUMMARY")
    println("="^80)
    println()
    println("Optimizer                    | Peak Force    | Evaluations | Success")
    println("-"^80)
    
    for (name, res) in results
        if res.success
            force = res.force
            evals = res.state.evaluation_count
            println(@sprintf("%-28s | %13.6f | %11d | ✓", name, force, evals))
        else
            println(@sprintf("%-28s | %13s | %11s | ✗", name, "N/A", "N/A"))
        end
    end
    
    println()
    
    # Find best overall
    successful_results = [(name, res) for (name, res) in results if res.success]
    if !isempty(successful_results)
        best_idx = argmin([res.force for (_, res) in successful_results])
        best_name, best_res = successful_results[best_idx]
        println("Best Overall: $best_name with peak force = $(best_res.force)")
        println("Best Permutation: $(best_res.permutation)")
    end
end


########################################################################
#  SAVE COMPARISON RESULTS
########################################################################
function save_comparison_results(results, filename)
    """
    Save comparison results to a text file.
    """
    open(filename, "w") do f
        println(f, "# Optimizer Comparison Results")
        println(f, "# Date: $(Dates.now())")
        println(f, "")
        
        for (name, res) in results
            println(f, "## $name")
            if res.success
                println(f, "Success: Yes")
                println(f, "Best Peak Force: $(res.force)")
                println(f, "Best Permutation: $(res.permutation)")
                println(f, "Evaluations: $(res.state.evaluation_count)")
                println(f, "")
            else
                println(f, "Success: No")
                println(f, "Error: $(res.error)")
                println(f, "")
            end
        end
    end
    @printf("Comparison results saved to %s\n", filename)
end


########################################################################
#  GENERATE CONVERGENCE PLOTS
########################################################################
function generate_convergence_plot(results, output_dir)
    """
    Generate convergence plots showing force vs evaluations for all optimizers.
    """
    fig = Figure(size = (1200, 800), fontsize = 14)
    ax = Axis(fig[1, 1],
              title = "Optimizer Convergence Comparison: Peak Force vs Evaluations",
              xlabel = "Function Evaluations",
              ylabel = "Peak Force (N)",
              xscale = log10)
    
    # Color palette for different optimizers
    colors = Dict(
        :BlackBoxOptim => :blue,
        :Optim => :red,
        :Evolutionary => :green,
        :Metaheuristics_PSO => :orange,
        :Metaheuristics_DE => :purple,
        :Metaheuristics_ECA => :brown,
        :Metaheuristics_ES => :pink
    )
    
    # Plot convergence history for each successful optimizer
    for (name, res) in results
        if res.success && !isempty(res.state.force_history)
            color = get(colors, name, :black)
            evaluations = 1:length(res.state.force_history)
            lines!(ax, evaluations, res.state.force_history,
                   color = color, linewidth = 2, label = string(name))
            # Add marker at best point
            best_idx = argmin(res.state.force_history)
            scatter!(ax, [evaluations[best_idx]], [res.state.force_history[best_idx]],
                    color = color, markersize = 10, marker = :circle)
        end
    end
    
    axislegend(ax, position = :rt)
    
    # Save figure
    output_path = joinpath(output_dir, "optimizer_convergence_comparison.png")
    try
        save(output_path, fig)
        @printf("Convergence plot saved to %s\n", output_path)
    catch e
        @warn "Failed to save convergence plot: $e"
    end
    
    return fig
end


########################################################################
#  GENERATE LATEX TABLE
########################################################################
function generate_latex_table(results, filename)
    """
    Generate a LaTeX table comparing all optimizers.
    """
    open(filename, "w") do f
        println(f, "% Optimizer Comparison Table")
        println(f, "% Generated: $(Dates.now())")
        println(f, "")
        println(f, "\\begin{table}[h]")
        println(f, "\\centering")
        println(f, "\\caption{Comparison of Optimization Algorithms for Material Ordering}")
        println(f, "\\label{tab:optimizer_comparison}")
        println(f, "\\begin{tabular}{lcccc}")
        println(f, "\\toprule")
        println(f, "Optimizer & Peak Force (N) & Evaluations & Time (s) & Success \\\\")
        println(f, "\\midrule")
        
        for (name, res) in results
            name_str = string(name)
            if res.success
                force = res.force
                evals = res.state.evaluation_count
                time_elapsed = time() - res.state.start_time
                println(f, "$name_str & $(@sprintf("%.2f", force)) & $evals & $(@sprintf("%.1f", time_elapsed)) & \\checkmark \\\\")
            else
                println(f, "$name_str & --- & --- & --- & \\times \\\\")
            end
        end
        
        println(f, "\\bottomrule")
        println(f, "\\end{tabular}")
        println(f, "\\end{table}")
    end
    @printf("LaTeX table saved to %s\n", filename)
end


########################################################################
#  EXPORT CSV DATA
########################################################################
function export_csv_data(results, output_dir)
    """
    Export comparison data to CSV files for easy analysis.
    """
    # Summary CSV
    summary_file = joinpath(output_dir, "optimizer_summary.csv")
    open(summary_file, "w") do f
        println(f, "Optimizer,Peak_Force_N,Evaluations,Time_s,Success")
        for (name, res) in results
            if res.success
                time_elapsed = time() - res.state.start_time
                println(f, "$(name),$(res.force),$(res.state.evaluation_count),$time_elapsed,true")
            else
                println(f, "$(name),,,,false")
            end
        end
    end
    @printf("Summary CSV saved to %s\n", summary_file)
    
    # Convergence history CSV (one file per optimizer)
    for (name, res) in results
        if res.success && !isempty(res.state.force_history)
            history_file = joinpath(output_dir, "convergence_$(name).csv")
            open(history_file, "w") do f
                println(f, "Evaluation,Peak_Force_N")
                for (i, force) in enumerate(res.state.force_history)
                    println(f, "$i,$force")
                end
            end
        end
    end
    @printf("Convergence history CSVs saved to %s\n", output_dir)
end


########################################################################
#  GENERATE SUMMARY STATISTICS
########################################################################
function generate_summary_statistics(results)
    """
    Calculate and print summary statistics.
    """
    successful_results = [(name, res) for (name, res) in results if res.success]
    
    if isempty(successful_results)
        println("No successful optimizations to analyze.")
        return
    end
    
    forces = [res.force for (_, res) in successful_results]
    evaluations = [res.state.evaluation_count for (_, res) in successful_results]
    times = [time() - res.state.start_time for (_, res) in successful_results]
    
    println("\n" * "="^80)
    println("SUMMARY STATISTICS")
    println("="^80)
    println("Successful Optimizations: $(length(successful_results))")
    println()
    println("Peak Force Statistics:")
    println("  Best: $(minimum(forces)) N")
    println("  Worst: $(maximum(forces)) N")
    println("  Mean: $(Statistics.mean(forces)) N")
    println("  Std Dev: $(Statistics.std(forces)) N")
    println()
    println("Evaluation Statistics:")
    println("  Min: $(minimum(evaluations))")
    println("  Max: $(maximum(evaluations))")
    println("  Mean: $(Statistics.mean(evaluations))")
    println()
    println("Time Statistics:")
    println("  Min: $(@sprintf("%.1f", minimum(times))) s")
    println("  Max: $(@sprintf("%.1f", maximum(times))) s")
    println("  Mean: $(@sprintf("%.1f", Statistics.mean(times))) s")
    println()
    
    # Find best optimizer
    best_idx = argmin([res.force for (_, res) in successful_results])
    best_name, best_res = successful_results[best_idx]
    println("Best Overall: $best_name")
    println("  Peak Force: $(best_res.force) N")
    println("  Permutation: $(best_res.permutation)")
    println("  Evaluations: $(best_res.state.evaluation_count)")
    println("  Time: $(@sprintf("%.1f", time() - best_res.state.start_time)) s")
end


########################################################################
#  HELPER FUNCTIONS FOR STATISTICS
########################################################################
# Using Statistics.mean and Statistics.std from the Statistics package


########################################################################
#  RUN COMPARISON
########################################################################
# Run if executed directly OR if included with MAX_EVALUATIONS set
# When included via -e 'MAX_EVALUATIONS=50; include(...)', PROGRAM_FILE is empty string ""
const CURRENT_FILE = @__FILE__
should_run = abspath(PROGRAM_FILE) == CURRENT_FILE || 
             ((PROGRAM_FILE == "" || PROGRAM_FILE == "none") && isdefined(Main, :MAX_EVALUATIONS))
if should_run
    # Create output directory
    output_dir = joinpath(@__DIR__, "comparison_output")
    mkpath(output_dir)
    
    # Get max evaluations from command line or environment
    max_evals = get_max_evaluations()
    println("Using max_evaluations = $max_evals")
    println("(Override with: julia compare_optimizers.jl <number> or set MAX_EVALUATIONS env var)")
    println()
    
    # Run comparison
    results = compare_all_optimizers(max_evaluations=max_evals)
    
    # Print summary
    print_comparison_summary(results)
    
    # Generate summary statistics
    generate_summary_statistics(results)
    
    # Save text results
    result_file = joinpath(output_dir, "comparison_results.txt")
    save_comparison_results(results, result_file)
    
    # Generate convergence plot
    println("\nGenerating convergence plot...")
    generate_convergence_plot(results, output_dir)
    
    # Generate LaTeX table
    println("\nGenerating LaTeX table...")
    latex_file = joinpath(output_dir, "optimizer_comparison_table.tex")
    generate_latex_table(results, latex_file)
    
    # Export CSV data
    println("\nExporting CSV data...")
    export_csv_data(results, output_dir)
    
    println("\n" * "="^80)
    println("Comparison complete!")
    println("="^80)
    println("All outputs saved to: $output_dir")
    println("  - Text summary: comparison_results.txt")
    println("  - Convergence plot: optimizer_convergence_comparison.png")
    println("  - LaTeX table: optimizer_comparison_table.tex")
    println("  - CSV data: optimizer_summary.csv + convergence_*.csv")
end

