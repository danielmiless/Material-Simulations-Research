# Material Ordering Optimization

This directory contains scripts for optimizing the ordering of materials across columns to minimize force transferred to the backplate.

## Quick Start

### Run Optimizer Comparison

The easiest way to run all optimizers and compare results:

```bash
# Quick test (2-5 evaluations, takes ~5-15 minutes)
./scripts/optimization/run_overnight.sh 5

# Full run (50-100 evaluations, takes ~4-7 hours)
./scripts/optimization/run_overnight.sh 50
```

Or directly with Julia:

```bash
# Set max evaluations
MAX_EVALUATIONS=20 julia --project=. -e 'include("scripts/optimization/compare_optimizers.jl")'

# Or pass as argument (if script supports it)
julia --project=. scripts/optimization/compare_optimizers.jl
```

### Run Individual Optimizer

```bash
# BlackBoxOptim
julia --project=. scripts/optimization/optimize_blackbox.jl

# Optim.jl
julia --project=. scripts/optimization/optimize_optim.jl

# Evolutionary.jl
julia --project=. scripts/optimization/optimize_evolutionary.jl

# Metaheuristics.jl (PSO, DE, ECA, ES)
julia --project=. scripts/optimization/optimize_metaheuristics.jl
```

## Files

### Core Scripts

- **`compare_optimizers.jl`**: Main comparison script that runs all optimizers in parallel and generates comparison reports
- **`optimization_utils.jl`**: Shared utilities (objective functions, permutation handling, result tracking)
- **`run_overnight.sh`**: Automated script for running comparisons with logging and output management

### Individual Optimizers

- **`optimize_blackbox.jl`**: BlackBoxOptim.jl optimizer
- **`optimize_optim.jl`**: Optim.jl optimizer (SimulatedAnnealing/NelderMead)
- **`optimize_evolutionary.jl`**: Evolutionary.jl genetic algorithm optimizer
- **`optimize_metaheuristics.jl`**: Metaheuristics.jl optimizer (PSO, DE, ECA, ES)

### Utilities

- **`test_basic.jl`**: Basic setup test to verify installation
- **`check_progress.sh`**: Script to check status of running optimizations

## Configuration

### Setting Max Evaluations

The number of function evaluations can be set in several ways:

1. **Via run_overnight.sh argument**:
   ```bash
   ./scripts/optimization/run_overnight.sh 50
   ```

2. **Via environment variable**:
   ```bash
   MAX_EVALUATIONS=50 julia --project=. -e 'include("scripts/optimization/compare_optimizers.jl")'
   ```

3. **Edit constants in scripts** (not recommended):
   - Each optimizer script has a `MAX_EVALUATIONS` or `MAX_ITERATIONS` constant at the top

### Algorithm Requirements

Different algorithms have different minimum evaluation requirements:

- **PSO/ES**: Minimum 2 evaluations
- **DE/ECA**: Minimum 5 evaluations (4 for population + 1 for iteration)
- **BlackBoxOptim**: No minimum (but needs at least 1)
- **Optim.jl**: No minimum (but needs at least 1)
- **Evolutionary.jl**: Depends on population size

**Recommendation**: Use at least 5 evaluations to test all algorithms, or 20-50 for meaningful comparisons.

## Output

Results are saved in `comparison_output/`:

- **`comparison_results.txt`**: Text summary of all optimizers
- **`optimizer_summary.csv`**: CSV summary for analysis
- **`optimizer_comparison_table.tex`**: LaTeX table for papers
- **`optimizer_convergence_comparison.png`**: Convergence plot
- **`convergence_*.csv`**: Individual convergence histories
- **`overnight_run_YYYYMMDD_HHMMSS/`**: Timestamped directories with full run logs

## Parallel Execution

All optimizers run in **parallel by default** using Julia's `@async` tasks. This provides ~5-7x speedup compared to sequential execution.

- All 7 optimizers run concurrently
- Progress updates appear in real-time
- Total time ≈ time of longest optimizer (not sum of all)

## Dependencies

Required packages (install in your Julia environment):

```julia
using Pkg
Pkg.add("Optim")
Pkg.add("Evolutionary")
Pkg.add("Metaheuristics")
Pkg.add("BlackBoxOptim")  # Optional - may need: Pkg.add(url="https://github.com/robertfeldt/BlackBoxOptim.jl")
```

**Note**: If you encounter registry errors, see `INSTALL_OPTIMIZATION_PACKAGES.md` in the project root for alternative installation methods.

## Testing

Before running optimizers, test the basic setup:

```bash
julia --project=. scripts/optimization/test_basic.jl
```

This verifies:
- Simulation loads correctly
- Material ordering function works
- Which optimizers are available

## Troubleshooting

### Optimizers Exceeding Evaluation Limits

If DE or ECA run more evaluations than specified, this has been fixed in the latest version. The fix enforces a hard limit in the objective function.

### Progress Updates Not Appearing

Progress updates are printed in real-time. If they appear batched, check that stdout is not being buffered.

### Algorithms Failing

- **ES/CMA_ES**: May not be available in all Metaheuristics.jl versions - this is expected
- **Optim.jl**: May fail if SimulatedAnnealing is not available - falls back to NelderMead
- **Evolutionary.jl**: Uses default GA operators if custom ones fail

## Documentation

For detailed information, see:

- **`RUN_INSTRUCTIONS.md`**: Detailed usage instructions
- **`PARALLEL_USAGE.md`**: Parallel execution details (now default)

## Performance Tips

1. **For testing**: Use 2-5 evaluations (completes in ~5-15 minutes)
2. **For comparisons**: Use 20-50 evaluations (completes in ~1-4 hours)
3. **For final results**: Use 50-100+ evaluations (completes in ~4-7+ hours)

The parallel execution means all optimizers run simultaneously, so total time is approximately the time of the slowest optimizer.
