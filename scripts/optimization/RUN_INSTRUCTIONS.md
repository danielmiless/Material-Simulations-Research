# Running Optimizer Comparisons

This guide explains how to run optimizer comparisons efficiently using the automated scripts.

## Quick Start

### Using the Overnight Script (Recommended)

The easiest way to run comparisons with full logging and output management:

```bash
# Quick test (5 evaluations, takes ~5-15 minutes)
./scripts/optimization/run_overnight.sh 5

# Standard run (20 evaluations, takes ~1-2 hours)
./scripts/optimization/run_overnight.sh 20

# Full run (50-100 evaluations, takes ~4-7 hours)
./scripts/optimization/run_overnight.sh 50
```

This script:
- Creates timestamped output directories
- Logs everything to a file
- Runs all optimizers in parallel (default)
- Saves all results with timestamps
- Provides real-time progress updates

### Direct Julia Execution

Run directly with Julia:

```bash
# Via environment variable
MAX_EVALUATIONS=20 julia --project=. -e 'include("scripts/optimization/compare_optimizers.jl")'

# Or edit the script to set MAX_EVALUATIONS constant
julia --project=. scripts/optimization/compare_optimizers.jl
```

## Parallel Execution

**All optimizers run in parallel by default** using Julia's `@async` tasks. This provides ~5-7x speedup compared to sequential execution.

- All 7 optimizers run concurrently
- Progress updates appear in real-time
- Total time ≈ time of longest optimizer (not sum of all)

No configuration needed - parallel execution is automatic!

## Algorithm Requirements

Different algorithms have different minimum evaluation requirements:

- **PSO/ES**: Minimum 2 evaluations
- **DE/ECA**: Minimum 5 evaluations (4 for population + 1 for iteration)
- **BlackBoxOptim**: No minimum (but needs at least 1)
- **Optim.jl**: No minimum (but needs at least 1)
- **Evolutionary.jl**: Depends on population size

**Recommendation**: Use at least 5 evaluations to test all algorithms, or 20-50 for meaningful comparisons.

## Output Files

All results are saved to `scripts/optimization/comparison_output/`:

### Main Output Files
- `comparison_results.txt` - Text summary of all optimizers
- `optimizer_summary.csv` - CSV summary for analysis
- `optimizer_comparison_table.tex` - LaTeX table for papers
- `optimizer_convergence_comparison.png` - Convergence plot
- `convergence_*.csv` - Individual convergence histories

### Timestamped Run Directories
- `overnight_run_YYYYMMDD_HHMMSS/` - Full run logs and results

## Time Estimates

Based on ~30-50 seconds per evaluation per optimizer:

| Evaluations | Sequential (if used) | Parallel (default) |
|-------------|---------------------|-------------------|
| 5           | ~15-25 minutes      | ~3-5 minutes      |
| 20          | ~1-2 hours          | ~10-20 minutes    |
| 50          | ~4-7 hours          | ~30-60 minutes    |
| 100         | ~8-14 hours         | ~1-2 hours        |

**Note**: Parallel execution time is approximately the time of the slowest optimizer, not the sum.

## Long-Running Jobs

For runs that take hours, use one of these methods:

### Option A: Using `nohup` (Unix/Linux/macOS)

Run in background and detach from terminal:

```bash
nohup ./scripts/optimization/run_overnight.sh 50 > optimization.log 2>&1 &
```

Check progress:
```bash
tail -f optimization.log
```

### Option B: Using `screen` or `tmux` (Recommended for Servers)

**Using screen:**
```bash
screen -S optimization
./scripts/optimization/run_overnight.sh 50
# Press Ctrl+A then D to detach
# Reattach with: screen -r optimization
```

**Using tmux:**
```bash
tmux new -s optimization
./scripts/optimization/run_overnight.sh 50
# Press Ctrl+B then D to detach
# Reattach with: tmux attach -t optimization
```

### Option C: Check Progress

Use the progress check script:

```bash
./scripts/optimization/check_progress.sh
```

This shows:
- Running processes
- Screen sessions
- Latest output directories
- Recent log entries

## Tips

1. **Start small**: Test with 5 evaluations first to verify everything works
2. **Monitor progress**: Progress updates appear in real-time
3. **Check disk space**: Results can be several MB per optimizer
4. **Use server**: For long runs, use a server or machine that won't sleep
5. **Review logs**: Check the timestamped run directories for detailed logs

## Troubleshooting

### Progress Updates Not Appearing

Progress updates are printed in real-time. If they appear batched, check that stdout is not being buffered.

### Algorithms Failing

- **ES/CMA_ES**: May not be available in all Metaheuristics.jl versions - this is expected
- **Optim.jl**: May fail if SimulatedAnnealing is not available - falls back to NelderMead
- **Evolutionary.jl**: Uses default GA operators if custom ones fail

### Evaluation Limits

If DE or ECA run more evaluations than specified, ensure you're using the latest version which enforces hard limits.

## See Also

- **`README.md`**: Overview and quick start
- **`PARALLEL_USAGE.md`**: Parallel execution details (now default)
