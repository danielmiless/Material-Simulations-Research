# Running Optimizer Comparisons

This guide explains how to run optimizer comparisons efficiently, including parallel execution and overnight runs.

## Quick Start

### Standard Comparison (Sequential)

Run all optimizers sequentially with default settings (20 evaluations):

```bash
julia --project=. scripts/optimization/compare_optimizers.jl
```

### Custom Max Evaluations

Specify the number of evaluations:

```bash
# Via command line argument
julia --project=. scripts/optimization/compare_optimizers.jl 50

# Via environment variable
MAX_EVALUATIONS=50 julia --project=. scripts/optimization/compare_optimizers.jl
```

## Recommendation 1: Reduced Evaluations for Testing

The default `MAX_EVALUATIONS` is now set to **20** (reduced from 50) for faster testing.

- **For testing**: Use 20 evaluations (default) - takes ~1-2 hours
- **For final runs**: Use 50-100 evaluations - takes ~4-7 hours

## Recommendation 2: Parallel Execution

For faster execution, you can run optimizers in separate terminal windows/tabs:

### Option A: Manual Parallel Execution

Open multiple terminals and run each optimizer separately:

**Terminal 1:**
```bash
julia --project=. scripts/optimization/optimize_optim.jl 50
```

**Terminal 2:**
```bash
julia --project=. scripts/optimization/optimize_evolutionary.jl 50
```

**Terminal 3:**
```bash
julia --project=. scripts/optimization/optimize_metaheuristics.jl 50 PSO
```

**Terminal 4:**
```bash
julia --project=. scripts/optimization/optimize_metaheuristics.jl 50 DE
```

Then combine results manually or use the comparison script to generate final plots.

### Option B: Using the Parallel Script

```bash
./scripts/optimization/run_parallel_comparison.sh 50
```

Note: This script runs optimizers in background processes. Monitor progress with:
```bash
tail -f scripts/optimization/comparison_output/temp_parallel_*/optim.log
```

## Recommendation 4: Overnight/Server Runs

### Option A: Using the Overnight Script

Run with higher evaluations and full logging:

```bash
# Default: 50 evaluations
./scripts/optimization/run_overnight.sh

# Custom: 100 evaluations
./scripts/optimization/run_overnight.sh 100
```

This script:
- Creates timestamped output directories
- Logs everything to a file
- Has a 24-hour timeout
- Saves all results with timestamps

### Option B: Using `nohup` (Unix/Linux/macOS)

Run in background and detach from terminal:

```bash
nohup julia --project=. scripts/optimization/compare_optimizers.jl 50 > optimization.log 2>&1 &
```

Check progress:
```bash
tail -f optimization.log
```

### Option C: Using `screen` or `tmux` (Recommended for Servers)

**Using screen:**
```bash
screen -S optimization
julia --project=. scripts/optimization/compare_optimizers.jl 50
# Press Ctrl+A then D to detach
# Reattach with: screen -r optimization
```

**Using tmux:**
```bash
tmux new -s optimization
julia --project=. scripts/optimization/compare_optimizers.jl 50
# Press Ctrl+B then D to detach
# Reattach with: tmux attach -t optimization
```

### Option D: Using `at` (Schedule for Later)

Schedule to run at a specific time:

```bash
# Run at 10 PM tonight
echo "cd $(pwd) && julia --project=. scripts/optimization/compare_optimizers.jl 50" | at 22:00

# Run tomorrow at 2 AM
echo "cd $(pwd) && julia --project=. scripts/optimization/compare_optimizers.jl 50" | at 2:00 tomorrow
```

## Output Files

All results are saved to `scripts/optimization/comparison_output/`:

- `comparison_results.txt` - Text summary
- `optimizer_convergence_comparison.png` - Convergence plot
- `optimizer_comparison_table.tex` - LaTeX table
- `optimizer_summary.csv` - Summary CSV
- `convergence_*.csv` - Individual convergence histories

## Time Estimates

Based on ~30-50 seconds per evaluation:

| Evaluations | Estimated Time |
|-------------|----------------|
| 20          | 1-2 hours      |
| 50          | 4-7 hours      |
| 100         | 8-14 hours     |

With parallel execution, you can reduce total time by running optimizers simultaneously.

## Tips

1. **Start small**: Test with 10-20 evaluations first
2. **Monitor progress**: Check logs regularly to ensure everything is working
3. **Save frequently**: Results are saved incrementally
4. **Use server**: For long runs, use a server or machine that won't sleep
5. **Check disk space**: Results can be several MB per optimizer

