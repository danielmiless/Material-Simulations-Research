# Optimizer Comparison - Parallel Execution

## Overview

The optimizer comparison script uses **parallel execution by default**, running all optimizers concurrently to significantly speed up the comparison.

## Performance

- **Total time**: ≈ time of longest optimizer (~1-3 minutes for 2 evaluations)
- **Speedup**: ~5-7x faster than sequential execution would be
- All 7 optimizers run simultaneously using Julia's `@async` tasks

## Usage

Simply run the script - parallel execution is the default and only mode:

```bash
# Run comparison with parallel execution
./scripts/optimization/run_overnight.sh 2

# Or directly with Julia
julia --project=. -e 'MAX_EVALUATIONS=2; include("scripts/optimization/compare_optimizers.jl")'
```

## How It Works

1. **Parallel Execution**: Uses Julia's `@async` tasks to launch all 7 optimizers concurrently
2. **Thread-Safe Progress**: Uses `Atomic` counters and `ReentrantLock` for thread-safe progress tracking
3. **Result Collection**: Waits for all tasks to complete using `fetch()` and collects results

## Implementation Details

- All optimizers run in separate async tasks
- Progress updates are synchronized with locks
- Each optimizer's output may interleave (this is expected)
- Final results are collected and displayed in a consistent format

