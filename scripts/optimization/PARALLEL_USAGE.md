# Parallel Execution (Default)

## Overview

The optimizer comparison script uses **parallel execution by default**, running all optimizers concurrently to significantly speed up the comparison.

## Performance

- **Total time**: ≈ time of longest optimizer
- **Speedup**: ~5-7x faster than sequential execution would be
- All 7 optimizers run simultaneously using Julia's `@async` tasks

## Usage

Parallel execution is automatic - no configuration needed:

```bash
# Run comparison (parallel by default)
./scripts/optimization/run_overnight.sh 20

# Or directly with Julia
MAX_EVALUATIONS=20 julia --project=. -e 'include("scripts/optimization/compare_optimizers.jl")'
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

## Time Estimates

| Evaluations | Parallel Time (approx) |
|-------------|----------------------|
| 5           | ~3-5 minutes         |
| 20          | ~10-20 minutes       |
| 50          | ~30-60 minutes       |
| 100         | ~1-2 hours           |

Time is approximately the time of the slowest optimizer, not the sum of all.

