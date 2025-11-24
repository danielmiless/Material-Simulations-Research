# Material Ordering Optimization Scripts

This directory contains scripts for optimizing the ordering of materials across columns to minimize force transferred to the backplate.

## Files

- `optimization_utils.jl`: Shared utilities for optimization (objective functions, permutation handling, result tracking)
- `optimize_blackbox.jl`: BlackBoxOptim.jl optimizer
- `optimize_optim.jl`: Optim.jl optimizer
- `optimize_evolutionary.jl`: Evolutionary.jl genetic algorithm optimizer
- `optimize_metaheuristics.jl`: Metaheuristics.jl optimizer (tests multiple algorithms)
- `compare_optimizers.jl`: Comparison script that runs all optimizers and compares results

## Quick Start

### Run a single optimizer:

```julia
julia scripts/optimization/optimize_blackbox.jl
```

### Compare all optimizers:

```julia
julia scripts/optimization/compare_optimizers.jl
```

## Configuration

Edit the constants at the top of each script to adjust:
- Maximum evaluations/iterations
- Verbosity
- History tracking

## Results

Results are saved in this directory:
- Individual results: `results_<optimizer>.txt`
- Comparison results: `comparison_results.txt`

## Dependencies

Make sure all optimization packages are installed. Due to potential registry issues, install them in your default Julia environment:

```julia
using Pkg
Pkg.add("Optim")
Pkg.add("Evolutionary")
Pkg.add("Metaheuristics")
Pkg.add("BlackBoxOptim")  # Optional - may need: Pkg.add(url="https://github.com/robertfeldt/BlackBoxOptim.jl")
```

**Note**: If you encounter registry errors, see `INSTALL_OPTIMIZATION_PACKAGES.md` in the project root for alternative installation methods.

## Quick Test

Before running optimizers, test the basic setup:

```bash
julia scripts/optimization/test_basic.jl
```

This will verify:
- Simulation loads correctly
- Material ordering function works
- Which optimizers are available

