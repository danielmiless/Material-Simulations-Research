# Installing Optimization Packages

This document provides instructions for installing the optimization packages required for the material ordering optimization scripts.

## Standard Installation (Recommended)

The project's `Project.toml` file already includes all required dependencies. To install them:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

This will install all packages listed in `Project.toml`, including:
- `Optim` - For Simulated Annealing and Nelder-Mead
- `Evolutionary` - For Genetic Algorithms
- `Metaheuristics` - For PSO, DE, ECA, and ES algorithms
- `BlackBoxOptim` - For adaptive differential evolution
- Other dependencies (CairoMakie, Dates, DelimitedFiles, Statistics, etc.)

## Alternative Installation Methods

If you encounter registry issues or need to install packages manually:

### Method 1: Install in Default Environment

Install the packages in your default Julia environment:

```julia
using Pkg
Pkg.add("Optim")
Pkg.add("Evolutionary")
Pkg.add("Metaheuristics")
Pkg.add("BlackBoxOptim")  # Optional - may need to install from GitHub
```

If BlackBoxOptim fails, try:
```julia
Pkg.add(url="https://github.com/robertfeldt/BlackBoxOptim.jl")
```

### Method 2: Manual Installation Script

Create a file `install_packages.jl`:

```julia
using Pkg

# Activate project
Pkg.activate(".")

# Add packages one by one
try
    Pkg.add("Optim")
    println("✓ Optim installed")
catch e
    println("✗ Optim failed: $e")
end

try
    Pkg.add("Evolutionary")
    println("✓ Evolutionary installed")
catch e
    println("✗ Evolutionary failed: $e")
end

try
    Pkg.add("Metaheuristics")
    println("✓ Metaheuristics installed")
catch e
    println("✗ Metaheuristics failed: $e")
end

try
    Pkg.add("BlackBoxOptim")
    println("✓ BlackBoxOptim installed")
catch e
    println("✗ BlackBoxOptim failed, trying GitHub: $e")
    try
        Pkg.add(url="https://github.com/robertfeldt/BlackBoxOptim.jl")
        println("✓ BlackBoxOptim installed from GitHub")
    catch e2
        println("✗ BlackBoxOptim completely failed: $e2")
    end
end
```

Run with: `julia install_packages.jl`

### Method 3: Use Without Project File

You can run the optimization scripts without the project file - they'll use packages from your default environment:

```bash
julia --project=@. scripts/optimization/optimize_optim.jl
```

## Testing Installation

After installation, verify the setup:

```bash
julia --project=. scripts/optimization/test_basic.jl
```

This will test:
- Simulation loading
- Material ordering function
- Optimization utilities
- Available optimizers

## Testing Without All Packages

The comparison script handles missing packages gracefully. You can test with just the packages that install successfully:

1. Test Optim.jl: `julia --project=. scripts/optimization/optimize_optim.jl`
2. Test Evolutionary.jl: `julia --project=. scripts/optimization/optimize_evolutionary.jl`
3. Test Metaheuristics.jl: `julia --project=. scripts/optimization/optimize_metaheuristics.jl`

The comparison script will skip any unavailable optimizers automatically.

## Troubleshooting

If you encounter registry errors or package installation issues:

1. **Try standard installation first**: `Pkg.activate("."); Pkg.instantiate()`
2. **Check Julia version**: Requires Julia 1.9 or later
3. **Update package registry**: `using Pkg; Pkg.update()`
4. **Clear package cache**: Delete `~/.julia/packages` and reinstall (last resort)

