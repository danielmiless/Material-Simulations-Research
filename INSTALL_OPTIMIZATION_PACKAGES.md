# Installing Optimization Packages

There appears to be a registry issue with Julia 1.11.6. Here are several ways to install the optimization packages:

## Method 1: Install in Default Environment (Recommended)

Install the packages in your default Julia environment, then they'll be available when running scripts:

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

## Method 2: Remove UUIDs from Project.toml

Edit `Project.toml` to remove the UUIDs and let Pkg resolve them automatically:

```toml
[deps]
Optim = ""
Evolutionary = ""
Metaheuristics = ""
# ... other packages
```

Then run:
```julia
using Pkg
Pkg.activate(".")
Pkg.resolve()
Pkg.instantiate()
```

## Method 3: Manual Installation Script

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

## Method 4: Use Without Project File

You can run the optimization scripts without the project file - they'll use packages from your default environment:

```bash
julia --project=@. scripts/optimization/optimize_optim.jl
```

## Testing Without All Packages

The comparison script handles missing packages gracefully. You can test with just the packages that install successfully:

1. Test Optim.jl: `julia scripts/optimization/optimize_optim.jl`
2. Test Evolutionary.jl: `julia scripts/optimization/optimize_evolutionary.jl`
3. Test Metaheuristics.jl: `julia scripts/optimization/optimize_metaheuristics.jl`

The comparison script will skip BlackBoxOptim if it's not available.

