# Getting Started

This guide will help you set up the development environment and run your first simulation.

## Prerequisites

Before you begin, ensure you have the following installed:

- **Julia 1.9 or later**: Download from [julialang.org](https://julialang.org/downloads/)
- **Git**: For version control (usually pre-installed on macOS/Linux)
- **LaTeX** (optional): For compiling papers in the `papers/` directory

## Installation

### 1. Clone the Repository

```bash
git clone <repository-url>
cd Material-Simulations-Research
```

### 2. Activate the Julia Environment

The project includes a `Project.toml` file that specifies all dependencies. Activate the project environment:

```bash
julia --project=.
```

Or from within Julia:

```julia
using Pkg
Pkg.activate(".")
```

### 3. Install Dependencies

Install all required packages:

```julia
using Pkg
Pkg.instantiate()
```

This will install all required packages, including:
- `DifferentialEquations` - For solving ODEs
- `GLMakie` / `CairoMakie` - For visualization and animation
- `Optim`, `Evolutionary`, `Metaheuristics`, `BlackBoxOptim` - For optimization
- `LinearAlgebra`, `Printf`, `Statistics`, `Dates`, `DelimitedFiles` - Standard libraries

### 4. Verify Installation

Test that everything is working:

```julia
using DifferentialEquations
using GLMakie
using LinearAlgebra

println("Installation successful!")
```

## Running Simulations

### Basic Usage

The main simulation script is located in `src/`:

#### 11×11 Lattice Simulation

```bash
julia --project=. src/lattice_simulation_11x11.jl
```

This runs the 11×11 lattice simulation with:
- Exponential spring force law
- Damping on nearest-neighbor springs
- Backplate boundary conditions
- Column-based material property scaling
- Distributed load application
- Real-time animation

**Note**: Older simulation versions (5×5 and 10×10 systems) are available in `Deprecated Scripts/` for reference.

### Configuration

Both scripts contain configuration constants at the top that can be modified:

- **Physical Parameters**: Mass, spring constants, damping coefficients
- **Lattice Size**: Grid dimensions (currently 5×5)
- **Force Configuration**: Magnitude, angle, application point
- **Numerical Parameters**: Time span, tolerances, output intervals
- **Animation Parameters**: Speed, node sizes, colors

### Example: Changing Force Configuration

Edit the constants in the script:

```julia
const FORCE_ANGLE_DEGREES = 45.0    # Change angle
const FORCE_TARGET_ROW = 2          # Change row
const FORCE_TARGET_COL = 3          # Change column
const F_MAG = 500.0                  # Change magnitude
```

## Development Workflow

### Using the Julia REPL

For interactive development:

```bash
julia --project=.
```

Then load modules:

```julia
include("src/lattice_simulation_11x11.jl")
```

### Running Tests

Run the test suite:

```bash
julia --project=. tests/runtests.jl
```

### Creating Animations

Animations are automatically saved to `animations/` when running simulations. The `save_animation_to_file()` function creates MP4 files.

## Troubleshooting

### Common Issues

**Issue**: `GLMakie` fails to load
- **Solution**: Ensure you have OpenGL support. On macOS, this is usually available. On Linux, you may need to install OpenGL libraries.

**Issue**: Dependencies fail to install
- **Solution**: Try updating Julia's package registry:
  ```julia
  using Pkg
  Pkg.update()
  Pkg.instantiate()
  ```

**Issue**: Animation window doesn't appear
- **Solution**: Check that your display server is running (X11 on Linux, or standard macOS display).

### Running Optimization

For material ordering optimization:

```bash
# Quick test
./scripts/optimization/run_overnight.sh 5

# Full optimization run
./scripts/optimization/run_overnight.sh 50
```

See the [Optimization Documentation](optimization.md) for details.

### Getting Help

- Check the [API Reference](api-reference.md) for function documentation
- Review the [Models](models.md) documentation for theoretical background
- See [Optimization](optimization.md) for optimization usage
- Contact the research team for assistance

## Next Steps

- Read the [Methodology](methodology.md) to understand the research approach
- Explore the [Models](models.md) documentation for detailed model descriptions
- Check the [API Reference](api-reference.md) for code documentation
- Review [Optimization](optimization.md) for material ordering optimization

