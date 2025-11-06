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

This will install:
- `DifferentialEquations` - For solving ODEs
- `GLMakie` - For visualization and animation
- `LinearAlgebra` - Standard library for linear algebra
- `Printf` - Standard library for formatted output

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

The main simulation scripts are located in `src/`:

#### Exponential Springs Model

```bash
julia --project=. src/exponential-springs.jl
```

This runs the exponential spring model with:
- Exponential spring force law
- Damping on nearest-neighbor springs
- Configurable external force
- Real-time animation

#### Linear Springs Model

```bash
julia --project=. src/lattice_w_anim.jl
```

This runs the linear (Hooke's law) spring model with:
- Standard linear springs
- Diagonal spring connections
- Animation output

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
include("src/exponential-springs.jl")
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

### Getting Help

- Check the [API Reference](api-reference.md) for function documentation
- Review the [Models](models.md) documentation for theoretical background
- Contact the research team for assistance

## Next Steps

- Read the [Methodology](methodology.md) to understand the research approach
- Explore the [Models](models.md) documentation for detailed model descriptions
- Check the [API Reference](api-reference.md) for code documentation

