# Material Simulations Research

**Work in Progress**

This repository contains computational research on mass-spring lattice systems with exponential spring models, conducted in collaboration with Dr. Romesh Batra, Distinguished Professor and Clifton C. Garvin Professor, Department of Mechanical Engineering, Virginia Tech.

## Affiliation

- **Researcher**: Daniel Miles (Undergraduate Student)
- **Advisor**: Dr. Romesh Batra, Department of Mechanical Engineering, Virginia Tech
- **Institution**: Virginia Tech

## Overview

This project investigates the dynamics of 2D mass-spring lattice systems using exponential spring models. The simulations explore:

- Exponential spring force models with configurable decay rates
- 2D lattice systems with nearest-neighbor and diagonal spring connections
- Viscous damping on nearest-neighbor springs (energy dissipation)
- Immovable boundary conditions (backplate constraints)
- Energy conservation, dissipation, and work calculations
- Configurable external force application
- Real-time visualization and animation

## Quick Start

### Prerequisites

- Julia 1.9 or later
- Git

### Installation

1. Clone this repository:
```bash
git clone <repository-url>
cd Material-Simulations-Research
```

2. Activate the Julia environment:
```bash
julia --project=.
```

3. Install dependencies (if not already installed):
```julia
using Pkg
Pkg.instantiate()
```

### Running Simulations

The main simulation scripts are located in `src/`:

- **`lattice_simulation.jl`**: Main simulation with exponential springs, viscous damping, and configurable forces
- **`lattice_simulation_with_backplate.jl`**: Same as above, with an immovable backplate constraint

To run a simulation:

```bash
julia --project=. src/lattice_simulation.jl
```

or

```bash
julia --project=. src/lattice_simulation_with_backplate.jl
```

**Note**: Older deprecated scripts are available in `Deprecated Scripts/` for reference.

## Repository Structure

```
Material-Simulations-Research/
├── src/                    # Main simulation code
│   ├── lattice_simulation.jl
│   └── lattice_simulation_with_backplate.jl
├── docs/                    # Comprehensive documentation
├── papers/                  # LaTeX papers and research updates
│   ├── updates/            # Progress updates
│   ├── figures/            # Figures for papers
│   └── build/              # LaTeX build artifacts (gitignored)
├── scripts/                 # Utility scripts
│   └── generate_comparison_data.jl  # Script for generating comparison figures
├── Deprecated Scripts/      # Older versions for reference
├── tests/                   # Test suite
├── animations/              # Generated animations
└── data/                    # Data files
```

## Documentation

For detailed documentation, see the [docs/](docs/) directory:

- [Getting Started](docs/getting-started.md) - Setup and installation guide
- [Methodology](docs/methodology.md) - Research methodology and theoretical background
- [Models](docs/models.md) - Detailed model descriptions
- [Results](docs/results.md) - Results and analysis (work in progress)
- [API Reference](docs/api-reference.md) - Code documentation

## Current Status

⚠️ **Work in Progress**: This repository is actively under development. The codebase and documentation are being refined as research progresses.

## Research Objectives

- Investigate dynamics of exponential spring-mass systems
- Analyze energy dissipation through viscous damping
- Study wave propagation and reflection at boundaries
- Analyze response to external forces and boundary conditions
- Develop computational tools for material simulation

## Papers

Research papers and progress updates will be added to the `papers/` directory as they are developed.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

For questions or collaboration inquiries, please contact:
- Daniel Miles: danielmiles@vt.edu
- Dr. Romesh Batra: Department of Mechanical Engineering, Virginia Tech

## Acknowledgments

This research is conducted under the guidance of Dr. Romesh Batra, Distinguished Professor and Clifton C. Garvin Professor, Department of Mechanical Engineering, Virginia Tech.

