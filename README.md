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
- Column-based material property scaling (heterogeneous materials)
- Energy conservation, dissipation, and work calculations
- Configurable external force application
- Parameter sweep analysis for wall properties
- Real-time visualization and animation
- Multiple system sizes (5×5, 10×10, and 11×11 lattices)

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

The main simulation script is located in `src/`:

- **`lattice_simulation_11x11.jl`**: 11×11 simulation with backplate, material scaling, and distributed load

To run the simulation:

```bash
julia --project=. src/lattice_simulation_11x11.jl
```

**Note**: Older simulation versions (5×5 and 10×10 systems) are available in `Deprecated Scripts/` for reference.

## Repository Structure

```
Material-Simulations-Research/
├── src/                    # Main simulation code
│   └── lattice_simulation_11x11.jl       # 11×11 with backplate, material scaling, and distributed load
├── docs/                    # Comprehensive documentation
├── papers/                  # LaTeX papers and research updates
│   ├── updates/            # Progress updates
│   ├── figures/            # Figures for papers
│   └── build/              # LaTeX build artifacts (gitignored)
├── scripts/                 # Utility scripts
│   ├── generate_comparison_data.jl        # Generate comparison figures
│   ├── generate_update_paper_materials.jl # Generate all materials for update papers
│   ├── generate_distributed_load_figures.jl # Generate figures for distributed load paper
│   ├── parameter_sweep_wall_properties.jl # Parameter sweep analysis
│   └── optimization/       # Material ordering optimization scripts
│       ├── compare_optimizers.jl          # Main comparison script (runs all optimizers)
│       ├── optimize_*.jl                  # Individual optimizer scripts
│       ├── run_overnight.sh               # Automated comparison runner
│       └── README.md                      # Optimization documentation
├── Deprecated Scripts/      # Older simulation versions (5×5, 10×10) for reference
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
- Model heterogeneous materials with column-based property scaling
- Perform parameter sensitivity analysis for boundary conditions
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

