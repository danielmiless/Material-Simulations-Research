# Material Simulations Research

**Work in Progress**

This repository contains computational research on mass-spring lattice systems with exponential spring models, conducted in collaboration with [Dr. Batra](https://www.me.vt.edu/people/faculty/rakesh-batra.html), Professor of Mechanical Engineering at Virginia Tech.

## Affiliation

- **Researcher**: Daniel Miles (Virginia Tech Student)
- **Advisor**: Dr. Rakesh Batra, Mechanical Engineering, Virginia Tech
- **Institution**: Virginia Tech

## Overview

This project investigates the dynamics of 2D mass-spring lattice systems using exponential spring models. The simulations explore:

- Exponential spring force models with configurable decay rates
- 2D lattice systems with nearest-neighbor and diagonal spring connections
- Energy conservation and work calculations
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

- **`exponential-springs.jl`**: Exponential spring model with damping
- **`lattice_w_anim.jl`**: Linear spring model with animation

To run a simulation:

```bash
julia --project=. src/exponential-springs.jl
```

or

```bash
julia --project=. src/lattice_w_anim.jl
```

## Repository Structure

```
Material-Simulations-Research/
├── src/              # Main simulation code
├── docs/             # Comprehensive documentation
├── papers/           # LaTeX papers and research updates
├── scripts/          # Utility and example scripts
├── tests/            # Test suite
├── animations/        # Generated animations
└── data/             # Data files
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
- Analyze energy conservation in lattice systems
- Study wave propagation and response to external forces
- Develop computational tools for material simulation

## Papers

Research papers and progress updates will be added to the `papers/` directory as they are developed.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

For questions or collaboration inquiries, please contact:
- Daniel Miles: danielmiles@vt.edu
- Dr. Rakesh Batra: [Department of Mechanical Engineering, Virginia Tech](https://www.me.vt.edu/)

## Acknowledgments

This research is conducted under the guidance of Dr. Rakesh Batra at Virginia Tech's Department of Mechanical Engineering.

