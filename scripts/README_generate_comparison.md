# Generate Comparison Data Script

This script generates comparison data and figures for the LaTeX update document (`Batra_Update_11_06_2025.tex`).

## What It Does

1. **Runs both simulations**:
   - Without backplate (using `lattice_simulation.jl`)
   - With backplate (using `lattice_simulation_with_backplate.jl`)

2. **Extracts energy statistics**:
   - Total work input
   - Final total energy
   - Energy dissipated
   - Dissipation percentage

3. **Generates figures**:
   - `energy_comparison_damping.png` - Energy evolution comparison
   - `backplate_initial_setup.png` - Initial configuration with backplate
   - `no_backplate_snapshot.png` - Snapshot without backplate
   - `with_backplate_snapshot.png` - Snapshot with backplate

4. **Outputs LaTeX table data**:
   - Prints formatted table data ready to copy into the LaTeX document

## Usage

From the repository root:

```bash
julia scripts/generate_comparison_data.jl
```

**Note**: The script will display animation windows for both simulations. You can close these windows - they won't affect the data extraction or figure generation.

## Output

- **Figures**: Saved to `papers/figures/`
- **Table data**: Printed to console in LaTeX format

## Requirements

- Both simulation files must be present in `src/`:
  - `lattice_simulation.jl`
  - `lattice_simulation_with_backplate.jl`
- Julia packages: DifferentialEquations, GLMakie, LinearAlgebra, Printf
- The script will activate the project environment automatically

## Troubleshooting

If you encounter errors:
1. Ensure both simulation files exist and are complete
2. Check that all required Julia packages are installed
3. Make sure the simulation files have the same parameters (force magnitude, etc.) except for the backplate

