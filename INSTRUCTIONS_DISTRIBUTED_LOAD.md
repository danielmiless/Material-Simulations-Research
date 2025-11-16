# Instructions: Running 11×11 Distributed Load Simulation and Generating Paper

This guide provides step-by-step instructions to:
1. Run the new 11×11 simulation with distributed load
2. Generate all figures for the update paper
3. Compile the LaTeX paper

## Prerequisites

Ensure you have:
- Julia installed (recommended: Julia 1.9+)
- LaTeX distribution installed (for compiling the paper)
- All Julia packages installed (run `Pkg.instantiate()` if needed)

## Step 1: Install/Update Dependencies

From the repository root directory:

```bash
cd /Users/danielmiles/Documents/School/Senior/Research/Material-Simulations-Research
julia --project=. -e "using Pkg; Pkg.instantiate()"
```

This ensures all required packages (DifferentialEquations, GLMakie, CairoMakie, etc.) are installed.

## Step 2: Run the 11×11 Simulation (Interactive)

To see the simulation in action with real-time visualization:

```bash
julia --project=. src/lattice_simulation_11x11.jl
```

**What to expect:**
- The simulation will run and display an interactive animation window
- You'll see the 11×11 lattice with distributed load (11 green force arrows on the left edge)
- The animation shows the lattice response to the distributed load
- Energy plots and statistics are displayed in real-time
- An animation file will be saved to `animations/lattice_anim_11x11_with_backplate.mp4`

**Note:** You can close the animation window after viewing - it won't affect the saved animation file.

## Step 3: Generate Figures for the Paper

To generate all comparison figures needed for the LaTeX paper:

```bash
julia --project=. scripts/generate_distributed_load_figures.jl
```

**What this does:**
1. Generates a snapshot of the 10×10 system with point load (for comparison)
2. Generates a snapshot of the 11×11 system with distributed load
3. Creates a force distribution diagram showing the trapezoidal pattern

**Output files** (saved to `papers/figures/`):
- `point_load_10x10_snapshot.png` - 10×10 system with point load (comparison)
- `distributed_load_11x11_snapshot.png` - 11×11 system with distributed load
- `force_distribution_pattern.png` - Bar chart showing force distribution

**Note:** If the 10×10 snapshot generation fails (e.g., if `FORCE_TARGET_ROW` is not defined), that's okay - the script will continue and generate the other figures.

## Step 4: Compile the LaTeX Paper

Navigate to the papers directory and compile the LaTeX document:

```bash
cd papers/updates
pdflatex Batra_Update_Distributed_Load.tex
pdflatex Batra_Update_Distributed_Load.tex  # Run twice for proper references
```

Or if you prefer using a LaTeX editor:
- Open `papers/updates/Batra_Update_Distributed_Load.tex` in your LaTeX editor
- Compile the document (usually Ctrl+B or Cmd+B)

**Output:** `Batra_Update_Distributed_Load.pdf` will be generated in the same directory.

## Step 5: View the Results

1. **Simulation Animation:** 
   - Location: `animations/lattice_anim_11x11_with_backplate.mp4`
   - Open with any video player

2. **Paper:**
   - Location: `papers/updates/Batra_Update_Distributed_Load.pdf`
   - Open with any PDF viewer

## Quick Reference: All Commands in Sequence

```bash
# 1. Navigate to repository root
cd /Users/danielmiles/Documents/School/Senior/Research/Material-Simulations-Research

# 2. Install dependencies (if needed)
julia --project=. -e "using Pkg; Pkg.instantiate()"

# 3. Run interactive simulation (optional - to see it in action)
julia --project=. src/lattice_simulation_11x11.jl

# 4. Generate figures for paper
julia --project=. scripts/generate_distributed_load_figures.jl

# 5. Compile paper
cd papers/updates
pdflatex Batra_Update_Distributed_Load.tex
pdflatex Batra_Update_Distributed_Load.tex
```

## Troubleshooting

### Issue: "Package not found" errors
**Solution:** Run `julia --project=. -e "using Pkg; Pkg.instantiate()"` to install all required packages.

### Issue: Figure generation fails
**Solution:** 
- Ensure `CairoMakie` is installed: `julia --project=. -e "using Pkg; Pkg.add(\"CairoMakie\")"`
- Check that simulation files exist:
  - `src/lattice_simulation_11x11.jl` (main simulation)
  - `Deprecated Scripts/lattice_simulation_10x10.jl` (for comparison, optional)

### Issue: LaTeX compilation errors
**Solution:**
- Ensure all figure files exist in `papers/figures/`
- Check that the `subcaption` package is installed in your LaTeX distribution
- Run `pdflatex` twice to resolve cross-references

### Issue: Simulation runs but no visualization appears
**Solution:**
- This is normal if running on a headless system
- The animation file will still be saved to `animations/`
- Use the figure generation script to create static snapshots

## Key Features of the 11×11 Simulation

- **121 masses** arranged in an 11×11 lattice
- **Distributed load** on left edge (column 1): F/20, F/10×9, F/10, F/20
- **Global force angle** (default 0° = horizontal right, configurable)
- **Material scaling** across 11 columns
- **Immovable backplate** on the right side
- **Real-time visualization** with 11 force arrows showing distributed load

## File Locations

- **Simulation:** `src/lattice_simulation_11x11.jl`
- **Figure generation script:** `scripts/generate_distributed_load_figures.jl`
- **Paper source:** `papers/updates/Batra_Update_Distributed_Load.tex`
- **Generated figures:** `papers/figures/*.png`
- **Animation:** `animations/lattice_anim_11x11_with_backplate.mp4`

## Next Steps

After generating the paper, you may want to:
1. Review the generated PDF for any formatting issues
2. Adjust figure sizes or captions if needed
3. Run additional simulations with different force angles or parameters
4. Generate additional snapshots at different time points if needed

