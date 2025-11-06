# API Reference

This document provides detailed documentation for the codebase.

## Main Functions

### `spring_force_2d(pos1, pos2, k, alpha)`

Calculate 2D exponential spring force between two masses.

**Parameters**:
- `pos1`: Position vector of mass 1 `[x, y]`
- `pos2`: Position vector of mass 2 `[x, y]`
- `k`: Spring constant (N)
- `alpha`: Exponential decay rate (m⁻¹)

**Returns**: Force vector on mass 1 due to mass 2 `[Fx, Fy]`

**Force Law**: $F = k \frac{\mathbf{r}}{|\mathbf{r}|} (e^{\alpha |\mathbf{r}|} - 1)$

---

### `spring_force_2d(pos1, pos2, k)`

Calculate 2D linear spring force (Hooke's law).

**Parameters**:
- `pos1`: Position vector of mass 1 `[x, y]`
- `pos2`: Position vector of mass 2 `[x, y]`
- `k`: Spring constant (N/m)

**Returns**: Force vector on mass 1 due to mass 2 `[Fx, Fy]`

**Force Law**: $F = k \mathbf{r}$

---

### `damping_force_2d(vel1, vel2, c)`

Calculate 2D damping force between two masses.

**Parameters**:
- `vel1`: Velocity vector of mass 1 `[vx, vy]`
- `vel2`: Velocity vector of mass 2 `[vx, vy]`
- `c`: Damping coefficient (N·s/m)

**Returns**: Damping force vector on mass 1 `[Fx, Fy]`

**Force Law**: $F = -c (\mathbf{v}_1 - \mathbf{v}_2)$

---

### `calculate_force_components(magnitude, angle_degrees)`

Calculate x and y components of force from magnitude and angle.

**Parameters**:
- `magnitude`: Force magnitude (N)
- `angle_degrees`: Angle in degrees (0-360)
  - 0° = +x direction (right)
  - 90° = +y direction (up)
  - 180° = -x direction (left)
  - 270° = -y direction (down)

**Returns**: Tuple `(fx, fy)` - Force components

---

### `lattice_2d_rhs_with_diagonals!(du, u, p, t)`

Right-hand side function for the ODE system. This is the main dynamics function.

**Parameters**:
- `du`: Derivative vector (output)
- `u`: State vector (input)
- `p`: Parameters (unused)
- `t`: Time

**State Vector Structure**:
```
u = [x1, y1, x2, y2, ..., xN, yN, vx1, vy1, vx2, vy2, ..., vxN, vyN]
```

**Modifies**: `du` in-place

---

### `kinetic_energy_2d(vel_matrix)`

Calculate total kinetic energy for 2D motion.

**Parameters**:
- `vel_matrix`: 2×N matrix where each column is `[vx, vy]` for one mass

**Returns**: Total kinetic energy (J)

**Formula**: $T = \frac{1}{2} m \sum |\mathbf{v}_i|^2$

---

### `potential_energy_2d_with_diagonals(pos_matrix)`

Calculate total potential energy for 2D spring system.

**Parameters**:
- `pos_matrix`: 2×N matrix where each column is `[x, y]` for one mass

**Returns**: Total potential energy (J)

**Note**: Implementation depends on spring model (exponential vs. linear)

---

### `work_done_2d_configurable(sol)`

Calculate total work done by external force.

**Parameters**:
- `sol`: ODE solution from DifferentialEquations.jl

**Returns**: Total work (J)

**Formula**: $W = \int \mathbf{F} \cdot d\mathbf{r}$

---

### `create_equilibrium_grid()`

Create equilibrium positions for the lattice.

**Returns**: 2×N matrix of equilibrium positions

**Grid**: Square grid with spacing `GRID_SPACING`

---

### `create_spring_connections_with_diagonals()`

Create list of spring connections for visualization.

**Returns**: Tuple `(nearest_neighbor_connections, diagonal_connections)`

Each connection is a tuple `(k1, k2)` of mass indices.

---

### `run_2d_simulation_with_configurable_force()`

Main simulation function with configurable force angle and application point.

**Returns**: Tuple `(fig, sol, total_work, final_energy)`
- `fig`: GLMakie figure object
- `sol`: ODE solution
- `total_work`: Work done by external force
- `final_energy`: Final total energy

**Side Effects**: 
- Prints simulation parameters and results
- Displays interactive animation window

---

### `save_animation_to_file(filename)`

Save the animation to a video file.

**Parameters**:
- `filename`: Output filename (default: `"lattice_anim.mp4"`)

**Side Effects**: Creates MP4 file in animations directory

---

## Utility Functions

### `lattice_idx(i, j)`

Convert 2D lattice coordinates to 1D index.

**Parameters**:
- `i`: Row index (1 to N)
- `j`: Column index (1 to N)

**Returns**: Mass index (1 to N²)

---

### `lattice_i(k)`

Convert 1D index to row coordinate.

**Parameters**:
- `k`: Mass index (1 to N²)

**Returns**: Row index (1 to N)

---

### `lattice_j(k)`

Convert 1D index to column coordinate.

**Parameters**:
- `k`: Mass index (1 to N²)

**Returns**: Column index (1 to N)

---

## Constants

### Physical Parameters

- `MASS`: Mass per particle (kg)
- `K_COUPLING`: Nearest neighbor spring constant
- `K_DIAGONAL`: Diagonal spring constant
- `ALPHA_COUPLING`: Exponential decay rate for nearest neighbors (m⁻¹)
- `ALPHA_DIAGONAL`: Exponential decay rate for diagonals (m⁻¹)
- `C_DAMPING`: Damping coefficient (N·s/m)

### Force Configuration

- `FORCE_ANGLE_DEGREES`: Force direction angle
- `FORCE_TARGET_ROW`: Target mass row
- `FORCE_TARGET_COL`: Target mass column
- `F_MAG`: Force magnitude (N)
- `F_ACTIVE_TIME`: Force duration (s)

### Lattice Parameters

- `N`: Lattice size (N×N grid)
- `TOTAL_MASSES`: Total number of masses
- `DOF_PER_MASS`: Degrees of freedom per mass (2 for 2D)
- `TOTAL_DOF`: Total degrees of freedom

### Numerical Parameters

- `T_END`: Simulation end time (s)
- `OUTPUT_INTERVAL`: Time between output points (s)
- `REL_TOL`: Relative tolerance for solver
- `ABS_TOL`: Absolute tolerance for solver

### Animation Parameters

- `ANIMATION_SPEED`: Playback speed multiplier
- `NODE_SIZE`: Size of lattice nodes
- `SPRING_WIDTH`: Width of spring lines
- `DIAGONAL_SPRING_WIDTH`: Width of diagonal springs
- `GRID_SPACING`: Spacing between equilibrium positions

## Code Organization

### File Structure

- `src/exponential-springs.jl`: Exponential spring model implementation
- `src/lattice_w_anim.jl`: Linear spring model implementation
- `scripts/run_simulation.jl`: Example usage script
- `tests/runtests.jl`: Test suite

## Dependencies

### Required Packages

- `DifferentialEquations`: ODE solving
- `GLMakie`: Visualization and animation
- `LinearAlgebra`: Linear algebra operations (standard library)
- `Printf`: Formatted output (standard library)

### Standard Library

- `LinearAlgebra`
- `Printf`

## Examples

### Running a Simulation

```julia
using Pkg
Pkg.activate(".")
include("src/exponential-springs.jl")

# Modify constants as needed
const FORCE_ANGLE_DEGREES = 45.0
const F_MAG = 500.0

# Run simulation
fig, sol, work, energy = run_2d_simulation_with_configurable_force()
```

### Saving an Animation

```julia
include("src/exponential-springs.jl")
save_animation_to_file("animations/my_animation.mp4")
```

## Notes

- All functions use SI units (meters, kilograms, seconds, Newtons)
- Angles are specified in degrees for user convenience
- The state vector uses a flattened representation for efficiency
- Visualization requires OpenGL support

