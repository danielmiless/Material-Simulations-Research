# Research Methodology

This document describes the research methodology, theoretical background, and computational approach used in this project.

## Research Objectives

The primary objectives of this research are:

1. **Model Development**: Develop and implement exponential spring models for mass-spring lattice systems
2. **Dynamics Analysis**: Investigate the dynamic behavior of 2D lattice systems under various loading conditions
3. **Energy Conservation**: Study energy conservation and work-energy relationships in these systems
4. **Wave Propagation**: Analyze wave propagation and response characteristics

## Theoretical Background

### Mass-Spring Systems

Mass-spring systems are fundamental models in mechanics and materials science. They consist of:
- **Masses**: Point particles with mass $m$
- **Springs**: Connections between masses that exert restoring forces
- **External Forces**: Applied loads that drive the system

### Spring Force Models

#### Linear Springs (Hooke's Law)

The standard linear spring model follows Hooke's law:

$$F = -k \Delta x$$

where:
- $F$ is the spring force
- $k$ is the spring constant
- $\Delta x$ is the displacement from equilibrium

#### Exponential Springs

The exponential spring model uses a force law of the form:

$$F = k \frac{\mathbf{r}}{|\mathbf{r}|} \left( e^{\alpha |\mathbf{r}|} - 1 \right)$$

where:
- $k$ is the spring constant
- $\alpha$ is the exponential decay rate (units: m⁻¹)
- $\mathbf{r}$ is the displacement vector
- $|\mathbf{r}|$ is the magnitude of displacement

The potential energy for an exponential spring is:

$$U = \frac{k}{\alpha} \left( e^{\alpha |\mathbf{r}|} - \alpha |\mathbf{r}| - 1 \right)$$

### Lattice Structure

The simulations use a 2D square lattice with:
- **Nearest Neighbors**: Horizontal and vertical connections
- **Diagonal Connections**: Next-nearest neighbors in an X-pattern
- **Boundary Conditions**: Free boundaries (masses at edges have fewer connections)

### Damping

Viscous damping is applied to nearest-neighbor springs:

$$F_{\text{damping}} = -c (\mathbf{v}_1 - \mathbf{v}_2)$$

where:
- $c$ is the damping coefficient
- $\mathbf{v}_1$ and $\mathbf{v}_2$ are the velocities of connected masses

## Computational Approach

### Numerical Integration

The equations of motion are solved using Julia's `DifferentialEquations.jl` package:

$$\frac{d\mathbf{x}}{dt} = \mathbf{v}$$

$$\frac{d\mathbf{v}}{dt} = \frac{1}{m} \sum \mathbf{F}$$

The solver uses:
- **Algorithm**: Vern9 (9th order Runge-Kutta method)
- **Tolerance**: Relative tolerance $10^{-12}$, absolute tolerance $10^{-14}$
- **Time Stepping**: Adaptive with output at fixed intervals

### State Vector

The system state is represented as a vector:

$$\mathbf{u} = \begin{bmatrix} \mathbf{x}_1 \\ \mathbf{x}_2 \\ \vdots \\ \mathbf{x}_N \\ \mathbf{v}_1 \\ \mathbf{v}_2 \\ \vdots \\ \mathbf{v}_N \end{bmatrix}$$

where $N$ is the number of masses.

### Energy Calculations

The system energy is tracked:

- **Kinetic Energy**: $T = \frac{1}{2} m \sum_i |\mathbf{v}_i|^2$
- **Potential Energy**: Sum over all springs
- **Total Energy**: $E = T + U$
- **Work Input**: Calculated from external force application

## Simulation Parameters

### Default Configuration

- **Lattice Size**: 5×5 = 25 masses
- **Mass**: 1.0 kg per mass
- **Spring Constants**: 
  - Nearest neighbor: 100.0 N (or exponential parameters)
  - Diagonal: 50.0 N (or exponential parameters)
- **Damping**: 5.0 N·s/m (nearest neighbors only)
- **Time Span**: 0 to 5.0 seconds
- **Output Interval**: 0.01 seconds

### Force Configuration

External forces can be configured:
- **Magnitude**: Force strength
- **Angle**: Direction (0° = right, 90° = up, etc.)
- **Application Point**: Which mass receives the force
- **Duration**: How long the force is applied

## Validation

### Energy Conservation

Energy conservation is verified by comparing:
- Total work input from external force
- Final total energy (kinetic + potential)

The error should be small (typically < 1%) for conservative systems.

### Numerical Accuracy

The high-order solver and tight tolerances ensure:
- Accurate trajectory calculations
- Proper energy conservation
- Stable long-time integration

## Future Work

Areas for future investigation:
- Larger lattice sizes
- Different boundary conditions
- Frequency domain analysis
- Comparison with analytical solutions
- Material property relationships

## References

- DifferentialEquations.jl documentation
- Classical mechanics textbooks on mass-spring systems
- Materials science literature on lattice dynamics

