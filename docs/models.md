# Model Descriptions

This document provides detailed descriptions of the mathematical models implemented in this research.

## Exponential Spring Model

### Force Law

The exponential spring model uses a force law that grows exponentially with displacement:

$$F = k \frac{\mathbf{r}}{|\mathbf{r}|} \left( e^{\alpha |\mathbf{r}|} - 1 \right)$$

**Parameters**:
- $k$: Spring constant (units: N)
- $\alpha$: Exponential decay rate (units: m⁻¹)
- $\mathbf{r} = \mathbf{x}_2 - \mathbf{x}_1$: Displacement vector from mass 1 to mass 2

**Characteristics**:
- For small displacements: Approximates linear spring ($F \approx k \alpha |\mathbf{r}|$)
- For large displacements: Force grows exponentially
- Direction: Force points along the displacement vector

### Potential Energy

The potential energy stored in an exponential spring is:

$$U = \frac{k}{\alpha} \left( e^{\alpha |\mathbf{r}|} - \alpha |\mathbf{r}| - 1 \right)$$

This ensures that $F = -\nabla U$.

### Implementation

```julia
function spring_force_2d(pos1, pos2, k, alpha)
    displacement = pos2 - pos1
    distance = norm(displacement)
    
    if distance < 1e-12
        return zeros(2)  # Avoid division by zero
    end
    
    direction = displacement / distance
    force_magnitude = k * (exp(alpha * distance) - 1.0)
    
    return force_magnitude * direction
end
```

## Linear Spring Model

### Force Law (Hooke's Law)

The linear spring model follows Hooke's law:

$$F = -k \mathbf{r}$$

**Parameters**:
- $k$: Spring constant (units: N/m)

**Potential Energy**:

$$U = \frac{1}{2} k |\mathbf{r}|^2$$

### Implementation

```julia
function spring_force_2d(pos1, pos2, k)
    displacement = pos2 - pos1
    return k * displacement  # Hooke's law
end
```

## Damping Model

### Viscous Damping

Damping is applied to nearest-neighbor springs:

$$F_{\text{damping}} = -c (\mathbf{v}_1 - \mathbf{v}_2)$$

**Parameters**:
- $c$: Damping coefficient (units: N·s/m)
- $\mathbf{v}_1, \mathbf{v}_2$: Velocities of connected masses

**Characteristics**:
- Dissipates energy
- Proportional to relative velocity
- Applied only to nearest-neighbor springs in current implementation

### Implementation

```julia
function damping_force_2d(vel1, vel2, c)
    relative_velocity = vel1 - vel2
    return -c * relative_velocity
end
```

## Lattice Structure

### 2D Square Lattice

The system uses a square lattice with:

- **Grid**: $N \times N$ masses arranged in a square
- **Nearest Neighbors**: 
  - Horizontal: $(i, j) \leftrightarrow (i, j \pm 1)$
  - Vertical: $(i, j) \leftrightarrow (i \pm 1, j)$
- **Diagonal Neighbors**:
  - $(i, j) \leftrightarrow (i \pm 1, j \pm 1)$

### Connectivity

For a 5×5 lattice:
- **Interior masses**: 8 connections (4 nearest + 4 diagonal)
- **Edge masses**: 5 connections
- **Corner masses**: 3 connections

### Index Mapping

1D index from 2D coordinates:
```julia
lattice_idx(i, j) = (i-1)*N + j
```

2D coordinates from 1D index:
```julia
lattice_i(k) = 1 + div(k-1, N)
lattice_j(k) = 1 + mod(k-1, N)
```

## External Force Model

### Configurable Force

External forces can be applied with:

- **Magnitude**: $F_{\text{mag}}$ (units: N)
- **Angle**: $\theta$ (degrees: 0° = +x, 90° = +y)
- **Application Point**: Mass at position $(i, j)$
- **Duration**: $t_{\text{active}}$ (seconds)

### Force Components

$$F_x = F_{\text{mag}} \cos(\theta)$$
$$F_y = F_{\text{mag}} \sin(\theta)$$

### Implementation

```julia
function calculate_force_components(magnitude, angle_degrees)
    angle_radians = deg2rad(angle_degrees)
    fx = magnitude * cos(angle_radians)
    fy = magnitude * sin(angle_radians)
    return fx, fy
end
```

## Equations of Motion

### System Dynamics

For each mass $i$:

$$\frac{d\mathbf{x}_i}{dt} = \mathbf{v}_i$$

$$\frac{d\mathbf{v}_i}{dt} = \frac{1}{m_i} \sum_j \mathbf{F}_{ij} + \frac{1}{m_i} \mathbf{F}_{\text{ext}, i}$$

where:
- $\mathbf{F}_{ij}$: Force from spring/damping connecting masses $i$ and $j$
- $\mathbf{F}_{\text{ext}, i}$: External force on mass $i$

### State Vector

The full system state:

$$\mathbf{u} = \begin{bmatrix}
\mathbf{x}_1 \\
\mathbf{x}_2 \\
\vdots \\
\mathbf{x}_N \\
\mathbf{v}_1 \\
\mathbf{v}_2 \\
\vdots \\
\mathbf{v}_N
\end{bmatrix}$$

For a 5×5 lattice with 2D motion: 50 position DOF + 50 velocity DOF = 100 total DOF

## Energy Calculations

### Kinetic Energy

$$T = \frac{1}{2} m \sum_{i=1}^N |\mathbf{v}_i|^2$$

### Potential Energy

Sum over all springs:

$$U = \sum_{\text{springs}} U_{\text{spring}}$$

For exponential springs:
$$U_{\text{spring}} = \frac{k}{\alpha} \left( e^{\alpha |\mathbf{r}|} - \alpha |\mathbf{r}| - 1 \right)$$

For linear springs:
$$U_{\text{spring}} = \frac{1}{2} k |\mathbf{r}|^2$$

### Work Done

Work done by external force:

$$W = \int \mathbf{F}_{\text{ext}} \cdot d\mathbf{r}$$

Discretized:
$$W \approx \sum_n \mathbf{F}_{\text{ext}} \cdot \Delta \mathbf{r}_n$$

## Numerical Parameters

### Solver Settings

- **Algorithm**: Vern9 (9th order Runge-Kutta)
- **Relative Tolerance**: $10^{-12}$
- **Absolute Tolerance**: $10^{-14}$
- **Dense Output**: Disabled (saves memory)

### Time Stepping

- **Time Span**: $[0, T_{\text{end}}]$ (typically 5.0 seconds)
- **Output Interval**: $\Delta t_{\text{out}}$ (typically 0.01 seconds)
- **Adaptive**: Solver chooses internal time steps

## Model Comparison

| Feature | Linear Spring | Exponential Spring |
|---------|--------------|-------------------|
| Force Law | $F = -k r$ | $F = k (e^{\alpha r} - 1)$ |
| Small $r$ | Linear | Approximately linear |
| Large $r$ | Linear | Exponential growth |
| Potential | $\frac{1}{2} k r^2$ | $\frac{k}{\alpha}(e^{\alpha r} - \alpha r - 1)$ |
| Damping | Optional | Optional |

## Future Model Extensions

Potential additions:
- Nonlinear damping models
- Anisotropic springs
- Temperature effects
- Plastic deformation
- Fracture mechanics

