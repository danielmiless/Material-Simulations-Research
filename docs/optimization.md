# Material Ordering Optimization

## Overview

This document describes the material ordering optimization framework implemented to find the optimal arrangement of materials across columns in the 11×11 lattice system. The goal is to minimize the **peak force** transferred to the backplate by optimizing the order in which different materials are placed across the columns.

This optimization is designed for armor applications (e.g., Kevlar vest simulation), where peak force is the critical metric as it determines injury thresholds and penetration risk.

## Problem Formulation

### Objective
Minimize the peak force transferred to the backplate:

$$\text{minimize} \quad F_{peak} = \max_{t \in [0,T]} |F_{backplate}(t)|$$

where $F_{backplate}(t)$ is the force magnitude on the backplate at time $t$, and $T$ is the simulation end time.

**Why Peak Force?**
- **Time-independent**: Does not depend on simulation duration
- **Injury-relevant**: Peak force determines injury thresholds and penetration risk in armor applications
- **Worst-case focus**: Prevents critical failure scenarios
- **Physically meaningful**: Directly relates to maximum stress on the protected surface

### Decision Variable
The decision variable is a permutation of material indices $[1, 2, \ldots, N]$ where $N=11$ is the number of columns. Each column receives one material type, and the optimizer determines the optimal ordering.

### Constraints
- Each material must be used exactly once (permutation constraint)
- Material properties are fixed and predefined (not optimized)
- The simulation must complete successfully

### Material Properties
Each material is defined by a tuple of properties:
- `k_coupling`: Spring constant for nearest-neighbor springs (N)
- `k_diagonal`: Spring constant for diagonal springs (N)
- `c_damping`: Damping coefficient for nearest-neighbor springs (N·s/m)
- `alpha_coupling`: Exponential decay rate for nearest-neighbor springs (m⁻¹)
- `alpha_diagonal`: Exponential decay rate for diagonal springs (m⁻¹)

By default, materials are created using the scaling pattern from the original simulation, but they can be arbitrarily defined.

## Implementation

### Simulation Modifications

The simulation has been refactored to support material ordering:

1. **Material Property Lookup Functions**: Modified to accept a material ordering array and map materials to columns accordingly.

2. **ODE Right-Hand Side**: Updated to accept material ordering via parameters, allowing different material properties per column.

3. **Force Calculation**: New function `calculate_backplate_force()` computes the peak force transferred to the backplate by:
   - Extracting positions and velocities of rightmost column masses
   - Calculating wall force: $F = k_{wall} \cdot \text{penetration} + c_{wall} \cdot v_x$ (when penetration > 0)
   - Finding the maximum force magnitude over all time steps (time-independent metric)

4. **Optimization Wrapper**: Function `run_simulation_with_material_ordering()` provides a clean interface for optimizers to evaluate different material orderings.

### Optimization Utilities

The `optimization_utils.jl` module provides:

- **Permutation Validation**: Functions to validate and generate permutations
- **Objective Function Wrapper**: Creates objective functions that handle permutation encoding/decoding and track optimization progress
- **Continuous Encoding**: Utilities to convert between permutation and continuous representations (needed for some optimizers)
- **Result Tracking**: State tracking for optimization progress and history

### Optimizers Tested

Four optimization packages are tested, each with different approaches:

#### 1. BlackBoxOptim.jl
- **Algorithm**: Adaptive Differential Evolution (`adaptive_de_rand_1_bin_radiuslimited`)
- **Encoding**: Continuous variables converted to permutations via `argsort`
- **Strengths**: Robust for black-box optimization, handles noisy objectives well
- **File**: `scripts/optimization/optimize_blackbox.jl`

#### 2. Optim.jl
- **Algorithm**: Simulated Annealing (fallback to Nelder-Mead)
- **Encoding**: Continuous variables converted to permutations
- **Strengths**: Well-established optimization library, good for continuous problems
- **File**: `scripts/optimization/optimize_optim.jl`

#### 3. Evolutionary.jl
- **Algorithm**: Genetic Algorithm with custom permutation operators
- **Encoding**: Direct permutation representation
- **Operators**: 
  - Crossover: Order crossover (OX) preserving relative order
  - Mutation: Swap mutation (randomly swap two elements)
- **Strengths**: Naturally handles discrete permutation space
- **File**: `scripts/optimization/optimize_evolutionary.jl`

#### 4. Metaheuristics.jl
- **Algorithms Tested**: Particle Swarm Optimization (PSO), Differential Evolution (DE), Evolutionary Centers Algorithm (ECA), Evolution Strategy (ES)
- **Encoding**: Continuous variables converted to permutations
- **Strengths**: Multiple modern metaheuristic algorithms
- **File**: `scripts/optimization/optimize_metaheuristics.jl`

## Usage

### Running Individual Optimizers

Each optimizer can be run independently:

```julia
# BlackBoxOptim
julia scripts/optimization/optimize_blackbox.jl

# Optim.jl
julia scripts/optimization/optimize_optim.jl

# Evolutionary.jl
julia scripts/optimization/optimize_evolutionary.jl

# Metaheuristics.jl
julia scripts/optimization/optimize_metaheuristics.jl
```

### Comparing All Optimizers

Run the comparison script to test all optimizers:

```julia
julia scripts/optimization/compare_optimizers.jl
```

This will:
1. Run each optimizer with the same configuration
2. Compare results (best force, number of evaluations, computation time)
3. Generate a summary table
4. Save results to `scripts/optimization/comparison_results.txt`

### Custom Material Properties

To use custom material properties, modify the materials array:

```julia
# Define custom materials
custom_materials = [
    (100.0, 50.0, 5.0, 10.0, 10.0),  # Material 1
    (150.0, 75.0, 7.5, 10.0, 10.0),  # Material 2
    # ... define all 11 materials
]

# Run optimization with custom materials
best_perm, best_f, state, result = optimize_material_ordering_blackbox(
    materials=custom_materials,
    max_evaluations=100
)
```

## Questions and Answers

During the planning phase, several questions were asked to clarify the implementation requirements. The answers are documented below:

### Q1: How should we measure "force transferred to the backplate"?

**Answer**: Peak force (maximum force magnitude)

The force is measured as the peak force magnitude over the entire simulation: $\max_{t \in [0,T]} |F_{backplate}(t)|$. This time-independent metric focuses on the worst-case scenario, which is critical for armor applications where injury thresholds are determined by peak force/pressure rather than cumulative effects.

### Q2: Which material parameters should each column be able to have independently?

**Answer**: Material parameters should be set relatively arbitrarily and given as input (not to be optimized)

Each column can have independent values for all material parameters:
- `k_coupling` (spring constant for nearest-neighbor springs)
- `k_diagonal` (spring constant for diagonal springs)
- `c_damping` (damping coefficient)
- `alpha_coupling` (exponential decay rate for nearest-neighbor springs)
- `alpha_diagonal` (exponential decay rate for diagonal springs)

These parameters are predefined and fixed - they are not optimized. Only the ordering of materials (which material goes in which column) is optimized.

### Q3: For the material ordering optimization, should we have a fixed set of N distinct material types (one per column) that can be reordered?

**Answer**: Fixed set of N distinct materials (reorderable)

Yes, we have exactly N=11 distinct material types, one for each column. The optimizer can reorder these materials across the columns, but each material is used exactly once. This ensures a permutation optimization problem.

### Q4: What should the optimizer be allowed to change?

**Answer**: Only the order/permutation of materials (keeping material properties fixed)

The optimizer only changes the ordering of materials. The material properties themselves are fixed and predefined. For example, if Material 1 has properties (100, 50, 5, 10, 10) and Material 2 has properties (150, 75, 7.5, 10, 10), the optimizer can decide whether Material 1 goes in column 1 and Material 2 in column 2, or vice versa, but it cannot change the properties of Material 1 or Material 2.

### Q5: How should the material parameters be specified as input?

**Answer**: Define as constants/arrays in code

Material parameters are defined as arrays of tuples in the code. The default materials are created using the `create_default_material_properties()` function, which uses the existing scaling pattern. These can be modified or replaced with arbitrary material property sets.

### Q6: How should we test multiple optimizers?

**Answer**: Separate scripts for each optimizer

Each optimizer has its own script in `scripts/optimization/`:
- `optimize_blackbox.jl`
- `optimize_optim.jl`
- `optimize_evolutionary.jl`
- `optimize_metaheuristics.jl`

Additionally, a comparison script (`compare_optimizers.jl`) runs all optimizers and compares their performance.

### Q7: Which optimizers should we test?

**Answer**: BlackBoxOptim.jl, Optim.jl, Evolutionary.jl, and Metaheuristics.jl

All four optimization packages are implemented and tested:
- **BlackBoxOptim.jl**: For robust black-box optimization
- **Optim.jl**: For established optimization methods
- **Evolutionary.jl**: For genetic algorithms with permutation support
- **Metaheuristics.jl**: For modern metaheuristic algorithms (PSO, DE, ECA, ES)

### Q8: Where should the optimization code live?

**Answer**: Separate optimization script that calls the simulation

The optimization code is in separate scripts in `scripts/optimization/` that import and call the simulation functions. This keeps the simulation code clean and allows easy testing of different optimizers.

## Results

Results from optimization runs are saved in the `scripts/optimization/` directory:

- Individual optimizer results: `results_<optimizer>.txt`
- Comparison results: `comparison_results.txt`

Each result file contains:
- Best peak force value found (maximum force magnitude)
- Best permutation (material ordering)
- Number of function evaluations
- Elapsed time
- Configuration parameters
- Peak force history (if tracking enabled)

## Future Work

Potential improvements and extensions:

1. **Parallel Evaluation**: Run multiple simulations in parallel to speed up optimization
2. **Surrogate Models**: Use machine learning to approximate the objective function
3. **Multi-objective Optimization**: Optimize for multiple objectives (e.g., force and energy dissipation)
4. **Constraint Handling**: Add constraints on material ordering (e.g., certain materials must be adjacent)
5. **Sensitivity Analysis**: Analyze how sensitive the optimal solution is to material property variations

## References

- BlackBoxOptim.jl: https://github.com/robertfeldt/BlackBoxOptim.jl
- Optim.jl: https://github.com/JuliaNLSolvers/Optim.jl
- Evolutionary.jl: https://github.com/wildart/Evolutionary.jl
- Metaheuristics.jl: https://github.com/jmejia8/Metaheuristics.jl

