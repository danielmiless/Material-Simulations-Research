"""
Exponential Spring Mass-Spring Lattice Simulation: 11×11 System with Immovable Backplate, Material Scaling, and Distributed Load

This script implements a 2D mass-spring lattice system (11×11) with exponential spring models,
an immovable backplate to the right of the lattice, column-based material scaling, and a distributed load
applied along the left edge (x=0, column 1).
The system includes:
- Exponential spring force law: F = k * (r/|r|) * (exp(alpha*|r|) - 1)
- Viscous damping on nearest-neighbor springs (energy dissipation)
- Nearest-neighbor and diagonal spring connections
- Immovable backplate on the right side (fixed boundary condition)
- Column-based material property scaling (columns 2-11 scaled by MATERIAL_MULTIPLIER)
- Distributed external force application on left edge (trapezoidal distribution: F/20, F/10×9, F/20)
- Real-time visualization and animation
- Energy tracking (with dissipation)

Research conducted in collaboration with Dr. Romesh Batra, Distinguished Professor and Clifton C. Garvin Professor, Department of Mechanical Engineering, Virginia Tech.

Author: Daniel Miles
Institution: Virginia Tech
"""

using DifferentialEquations
using LinearAlgebra
using Printf
using GLMakie


########################################################################
#  PHYSICAL PARAMETERS
########################################################################
const MASS        = 1.0          # kg for every mass
const K_COUPLING  = 100.0        # N for nearest neighbor (horizontal/vertical) springs
const K_DIAGONAL  = 50.0         # N for diagonal springs (typically weaker)
const ALPHA_COUPLING = 10.0      # exponential decay rate for nearest neighbor springs (m⁻¹)
const ALPHA_DIAGONAL = 10.0      # exponential decay rate for diagonal springs (m⁻¹)
const C_DAMPING   = 5.0          # N·s/m damping coefficient for nearest-neighbor springs
const WALL_STIFFNESS = 10000.0   # N/m for wall repulsion force (very stiff to prevent penetration)
const WALL_DAMPING = 10.0        # N·s/m damping coefficient for wall collisions
const F_ACTIVE_TIME = 0.15       # duration of the driving pulse (s) - increased for larger system
const F_MAG       = 1200.0        # N, total magnitude of distributed external force
const BACKPLATE_DISTANCE = 1.0   # Distance from rightmost column to backplate (m)
const MATERIAL_MULTIPLIER = 2/3   # Material property scaling factor for columns (1.5 for first case, 2/3 for second case)


########################################################################
#  FORCE CONFIGURATION PARAMETERS - EASILY CHANGEABLE
########################################################################
const FORCE_ANGLE_DEGREES = 0.0    # Global angle in degrees (0-360) for all distributed forces: 0=right, 90=up, 180=left, 270=down


########################################################################
#  LATTICE SIZE  
########################################################################
const N = 11                     # ⇒ 11×11 = 121 masses
const TOTAL_MASSES = N * N       # 121 masses total
const DOF_PER_MASS = 2          # x and y components
const TOTAL_DOF = TOTAL_MASSES * DOF_PER_MASS  # 242 position DOF


########################################################################
#  NUMERICAL PARAMETERS
########################################################################
const T_END          = 8.0       # s - increased for larger system to see full material gradient effects
const OUTPUT_INTERVAL= 0.01      # s
const REL_TOL        = 1e-12
const ABS_TOL        = 1e-14


########################################################################
#  ANIMATION PARAMETERS
########################################################################
const ANIMATION_SPEED = 1.0      # Animation playback speed multiplier
const NODE_SIZE = 20             # Size of lattice nodes
const SPRING_WIDTH = 3           # Width of spring connections
const DIAGONAL_SPRING_WIDTH = 2  # Width of diagonal spring connections
const GRID_SPACING = 1.0         # Spacing between equilibrium positions


########################################################################
#  UTILITY: 2-D ↔ 1-D index conversion and state vector indexing
########################################################################
@inline lattice_idx(i,j) = (i-1)*N + j                    # (i,j) → mass index 1…121
@inline lattice_i(k) = 1 + div(k-1,N)                     # mass index → row
@inline lattice_j(k) = 1 + mod(k-1,N)                     # mass index → column


########################################################################
#  MATERIAL PROPERTY LOOKUP FUNCTIONS (COLUMN-BASED SCALING)
########################################################################
function get_column_k_coupling(j, material_multiplier)
    """
    Get spring constant for nearest-neighbor springs in column j.
    Column 1 uses base value, columns 2-11 use scaled values.
    """
    if j == 1
        return K_COUPLING
    else
        return K_COUPLING * (material_multiplier^(j-1))
    end
end

function get_column_k_diagonal(j, material_multiplier)
    """
    Get spring constant for diagonal springs in column j.
    Column 1 uses base value, columns 2-11 use scaled values.
    """
    if j == 1
        return K_DIAGONAL
    else
        return K_DIAGONAL * (material_multiplier^(j-1))
    end
end

function get_column_c_damping(j, material_multiplier)
    """
    Get damping coefficient for nearest-neighbor springs in column j.
    Column 1 uses base value, columns 2-11 use scaled values.
    """
    if j == 1
        return C_DAMPING
    else
        return C_DAMPING * (material_multiplier^(j-1))
    end
end


########################################################################
#  DISTRIBUTED LOAD FUNCTION
########################################################################
function get_distributed_force_magnitude(row_index)
    """
    Get force magnitude for distributed load at a given row.
    Distribution pattern: F/20 at edges (rows 1 and 11), F/10 for middle rows (2-10).
    
    Args:
        row_index: Row index (1 to N)
    
    Returns:
        Force magnitude for this row
    """
    if row_index == 1 || row_index == 11
        return F_MAG / 20.0
    else
        return F_MAG / 10.0
    end
end


########################################################################
#  BACKPLATE POSITIONS (IMMOVABLE)
########################################################################
function create_backplate_positions()
    """
    Create fixed positions for the immovable backplate.
    Returns a 2×N matrix where each column is [x, y] for one backplate point.
    Backplate is positioned to the right of the rightmost column.
    """
    backplate_positions = zeros(2, N)
    
    for i in 1:N
        # x position: rightmost column is at (N-1)*GRID_SPACING, backplate is BACKPLATE_DISTANCE further
        backplate_positions[1, i] = (N - 1) * GRID_SPACING + BACKPLATE_DISTANCE
        # y position: matches the row position
        backplate_positions[2, i] = (i - 1) * GRID_SPACING
    end
    
    return backplate_positions
end


########################################################################
#  FORCE VECTOR CALCULATION FROM ANGLE AND MAGNITUDE
########################################################################
function calculate_force_components(magnitude, angle_degrees)
    """
    Calculate x and y components of force from magnitude and angle.
    
    Args:
        magnitude: Force magnitude (N)
        angle_degrees: Angle in degrees (0-360)
                      0° = +x direction (right)
                      90° = +y direction (up)  
                      180° = -x direction (left)
                      270° = -y direction (down)
    
    Returns:
        (fx, fy): Force components in x and y directions
    """
    angle_radians = deg2rad(angle_degrees)
    fx = magnitude * cos(angle_radians)
    fy = magnitude * sin(angle_radians)
    return fx, fy
end


########################################################################
#  2D EXPONENTIAL SPRING FORCE CALCULATION
########################################################################
function spring_force_2d(pos1, pos2, k, alpha)
    """
    Calculate 2D exponential spring force between two masses.
    Returns force on mass 1 due to mass 2.
    
    For exponential spring: F = k * (r/|r|) * (exp(alpha*|r|) - 1)
    where r = pos2 - pos1 is the displacement vector
    """
    displacement = pos2 - pos1  # vector from mass1 to mass2
    distance = norm(displacement)
    
    if distance < 1e-12
        return zeros(2)
    end
    
    direction = displacement / distance
    force_magnitude = k * (exp(alpha * distance) - 1.0)
    
    return force_magnitude * direction
end


########################################################################
#  2D DAMPING FORCE CALCULATION
########################################################################
function damping_force_2d(vel1, vel2, c)
    """
    Calculate 2D viscous damping force between two masses.
    Returns damping force on mass 1 due to relative velocity with mass 2.
    
    F_damping = -c * (v1 - v2)
    where c is the damping coefficient and v1, v2 are velocities
    """
    relative_velocity = vel1 - vel2
    return -c * relative_velocity
end


########################################################################
#  ODE RIGHT-HAND SIDE FOR 2D MOTION WITH DIAGONAL SPRINGS AND BACKPLATE
########################################################################
function lattice_2d_rhs_with_diagonals_and_backplate!(du, u, p, t)
    # Extract positions and velocities from state vector
    pos = reshape(view(u, 1:TOTAL_DOF), 2, TOTAL_MASSES)
    vel = reshape(view(u, TOTAL_DOF+1:2*TOTAL_DOF), 2, TOTAL_MASSES)
    
    dpos = reshape(view(du, 1:TOTAL_DOF), 2, TOTAL_MASSES)
    dvel = reshape(view(du, TOTAL_DOF+1:2*TOTAL_DOF), 2, TOTAL_MASSES)

    # --- kinematics: dx/dt = v ---
    dpos .= vel

    # --- initialize forces to zero ---
    dvel .= 0.0

    # --- calculate spring forces between neighboring masses ---
    # Optimized: process each mass once, checking all neighbors
    for k in 1:TOTAL_MASSES
        i, j = lattice_i(k), lattice_j(k)
        pos_k = pos[:, k]
        vel_k = vel[:, k]
        
        # Helper function to add spring force from neighbor
        function add_spring_force!(neighbor_i, neighbor_j, k_spring, alpha_spring, add_damping=false, damping_coeff=C_DAMPING)
            if 1 <= neighbor_i <= N && 1 <= neighbor_j <= N
                neighbor_idx = lattice_idx(neighbor_i, neighbor_j)
                pos_neighbor = pos[:, neighbor_idx]
                vel_neighbor = vel[:, neighbor_idx]
                
                # Spring force
                force = spring_force_2d(pos_k, pos_neighbor, k_spring, alpha_spring)
                dvel[:, k] += force
                
                # Damping force (only for nearest neighbors)
                if add_damping
                    damping = damping_force_2d(vel_k, vel_neighbor, damping_coeff)
                    dvel[:, k] += damping
                end
            end
        end
        
        # NEAREST NEIGHBOR SPRINGS (horizontal and vertical) WITH DAMPING
        # Horizontal springs use base values
        add_spring_force!(i, j-1, K_COUPLING, ALPHA_COUPLING, true, C_DAMPING)  # left
        add_spring_force!(i, j+1, K_COUPLING, ALPHA_COUPLING, true, C_DAMPING)  # right
        
        # Vertical springs use column-based material properties
        k_col = get_column_k_coupling(j, MATERIAL_MULTIPLIER)
        c_col = get_column_c_damping(j, MATERIAL_MULTIPLIER)
        add_spring_force!(i-1, j, k_col, ALPHA_COUPLING, true, c_col)  # up
        add_spring_force!(i+1, j, k_col, ALPHA_COUPLING, true, c_col)  # down
        
        # DIAGONAL SPRINGS (next-nearest neighbors in X pattern) WITHOUT DAMPING
        add_spring_force!(i-1, j-1, K_DIAGONAL, ALPHA_DIAGONAL, false)  # upper-left
        add_spring_force!(i-1, j+1, K_DIAGONAL, ALPHA_DIAGONAL, false)  # upper-right
        add_spring_force!(i+1, j-1, K_DIAGONAL, ALPHA_DIAGONAL, false)  # lower-left
        add_spring_force!(i+1, j+1, K_DIAGONAL, ALPHA_DIAGONAL, false)  # lower-right
        
        # WALL CONSTRAINT (rightmost column cannot pass through backplate)
        if j == N  # Rightmost column
            # Equilibrium x position of rightmost column: (N-1)*GRID_SPACING
            # Backplate x position: (N-1)*GRID_SPACING + BACKPLATE_DISTANCE
            # Current x position: equilibrium_x + displacement_x
            equilibrium_x = (N - 1) * GRID_SPACING
            backplate_x = equilibrium_x + BACKPLATE_DISTANCE
            current_x = equilibrium_x + pos_k[1]  # absolute x position
            
            # If mass tries to move past the backplate, apply repulsive force
            if current_x > backplate_x
                penetration = current_x - backplate_x
                # Very stiff repulsive force (like a wall)
                wall_force_x = -WALL_STIFFNESS * penetration
                dvel[1, k] += wall_force_x
                
                # Damping when colliding with wall (only x-component)
                if vel_k[1] > 0  # Moving toward wall
                    dvel[1, k] -= WALL_DAMPING * vel_k[1]
                end
            end
        end
    end

    # --- DISTRIBUTED EXTERNAL DRIVING FORCE ON LEFT EDGE (x=0, column 1) ---
    if t <= F_ACTIVE_TIME
        # Calculate unit direction vector from global angle
        fx_unit, fy_unit = calculate_force_components(1.0, FORCE_ANGLE_DEGREES)
        
        # Apply distributed load to all masses in column 1 (left edge)
        for i in 1:N  # All rows in column 1
            force_mag = get_distributed_force_magnitude(i)
            target_idx = lattice_idx(i, 1)  # Column 1, row i
            dvel[1, target_idx] += force_mag * fx_unit  # x-component
            dvel[2, target_idx] += force_mag * fy_unit  # y-component
        end
    end

    # --- divide by mass to obtain accelerations ---
    dvel ./= MASS
    
    return nothing
end


########################################################################
#  ENERGY & WORK CALCULATIONS FOR 2D MOTION WITH DIAGONAL SPRINGS AND BACKPLATE
########################################################################
function kinetic_energy_2d(vel_matrix)
    """
    Calculate total kinetic energy for 2D motion.
    vel_matrix is 2×N matrix where each column is [vx, vy] for one mass.
    """
    return 0.5 * MASS * sum(vel_matrix.^2)
end


function potential_energy_2d_with_diagonals_and_backplate(pos_matrix)
    """
    Calculate total potential energy for 2D exponential spring system with diagonal springs and backplate.
    pos_matrix is 2×N matrix where each column is [x, y] for one mass.
    
    For exponential spring: U = (k/alpha) * (exp(alpha*|r|) - alpha*|r| - 1)
    """
    pe = 0.0
    
    # Helper function to calculate and add potential energy from spring
    function add_spring_pe(pos1, pos2, k_spring, alpha_spring)
        displacement = pos2 - pos1
        distance = norm(displacement)
        return (k_spring / alpha_spring) * (exp(alpha_spring * distance) - alpha_spring * distance - 1.0)
    end
    
    # Horizontal springs (nearest neighbors)
    for i in 1:N, j in 1:N-1
        k1 = lattice_idx(i, j)
        k2 = lattice_idx(i, j+1)
        pe += add_spring_pe(pos_matrix[:, k1], pos_matrix[:, k2], K_COUPLING, ALPHA_COUPLING)
    end
    
    # Vertical springs (nearest neighbors)
    for i in 1:N-1, j in 1:N
        k1 = lattice_idx(i, j)
        k2 = lattice_idx(i+1, j)
        pe += add_spring_pe(pos_matrix[:, k1], pos_matrix[:, k2], K_COUPLING, ALPHA_COUPLING)
    end
    
    # Diagonal springs (next-nearest neighbors)
    # Upper-left to lower-right diagonals
    for i in 1:N-1, j in 1:N-1
        k1 = lattice_idx(i, j)
        k2 = lattice_idx(i+1, j+1)
        pe += add_spring_pe(pos_matrix[:, k1], pos_matrix[:, k2], K_DIAGONAL, ALPHA_DIAGONAL)
    end
    
    # Upper-right to lower-left diagonals
    for i in 1:N-1, j in 2:N
        k1 = lattice_idx(i, j)
        k2 = lattice_idx(i+1, j-1)
        pe += add_spring_pe(pos_matrix[:, k1], pos_matrix[:, k2], K_DIAGONAL, ALPHA_DIAGONAL)
    end
    
    # Wall potential energy (penalty for penetration)
    # This is not a spring, but a penalty term for masses that penetrate the wall
    equilibrium_x = (N - 1) * GRID_SPACING
    backplate_x = equilibrium_x + BACKPLATE_DISTANCE
    for i in 1:N
        k = lattice_idx(i, N)  # Rightmost column
        current_x = equilibrium_x + pos_matrix[1, k]  # absolute x position
        if current_x > backplate_x
            penetration = current_x - backplate_x
            # Penalty potential energy (like a very stiff spring)
            pe += 0.5 * WALL_STIFFNESS * penetration^2
        end
    end
    
    return pe
end


function work_done_2d_distributed(sol)
    """
    Calculate total work done by distributed external forces.
    Work = sum over all forces of F · dr for each driven mass.
    """
    # Calculate unit direction vector from global angle
    fx_unit, fy_unit = calculate_force_components(1.0, FORCE_ANGLE_DEGREES)
    
    w = 0.0
    for i in 1:N
        force_mag = get_distributed_force_magnitude(i)
        target_idx = lattice_idx(i, 1)  # Column 1, row i
        
        for n in 2:length(sol.t)
            if sol.t[n-1] >= F_ACTIVE_TIME
                break
            end
            
            # Get positions at current and previous time (using views for efficiency)
            pos_prev = view(sol.u[n-1], (target_idx-1)*2+1:target_idx*2)
            pos_curr = view(sol.u[n], (target_idx-1)*2+1:target_idx*2)
            
            # Calculate displacement and work
            displacement = pos_curr .- pos_prev
            w += force_mag * (fx_unit * displacement[1] + fy_unit * displacement[2])
        end
    end
    
    return w
end


########################################################################
#  EQUILIBRIUM GRID POSITIONS
########################################################################
function create_equilibrium_grid()
    """
    Create equilibrium positions for the N×N lattice.
    """
    equilibrium_positions = zeros(2, TOTAL_MASSES)
    
    for i in 1:N, j in 1:N
        k = lattice_idx(i, j)
        equilibrium_positions[1, k] = (j - 1) * GRID_SPACING  # x coordinate
        equilibrium_positions[2, k] = (i - 1) * GRID_SPACING  # y coordinate
    end
    
    return equilibrium_positions
end


########################################################################
#  SPRING CONNECTIVITY INCLUDING DIAGONALS AND BACKPLATE
########################################################################
function create_spring_connections_with_diagonals()
    """
    Create list of spring connections for visualization, including diagonal springs.
    Returns separate arrays for nearest-neighbor and diagonal connections.
    Note: Backplate is a wall constraint, not connected by springs.
    """
    nearest_neighbor_connections = Tuple{Int, Int}[]
    diagonal_connections = Tuple{Int, Int}[]
    
    # Horizontal springs (nearest neighbors)
    for i in 1:N, j in 1:N-1
        push!(nearest_neighbor_connections, (lattice_idx(i, j), lattice_idx(i, j+1)))
    end
    
    # Vertical springs (nearest neighbors)
    for i in 1:N-1, j in 1:N
        push!(nearest_neighbor_connections, (lattice_idx(i, j), lattice_idx(i+1, j)))
    end
    
    # Diagonal springs (next-nearest neighbors)
    # Upper-left to lower-right diagonals
    for i in 1:N-1, j in 1:N-1
        push!(diagonal_connections, (lattice_idx(i, j), lattice_idx(i+1, j+1)))
    end
    
    # Upper-right to lower-left diagonals
    for i in 1:N-1, j in 2:N
        push!(diagonal_connections, (lattice_idx(i, j), lattice_idx(i+1, j-1)))
    end
    
    return nearest_neighbor_connections, diagonal_connections
end


########################################################################
#  FORCE VISUALIZATION HELPER
########################################################################
function create_distributed_force_arrows(equilibrium_grid)
    """
    Create points for visualizing all distributed forces as arrows.
    Returns array of arrow point pairs for all forces in column 1.
    """
    arrows = Point2f[]
    fx_unit, fy_unit = calculate_force_components(1.0, FORCE_ANGLE_DEGREES)
    
    # Scale for visualization (make arrows visible but not too long)
    scale_factor = 0.3
    
    for i in 1:N
        target_idx = lattice_idx(i, 1)  # Column 1, row i
        target_pos = equilibrium_grid[:, target_idx]
        force_mag = get_distributed_force_magnitude(i)
        
        # Scale arrow length proportional to force magnitude
        arrow_length = scale_factor * force_mag / F_MAG
        force_end = target_pos + arrow_length * [fx_unit, fy_unit]
        
        push!(arrows, Point2f(target_pos))
        push!(arrows, Point2f(force_end))
    end
    
    return arrows
end


########################################################################
#  MAIN SIMULATION AND ANIMATION FUNCTION
########################################################################
function run_2d_simulation_with_configurable_force_and_backplate()
    """
    Main simulation function with distributed load on left edge and immovable backplate.
    Returns: (figure, solution, total_work, final_energy)
    """
    println("11×11 Lattice Mass-Spring System with Distributed Load, Immovable Backplate, and Material Scaling")
    println("=" ^ 80)
    println("Physical Parameters:")
    println("  Total masses: $(TOTAL_MASSES)")
    println("  DOF per mass: $(DOF_PER_MASS) (x, y)")
    println("  Total DOF: $(TOTAL_DOF)")
    println("  Nearest neighbor spring constant: $(K_COUPLING) N")
    println("  Nearest neighbor alpha: $(ALPHA_COUPLING) m⁻¹")
    println("  Nearest neighbor damping: $(C_DAMPING) N·s/m")
    println("  Diagonal spring constant: $(K_DIAGONAL) N")
    println("  Diagonal alpha: $(ALPHA_DIAGONAL) m⁻¹")
    println("  Wall stiffness: $(WALL_STIFFNESS) N/m")
    println("  Wall damping: $(WALL_DAMPING) N·s/m")
    println("  Backplate distance: $(BACKPLATE_DISTANCE) m")
    println("  Mass: $(MASS) kg")
    println("")
    println("Spring Network:")
    println("  Nearest neighbors: horizontal & vertical connections (with damping)")
    println("  Diagonal springs: X-pattern connections in each square (no damping)")
    println("  Wall constraint: rightmost column cannot pass through immovable backplate")
    println("  Total connectivity: 8 neighbors per interior mass + wall constraint for rightmost column")
    println("  Spring type: Exponential")
    println("")
    println("Material Scaling (Column-Based):")
    println("  Material multiplier: $(MATERIAL_MULTIPLIER)×")
    println("  Column 1: base values (K=$(K_COUPLING) N, c=$(C_DAMPING) N·s/m)")
    # Show first few and last few columns for brevity
    for j in 2:4
        k_col = get_column_k_coupling(j, MATERIAL_MULTIPLIER)
        c_col = get_column_c_damping(j, MATERIAL_MULTIPLIER)
        println("  Column $j: $(round(MATERIAL_MULTIPLIER^(j-1), digits=3))× base (K=$(round(k_col, digits=1)) N, c=$(round(c_col, digits=2)) N·s/m)")
    end
    println("  ... (columns 5-10 with intermediate scaling) ...")
    for j in N:N
        k_col = get_column_k_coupling(j, MATERIAL_MULTIPLIER)
        c_col = get_column_c_damping(j, MATERIAL_MULTIPLIER)
        println("  Column $j: $(round(MATERIAL_MULTIPLIER^(j-1), digits=3))× base (K=$(round(k_col, digits=1)) N, c=$(round(c_col, digits=2)) N·s/m)")
    end
    println("  Visualization: Columns color-coded from light (column 1) to dark (column $N) to show gradient")
    println("")
    println("Distributed External Force (Left Edge, Column 1):")
    println("  Total magnitude: $(F_MAG) N")
    println("  Global angle: $(FORCE_ANGLE_DEGREES)° (applies to all forces, 0°=right, 90°=up, 180°=left, 270°=down)")
    println("  Distribution pattern:")
    println("    Row 1: $(F_MAG/20.0) N (F/20)")
    println("    Rows 2-10: $(F_MAG/10.0) N each (F/10)")
    println("    Row 11: $(F_MAG/20.0) N (F/20)")
    println("  Duration: $(F_ACTIVE_TIME) s")
    
    fx_unit, fy_unit = calculate_force_components(1.0, FORCE_ANGLE_DEGREES)
    println("  Force direction components: fx = $(round(fx_unit, digits=3)), fy = $(round(fy_unit, digits=3))")
    println("")

    # Initial state: all masses at rest at equilibrium positions
    u0 = zeros(2 * TOTAL_DOF)
    tspan = (0.0, T_END)

    # Create and solve ODE problem
    prob = ODEProblem(lattice_2d_rhs_with_diagonals_and_backplate!, u0, tspan)
    sol = solve(prob, Vern9();
                reltol = REL_TOL, 
                abstol = ABS_TOL,
                saveat = OUTPUT_INTERVAL,
                dense = false)

    if sol.retcode != :Success
        error("Solver failed: $(sol.retcode)")
    end

    println("Solver completed successfully!")
    println("Final time: $(sol.t[end]) s")
    println("Time steps: $(length(sol.t))")
    println("")

    # Calculate total work done by external forces
    total_work = work_done_2d_distributed(sol)
    
    # Create equilibrium grid, spring connections, and backplate positions
    equilibrium_grid = create_equilibrium_grid()
    nearest_neighbor_connections, diagonal_connections = create_spring_connections_with_diagonals()
    backplate_positions = create_backplate_positions()
    
    # Set up GLMakie figure and layout
    GLMakie.activate!()
    fig = Figure(size = (1400, 900), fontsize = 14)
    
    # Main animation axis
    ax = Axis(fig[1, 1:2], 
              title = "11×11 Mass-Spring Lattice with Distributed Load, Immovable Backplate, and Material Scaling",
              xlabel = "X Position",
              ylabel = "Y Position",
              aspect = DataAspect())
    
    # Energy tracking axis
    energy_ax = Axis(fig[2, 1],
                     title = "Energy Conservation",
                     xlabel = "Time (s)",
                     ylabel = "Energy")
    
    # Statistics display
    stats_ax = Axis(fig[2, 2],
                    title = "System Statistics",
                    xlabel = "",
                    ylabel = "")
    hidedecorations!(stats_ax)
    
    # Initialize observables for animation
    current_frame = Observable(1)
    
    # Extract position data for all time steps (pre-compute for efficiency)
    all_positions = [reshape(view(sol.u[i], 1:TOTAL_DOF), 2, TOTAL_MASSES) for i in 1:length(sol.t)]
    
    # Create observables for current positions
    current_positions = @lift(all_positions[$current_frame] .+ equilibrium_grid)
    
    # Plot nearest neighbor springs as lines
    for (k1, k2) in nearest_neighbor_connections
        spring_points = @lift(Point2f[($current_positions)[:, k1], ($current_positions)[:, k2]])
        lines!(ax, spring_points, color = :blue, linewidth = SPRING_WIDTH, label = "Nearest Neighbor")
    end
    
    # Plot diagonal springs as lines
    for (k1, k2) in diagonal_connections
        spring_points = @lift(Point2f[($current_positions)[:, k1], ($current_positions)[:, k2]])
        lines!(ax, spring_points, color = :red, linewidth = DIAGONAL_SPRING_WIDTH, label = "Diagonal")
    end
    
    # Plot backplate as a vertical line (wall, no springs)
    backplate_line_points = Point2f[backplate_positions[:, i] for i in 1:N]
    lines!(ax, backplate_line_points, color = :gray, linewidth = 8, label = "Backplate")
    
    # Plot masses as scatter points
    mass_points = @lift(Point2f[($current_positions)[:, k] for k in 1:TOTAL_MASSES])
    
    # Color-code masses by column to showcase material gradient
    # Use a color gradient: lighter (column 1) to darker (column N) to show increasing/decreasing stiffness
    colors = Vector{Symbol}(undef, TOTAL_MASSES)
    
    # Color scheme based on column index to visualize material gradient (11 columns)
    column_colors = [:lightblue, :lightcyan, :cyan, :blue, :mediumblue, :darkblue, :navy, :midnightblue, :darkslateblue, :darkblue, :black]
    for k in 1:TOTAL_MASSES
        i, j = lattice_i(k), lattice_j(k)
        # Map column index to color (1-based, so j=1..11)
        color_idx = min(j, length(column_colors))
        colors[k] = column_colors[color_idx]
    end
    
    scatter!(ax, mass_points, markersize = NODE_SIZE, color = colors, strokewidth = 2, strokecolor = :white)
    
    # Add distributed force arrows visualization (visible only during active time)
    force_active = @lift(sol.t[$current_frame] <= F_ACTIVE_TIME)
    force_arrows = @lift(begin
        if $force_active
            create_distributed_force_arrows(equilibrium_grid)
        else
            Point2f[]
        end
    end)
    
    # Plot force arrows (each arrow is two consecutive points)
    for i in 1:N
        arrow_start = @lift(begin
            if $force_active && length($force_arrows) >= 2*i
                [$force_arrows[2*i-1]]
            else
                Point2f[]
            end
        end)
        arrow_end = @lift(begin
            if $force_active && length($force_arrows) >= 2*i
                [$force_arrows[2*i]]
            else
                Point2f[]
            end
        end)
        
        # Draw arrow line
        arrow_line = @lift(begin
            if $force_active && length($force_arrows) >= 2*i
                [$force_arrows[2*i-1], $force_arrows[2*i]]
            else
                Point2f[]
            end
        end)
        lines!(ax, arrow_line, color = :green, linewidth = 3)
        
        # Draw arrow head
        scatter!(ax, arrow_end, markersize = 8, color = :green, marker = :circle)
    end
    
    # Add legend
    axislegend(ax, position = :rt)
    
    # Set up axis limits with some padding (extend to show backplate)
    grid_min = minimum(equilibrium_grid) - 0.5
    grid_max_x = maximum(backplate_positions[1, :]) + 0.5
    grid_max_y = maximum(backplate_positions[2, :]) + 0.5
    limits!(ax, grid_min, grid_max_x, grid_min, grid_max_y)
    
    # Energy conservation plot (pre-compute for efficiency)
    time_values = sol.t
    energy_values = [let
        pos = reshape(view(sol.u[i], 1:TOTAL_DOF), 2, TOTAL_MASSES)
        vel = reshape(view(sol.u[i], TOTAL_DOF+1:2*TOTAL_DOF), 2, TOTAL_MASSES)
        kinetic_energy_2d(vel) + potential_energy_2d_with_diagonals_and_backplate(pos)
    end for i in 1:length(sol.t)]
    
    work_line = fill(total_work, length(time_values))
    
    lines!(energy_ax, time_values, energy_values, color = :blue, label = "Total Energy")
    lines!(energy_ax, time_values, work_line, color = :red, linestyle = :dash, label = "Work Input")
    axislegend(energy_ax)
    
    # Current time indicator
    current_time = @lift(sol.t[$current_frame])
    vlines!(energy_ax, current_time, color = :green, linewidth = 2)
    
    # Statistics text
    stats_text = @lift begin
        t = sol.t[$current_frame]
        pos = all_positions[$current_frame]
        vel = reshape(view(sol.u[$current_frame], TOTAL_DOF+1:2*TOTAL_DOF), 2, TOTAL_MASSES)
        
        ke = kinetic_energy_2d(vel)
        pe = potential_energy_2d_with_diagonals_and_backplate(pos)
        total_e = ke + pe
        
        energy_dissipated = total_work - total_e
        dissipation_percent = total_work > 0 ? (energy_dissipated / total_work * 100) : 0.0
        
        """
        Time: $(round(t, digits=3)) s
        
        Distributed Force Configuration:
        - Total magnitude: $(F_MAG) N
        - Global angle: $(FORCE_ANGLE_DEGREES)°
        - Distribution: F/20, F/10×9, F/20 on column 1
        - Components: $(round(fx_unit, digits=2)), $(round(fy_unit, digits=2))
        - Active: $(t <= F_ACTIVE_TIME ? "Yes" : "No")
        
        Spring Network:
        - Nearest: $(length(nearest_neighbor_connections))
        - Diagonal: $(length(diagonal_connections))
        - Total: $(length(nearest_neighbor_connections) + length(diagonal_connections))
        - Wall constraint: prevents rightmost column from passing through backplate
        - Damping: $(C_DAMPING) N·s/m (nearest neighbors, column-based scaling)
        - Material scaling: $(MATERIAL_MULTIPLIER)× multiplier (columns 2-11 scaled)
        
        Energy:
        Kinetic: $(round(ke, digits=6))
        Potential: $(round(pe, digits=6))
        Total: $(round(total_e, digits=6))
        
        Work Input: $(round(total_work, digits=6))
        Dissipated: $(round(energy_dissipated, digits=6))
        Dissipation: $(round(dissipation_percent, digits=2))%
        """
    end
    
    text!(stats_ax, 0.05, 0.95, text = stats_text, space = :relative, fontsize = 11)
    
    # Animation controls
    play_button = Button(fig[3, 1], label = "Play/Pause")
    reset_button = Button(fig[3, 2], label = "Reset")
    
    # Animation state
    is_playing = Observable(false)
    
    # Play/Pause functionality
    on(play_button.clicks) do _
        is_playing[] = !is_playing[]
    end
    
    # Reset functionality
    on(reset_button.clicks) do _
        current_frame[] = 1
        is_playing[] = false
    end
    
    # Animation loop
    @async begin
        while isopen(fig.scene)
            if is_playing[]
                if current_frame[] < length(sol.t)
                    current_frame[] += 1
                else
                    current_frame[] = 1  # Loop animation
                end
                sleep(OUTPUT_INTERVAL / ANIMATION_SPEED)
            else
                sleep(0.01)  # Small delay when paused
            end
        end
    end
    
    # Start animation automatically
    is_playing[] = true
    
    # Display figure
    display(fig)
    
    # Final energy analysis
    final_energy = energy_values[end]
    
    println("Final Summary:")
    println("=" ^ 50)
    println("Distributed Force Configuration:")
    println("  Total magnitude: $(F_MAG) N")
    println("  Global angle: $(FORCE_ANGLE_DEGREES)°")
    println("  Distribution: F/20 at rows 1 and 11, F/10 at rows 2-10")
    println("  Direction components: fx = $(round(fx_unit, digits=3)), fy = $(round(fy_unit, digits=3))")
    println()
    println("Spring Network Statistics:")
    println("  Nearest neighbor springs: $(length(nearest_neighbor_connections))")
    println("  Diagonal springs: $(length(diagonal_connections))")
    println("  Total springs: $(length(nearest_neighbor_connections) + length(diagonal_connections))")
    println("  Damping coefficient: $(C_DAMPING) N·s/m (nearest neighbors)")
    println("  Wall constraint: prevents rightmost column from passing through backplate")
    println()
    println("Energy Analysis:")
    println("  Total work input: $(total_work)")
    println("  Final total energy: $(final_energy)")
    println("  Energy dissipated: $(total_work - final_energy)")
    dissipation_percent = total_work > 0 ? ((total_work - final_energy) / total_work * 100) : 0.0
    println("  Energy dissipation: $(dissipation_percent)%")
    
    return fig, sol, total_work, final_energy
end


########################################################################
#  OPTIONAL: SAVE ANIMATION TO FILE WITH DISTRIBUTED LOAD AND BACKPLATE
########################################################################
function save_animation_to_file(filename = "lattice_anim_11x11_with_backplate.mp4")
    """
    Save the animation to a video file with distributed load display and backplate.
    """
    println("Creating animation file: $filename")
    
    # Ensure directory exists
    dir = dirname(filename)
    if !isempty(dir) && !isdir(dir)
        mkpath(dir)
    end
    
    # Solve the system
    u0 = zeros(2 * TOTAL_DOF)
    tspan = (0.0, T_END)
    prob = ODEProblem(lattice_2d_rhs_with_diagonals_and_backplate!, u0, tspan)
    sol = solve(prob, Vern9();
                reltol = REL_TOL, 
                abstol = ABS_TOL,
                saveat = OUTPUT_INTERVAL,
                dense = false)
    
    # Create equilibrium grid, spring connections, and backplate positions
    equilibrium_grid = create_equilibrium_grid()
    nearest_neighbor_connections, diagonal_connections = create_spring_connections_with_diagonals()
    backplate_positions = create_backplate_positions()
    
    # Extract position data
    all_positions = [reshape(view(sol.u[i], 1:TOTAL_DOF), 2, TOTAL_MASSES) for i in 1:length(sol.t)]
    
    # Create figure for recording
    fig = Figure(size = (1400, 900), fontsize = 14)
    ax = Axis(fig[1, 1], 
              title = "11×11 Mass-Spring Lattice with Distributed Load, Immovable Backplate, and Material Scaling",
              xlabel = "X Position",
              ylabel = "Y Position",
              aspect = DataAspect())
    
    # Set axis limits
    grid_min = minimum(equilibrium_grid) - 0.5
    grid_max_x = maximum(backplate_positions[1, :]) + 0.5
    grid_max_y = maximum(backplate_positions[2, :]) + 0.5
    limits!(ax, grid_min, grid_max_x, grid_min, grid_max_y)
    
    # Precompute force arrows (they're static, based on equilibrium positions)
    force_arrows = create_distributed_force_arrows(equilibrium_grid)
    fx_unit, fy_unit = calculate_force_components(1.0, FORCE_ANGLE_DEGREES)
    
    # Precompute colors (static, based on column indices)
    colors = Vector{Symbol}(undef, TOTAL_MASSES)
    column_colors = [:lightblue, :lightcyan, :cyan, :blue, :mediumblue, :darkblue, :navy, :midnightblue, :darkslateblue, :darkblue, :black]
    for k in 1:TOTAL_MASSES
        i, j = lattice_i(k), lattice_j(k)
        color_idx = min(j, length(column_colors))
        colors[k] = column_colors[color_idx]
    end
    
    # Precompute backplate line points (static)
    backplate_line_points = Point2f[backplate_positions[:, i] for i in 1:N]
    
    # Record animation
    record(fig, filename, 1:length(sol.t); framerate = 30) do frame
        empty!(ax)
        
        current_positions = all_positions[frame] .+ equilibrium_grid
        
        # Plot nearest neighbor springs
        for (k1, k2) in nearest_neighbor_connections
            spring_points = Point2f[current_positions[:, k1], current_positions[:, k2]]
            lines!(ax, spring_points, color = :blue, linewidth = SPRING_WIDTH)
        end
        
        # Plot diagonal springs
        for (k1, k2) in diagonal_connections
            spring_points = Point2f[current_positions[:, k1], current_positions[:, k2]]
            lines!(ax, spring_points, color = :red, linewidth = DIAGONAL_SPRING_WIDTH)
        end
        
        # Plot backplate as a vertical line (wall, no springs)
        lines!(ax, backplate_line_points, color = :gray, linewidth = 8)
        
        # Plot masses with column-based color coding
        mass_points = Point2f[current_positions[:, k] for k in 1:TOTAL_MASSES]
        scatter!(ax, mass_points, markersize = NODE_SIZE, color = colors, strokewidth = 2, strokecolor = :white)
        
        # Add distributed force arrows if active (use precomputed arrows)
        if sol.t[frame] <= F_ACTIVE_TIME
            for i in 1:N
                if length(force_arrows) >= 2*i
                    arrow_points = [force_arrows[2*i-1], force_arrows[2*i]]
                    lines!(ax, arrow_points, color = :green, linewidth = 3)
                    scatter!(ax, [force_arrows[2*i]], markersize = 8, color = :green, marker = :circle)
                end
            end
        end
        
        # Add title with current configuration
        ax.title = "Distributed Load: $(F_MAG)N total @ $(FORCE_ANGLE_DEGREES)° on column 1 - Time: $(round(sol.t[frame], digits=3))s"
    end
    
    println("Animation saved to: $filename")
end


########################################################################
#  RUN SIMULATION
########################################################################
# Run the interactive simulation
# Use semicolon to suppress verbose output in REPL
fig, sol, work_total, final_energy = run_2d_simulation_with_configurable_force_and_backplate();

# Save animation to file (automatically updates on each run)
# Use semicolon to suppress verbose output
save_animation_to_file("animations/lattice_anim_11x11_with_backplate.mp4");

