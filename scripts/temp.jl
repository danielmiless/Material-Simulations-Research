using DifferentialEquations
using LinearAlgebra
using Printf
using GLMakie
using Observables


########################################################################
#  PHYSICAL PARAMETERS
########################################################################
const MASS        = 1.0          # kg for every mass
const K_COUPLING  = 100.0        # N for nearest neighbor (horizontal/vertical) springs
const K_DIAGONAL  = 50.0         # N for diagonal springs (typically weaker)
const ALPHA_COUPLING = 10.0      # exponential decay rate for nearest neighbor springs (m⁻¹)
const ALPHA_DIAGONAL = 10.0      # exponential decay rate for diagonal springs (m⁻¹)
const C_DAMPING   = 5.0          # N·s/m damping coefficient for nearest neighbor springs
const F_ACTIVE_TIME = 0.05       # duration of the driving pulse (s)
const F_MAG       = 250.0         # N, magnitude of external force


########################################################################
#  FORCE CONFIGURATION PARAMETERS - EASILY CHANGEABLE
########################################################################
const FORCE_ANGLE_DEGREES = 0.0    # Angle in degrees (0-360): 0=right, 90=up, 180=left, 270=down
const FORCE_TARGET_ROW = 3           # Row of target mass (1 to N)
const FORCE_TARGET_COL = 1           # Column of target mass (1 to N)


########################################################################
#  LATTICE SIZE  
########################################################################
const N = 5                      # ⇒ 5×5 = 25 masses
const TOTAL_MASSES = N * N       # 25 masses total
const DOF_PER_MASS = 2          # x and y components
const TOTAL_DOF = TOTAL_MASSES * DOF_PER_MASS  # 50 position DOF


########################################################################
#  NUMERICAL PARAMETERS
########################################################################
const T_END          = 5.0       # s
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
@inline lattice_idx(i,j) = (i-1)*N + j                    # (i,j) → mass index 1…25
@inline lattice_i(k) = 1 + div(k-1,N)                     # mass index → row
@inline lattice_j(k) = 1 + mod(k-1,N)                     # mass index → column


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
    Calculate 2D damping force between two masses.
    Returns damping force on mass 1 due to relative velocity with mass 2.
    
    F_damping = -c * (v1 - v2)
    """
    relative_velocity = vel1 - vel2
    return -c * relative_velocity
end


########################################################################
#  ODE RIGHT-HAND SIDE FOR 2D MOTION WITH DIAGONAL SPRINGS
########################################################################
function lattice_2d_rhs_with_diagonals!(du, u, p, t)
    # Extract positions and velocities from state vector
    pos = reshape(view(u, 1:TOTAL_DOF), 2, TOTAL_MASSES)           # 2×25 matrix
    vel = reshape(view(u, TOTAL_DOF+1:2*TOTAL_DOF), 2, TOTAL_MASSES) # 2×25 matrix
    
    dpos = reshape(view(du, 1:TOTAL_DOF), 2, TOTAL_MASSES)
    dvel = reshape(view(du, TOTAL_DOF+1:2*TOTAL_DOF), 2, TOTAL_MASSES)


    # --- kinematics: dx/dt = v ---
    dpos .= vel


    # --- initialize forces to zero ---
    dvel .= 0.0


    # --- calculate spring forces between neighboring masses ---
    for k in 1:TOTAL_MASSES
        i, j = lattice_i(k), lattice_j(k)
        pos_k = pos[:, k]  # position of mass k as [x, y]
        vel_k = vel[:, k]  # velocity of mass k as [vx, vy]
        
        # NEAREST NEIGHBOR SPRINGS (horizontal and vertical) WITH DAMPING
        # Horizontal springs (left-right neighbors)
        if j > 1  # left neighbor exists
            neighbor_idx = lattice_idx(i, j-1)
            pos_neighbor = pos[:, neighbor_idx]
            vel_neighbor = vel[:, neighbor_idx]
            force = spring_force_2d(pos_k, pos_neighbor, K_COUPLING, ALPHA_COUPLING)
            damping = damping_force_2d(vel_k, vel_neighbor, C_DAMPING)
            dvel[:, k] += force + damping
        end
        
        if j < N  # right neighbor exists
            neighbor_idx = lattice_idx(i, j+1)
            pos_neighbor = pos[:, neighbor_idx]
            vel_neighbor = vel[:, neighbor_idx]
            force = spring_force_2d(pos_k, pos_neighbor, K_COUPLING, ALPHA_COUPLING)
            damping = damping_force_2d(vel_k, vel_neighbor, C_DAMPING)
            dvel[:, k] += force + damping
        end
        
        # Vertical springs (up-down neighbors)
        if i > 1  # up neighbor exists
            neighbor_idx = lattice_idx(i-1, j)
            pos_neighbor = pos[:, neighbor_idx]
            vel_neighbor = vel[:, neighbor_idx]
            force = spring_force_2d(pos_k, pos_neighbor, K_COUPLING, ALPHA_COUPLING)
            damping = damping_force_2d(vel_k, vel_neighbor, C_DAMPING)
            dvel[:, k] += force + damping
        end
        
        if i < N  # down neighbor exists
            neighbor_idx = lattice_idx(i+1, j)
            pos_neighbor = pos[:, neighbor_idx]
            vel_neighbor = vel[:, neighbor_idx]
            force = spring_force_2d(pos_k, pos_neighbor, K_COUPLING, ALPHA_COUPLING)
            damping = damping_force_2d(vel_k, vel_neighbor, C_DAMPING)
            dvel[:, k] += force + damping
        end
        
        # DIAGONAL SPRINGS (next-nearest neighbors in X pattern) WITHOUT DAMPING
        # Upper-left diagonal
        if i > 1 && j > 1  # upper-left diagonal neighbor exists
            neighbor_idx = lattice_idx(i-1, j-1)
            pos_neighbor = pos[:, neighbor_idx]
            force = spring_force_2d(pos_k, pos_neighbor, K_DIAGONAL, ALPHA_DIAGONAL)
            dvel[:, k] += force
        end
        
        # Upper-right diagonal
        if i > 1 && j < N  # upper-right diagonal neighbor exists
            neighbor_idx = lattice_idx(i-1, j+1)
            pos_neighbor = pos[:, neighbor_idx]
            force = spring_force_2d(pos_k, pos_neighbor, K_DIAGONAL, ALPHA_DIAGONAL)
            dvel[:, k] += force
        end
        
        # Lower-left diagonal
        if i < N && j > 1  # lower-left diagonal neighbor exists
            neighbor_idx = lattice_idx(i+1, j-1)
            pos_neighbor = pos[:, neighbor_idx]
            force = spring_force_2d(pos_k, pos_neighbor, K_DIAGONAL, ALPHA_DIAGONAL)
            dvel[:, k] += force
        end
        
        # Lower-right diagonal
        if i < N && j < N  # lower-right diagonal neighbor exists
            neighbor_idx = lattice_idx(i+1, j+1)
            pos_neighbor = pos[:, neighbor_idx]
            force = spring_force_2d(pos_k, pos_neighbor, K_DIAGONAL, ALPHA_DIAGONAL)
            dvel[:, k] += force
        end
    end


    # --- CONFIGURABLE EXTERNAL DRIVING FORCE ---
    if t <= F_ACTIVE_TIME
        # Calculate force components from angle and magnitude
        fx, fy = calculate_force_components(F_MAG, FORCE_ANGLE_DEGREES)
        
        # Apply to specified target mass
        target_idx = lattice_idx(FORCE_TARGET_ROW, FORCE_TARGET_COL)
        dvel[1, target_idx] += fx  # x component
        dvel[2, target_idx] += fy  # y component
    end


    # --- divide by mass to obtain accelerations ---
    dvel ./= MASS
    
    return nothing
end


########################################################################
#  ENERGY & WORK CALCULATIONS FOR 2D MOTION WITH DIAGONAL SPRINGS
########################################################################
function kinetic_energy_2d(vel_matrix)
    """
    Calculate total kinetic energy for 2D motion.
    vel_matrix is 2×25 matrix where each column is [vx, vy] for one mass.
    """
    return 0.5 * MASS * sum(vel_matrix.^2)
end


function potential_energy_2d_with_diagonals(pos_matrix)
    """
    Calculate total potential energy for 2D exponential spring system with diagonal springs.
    pos_matrix is 2×25 matrix where each column is [x, y] for one mass.
    
    For exponential spring: U = (k/alpha) * (exp(alpha*|r|) - alpha*|r| - 1)
    """
    pe = 0.0
    
    # Horizontal springs (nearest neighbors)
    for i in 1:N, j in 1:N-1
        k1 = lattice_idx(i, j)
        k2 = lattice_idx(i, j+1)
        displacement = pos_matrix[:, k2] - pos_matrix[:, k1]
        distance = norm(displacement)
        pe += (K_COUPLING / ALPHA_COUPLING) * (exp(ALPHA_COUPLING * distance) - ALPHA_COUPLING * distance - 1.0)
    end
    
    # Vertical springs (nearest neighbors)
    for i in 1:N-1, j in 1:N
        k1 = lattice_idx(i, j)
        k2 = lattice_idx(i+1, j)
        displacement = pos_matrix[:, k2] - pos_matrix[:, k1]
        distance = norm(displacement)
        pe += (K_COUPLING / ALPHA_COUPLING) * (exp(ALPHA_COUPLING * distance) - ALPHA_COUPLING * distance - 1.0)
    end
    
    # Diagonal springs (next-nearest neighbors)
    # Upper-left to lower-right diagonals
    for i in 1:N-1, j in 1:N-1
        k1 = lattice_idx(i, j)
        k2 = lattice_idx(i+1, j+1)
        displacement = pos_matrix[:, k2] - pos_matrix[:, k1]
        distance = norm(displacement)
        pe += (K_DIAGONAL / ALPHA_DIAGONAL) * (exp(ALPHA_DIAGONAL * distance) - ALPHA_DIAGONAL * distance - 1.0)
    end
    
    # Upper-right to lower-left diagonals
    for i in 1:N-1, j in 2:N
        k1 = lattice_idx(i, j)
        k2 = lattice_idx(i+1, j-1)
        displacement = pos_matrix[:, k2] - pos_matrix[:, k1]
        distance = norm(displacement)
        pe += (K_DIAGONAL / ALPHA_DIAGONAL) * (exp(ALPHA_DIAGONAL * distance) - ALPHA_DIAGONAL * distance - 1.0)
    end
    
    return pe
end


function work_done_2d_configurable(sol)
    """
    Calculate total work done by configurable external force.
    Work = F · dr for the driven target mass.
    """
    w = 0.0
    target_idx = lattice_idx(FORCE_TARGET_ROW, FORCE_TARGET_COL)
    fx, fy = calculate_force_components(F_MAG, FORCE_ANGLE_DEGREES)
    
    for n in 2:length(sol.t)
        t₁, t₂ = sol.t[n-1], sol.t[n]
        if t₁ ≥ F_ACTIVE_TIME
            break
        end
        
        # Get positions at current and previous time
        pos_prev = reshape(view(sol.u[n-1], 1:TOTAL_DOF), 2, TOTAL_MASSES)
        pos_curr = reshape(view(sol.u[n], 1:TOTAL_DOF), 2, TOTAL_MASSES)
        
        # Calculate displacement of target mass
        displacement = pos_curr[:, target_idx] - pos_prev[:, target_idx]
        
        # Work = force dot displacement
        force_vector = [fx, fy]
        w += dot(force_vector, displacement)
    end
    
    return w
end


########################################################################
#  EQUILIBRIUM GRID POSITIONS
########################################################################
function create_equilibrium_grid()
    """
    Create equilibrium positions for the 5×5 lattice.
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
#  SPRING CONNECTIVITY INCLUDING DIAGONALS
########################################################################
function create_spring_connections_with_diagonals()
    """
    Create list of spring connections for visualization, including diagonal springs.
    Returns separate arrays for nearest-neighbor and diagonal connections.
    """
    nearest_neighbor_connections = Tuple{Int, Int}[]
    diagonal_connections = Tuple{Int, Int}[]
    
    # Horizontal springs (nearest neighbors)
    for i in 1:N, j in 1:N-1
        k1 = lattice_idx(i, j)
        k2 = lattice_idx(i, j+1)
        push!(nearest_neighbor_connections, (k1, k2))
    end
    
    # Vertical springs (nearest neighbors)
    for i in 1:N-1, j in 1:N
        k1 = lattice_idx(i, j)
        k2 = lattice_idx(i+1, j)
        push!(nearest_neighbor_connections, (k1, k2))
    end
    
    # Diagonal springs (next-nearest neighbors)
    # Upper-left to lower-right diagonals
    for i in 1:N-1, j in 1:N-1
        k1 = lattice_idx(i, j)
        k2 = lattice_idx(i+1, j+1)
        push!(diagonal_connections, (k1, k2))
    end
    
    # Upper-right to lower-left diagonals
    for i in 1:N-1, j in 2:N
        k1 = lattice_idx(i, j)
        k2 = lattice_idx(i+1, j-1)
        push!(diagonal_connections, (k1, k2))
    end
    
    return nearest_neighbor_connections, diagonal_connections
end


########################################################################
#  FORCE VISUALIZATION HELPER
########################################################################
function create_force_arrow_points(equilibrium_grid, target_idx, force_magnitude, angle_degrees)
    """
    Create points for visualizing the external force as an arrow.
    """
    # Get equilibrium position of target mass
    target_pos = equilibrium_grid[:, target_idx]
    
    # Calculate force components
    fx, fy = calculate_force_components(force_magnitude, angle_degrees)
    
    # Scale for visualization (make arrow visible but not too long)
    scale_factor = 0.3
    force_end = target_pos + scale_factor * [fx, fy] / force_magnitude
    
    return [Point2f(target_pos), Point2f(force_end)]
end


########################################################################
#  MAIN SIMULATION AND ANIMATION FUNCTION
########################################################################
function run_2d_simulation_with_configurable_force()
    """
    Main simulation function with configurable force angle and application point.
    """
    println("5×5 Lattice Mass-Spring System with Configurable Force")
    println("=" ^ 65)
    println("Physical Parameters:")
    println("  Total masses: $(TOTAL_MASSES)")
    println("  DOF per mass: $(DOF_PER_MASS) (x, y)")
    println("  Total DOF: $(TOTAL_DOF)")
    println("  Nearest neighbor spring constant: $(K_COUPLING) N")
    println("  Nearest neighbor alpha: $(ALPHA_COUPLING) m⁻¹")
    println("  Nearest neighbor damping: $(C_DAMPING) N·s/m")
    println("  Diagonal spring constant: $(K_DIAGONAL) N")
    println("  Diagonal alpha: $(ALPHA_DIAGONAL) m⁻¹")
    println("  Mass: $(MASS) kg")
    println("")
    println("Spring Network:")
    println("  Nearest neighbors: horizontal & vertical connections (with damping)")
    println("  Diagonal springs: X-pattern connections in each square (no damping)")
    println("  Total connectivity: 8 neighbors per interior mass")
    println("  Spring type: Exponential")
    println("")
    println("Configurable External Force:")
    println("  Magnitude: $(F_MAG) N")
    println("  Angle: $(FORCE_ANGLE_DEGREES)° (0°=right, 90°=up, 180°=left, 270°=down)")
    println("  Applied to: Mass at position ($(FORCE_TARGET_ROW), $(FORCE_TARGET_COL))")
    println("  Duration: $(F_ACTIVE_TIME) s")
    
    # Calculate and display force components
    fx, fy = calculate_force_components(F_MAG, FORCE_ANGLE_DEGREES)
    println("  Force components: Fx = $(round(fx, digits=3)) N, Fy = $(round(fy, digits=3)) N")
    println("")


    # Initial state: all masses at rest at equilibrium positions
    u0 = zeros(2 * TOTAL_DOF)
    tspan = (0.0, T_END)


    # Create and solve ODE problem
    prob = ODEProblem(lattice_2d_rhs_with_diagonals!, u0, tspan)
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


    # Calculate total work done by external force
    total_work = work_done_2d_configurable(sol)
    
    # Create equilibrium grid and spring connections
    equilibrium_grid = create_equilibrium_grid()
    nearest_neighbor_connections, diagonal_connections = create_spring_connections_with_diagonals()
    
    # Set up GLMakie figure and layout
    GLMakie.activate!()
    fig = Figure(size = (1400, 900), fontsize = 14)
    
    # Main animation axis
    ax = Axis(fig[1, 1:2], 
              title = "2D Mass-Spring Lattice with Configurable Force",
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
    
    # Extract position data for all time steps
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
    
    # Plot masses as scatter points
    mass_points = @lift(Point2f[($current_positions)[:, k] for k in 1:TOTAL_MASSES])
    
    # Color target mass differently to show where force is applied
    colors = fill(:black, TOTAL_MASSES)
    target_idx = lattice_idx(FORCE_TARGET_ROW, FORCE_TARGET_COL)
    colors[target_idx] = :orange
    
    scatter!(ax, mass_points, markersize = NODE_SIZE, color = colors, strokewidth = 2, strokecolor = :white)
    
    # Add force arrow visualization (visible only during active time)
    force_active = @lift(sol.t[$current_frame] <= F_ACTIVE_TIME)
    force_arrow_points = @lift(begin
        if $force_active
            create_force_arrow_points(equilibrium_grid, target_idx, F_MAG, FORCE_ANGLE_DEGREES)
        else
            Point2f[]
        end
    end)
    
    # Plot force arrow
    lines!(ax, force_arrow_points, color = :green, linewidth = 4, label = "External Force")
    
    # Add arrow head
    force_arrow_end = @lift(begin
        if $force_active && length($force_arrow_points) >= 2
            [$force_arrow_points[2]]
        else
            Point2f[]
        end
    end)
    
    scatter!(ax, force_arrow_end, markersize = 12, color = :green, marker = :circle)
    
    # Add legend
    axislegend(ax, position = :rt)
    
    # Set up axis limits with some padding
    grid_min = minimum(equilibrium_grid) - 0.5
    grid_max = maximum(equilibrium_grid) + 0.5
    limits!(ax, grid_min, grid_max, grid_min, grid_max)
    
    # Energy conservation plot
    time_values = sol.t
    energy_values = [let
        pos = reshape(view(sol.u[i], 1:TOTAL_DOF), 2, TOTAL_MASSES)
        vel = reshape(view(sol.u[i], TOTAL_DOF+1:2*TOTAL_DOF), 2, TOTAL_MASSES)
        kinetic_energy_2d(vel) + potential_energy_2d_with_diagonals(pos)
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
        pe = potential_energy_2d_with_diagonals(pos)
        total_e = ke + pe
        
        error_percent = abs(total_work - total_e) / abs(total_work) * 100
        
        target_pos = pos[:, target_idx]
        
        """
        Time: $(round(t, digits=3)) s
        
        Force Configuration:
        - Magnitude: $(F_MAG) N
        - Angle: $(FORCE_ANGLE_DEGREES)°
        - Target: ($(FORCE_TARGET_ROW), $(FORCE_TARGET_COL))
        - Components: $(round(fx, digits=2)), $(round(fy, digits=2)) N
        - Active: $(t <= F_ACTIVE_TIME ? "Yes" : "No")
        
        Spring Network:
        - Nearest: $(length(nearest_neighbor_connections))
        - Diagonal: $(length(diagonal_connections))
        - Total: $(length(nearest_neighbor_connections) + length(diagonal_connections))
        - Damping: $(C_DAMPING) N·s/m (nearest only)
        
        Energy:
        Kinetic: $(round(ke, digits=6))
        Potential: $(round(pe, digits=6))
        Total: $(round(total_e, digits=6))
        
        Work Input: $(round(total_work, digits=6))
        Error: $(round(error_percent, digits=4))%
        
        Target Mass:
        X: $(round(target_pos[1], digits=4))
        Y: $(round(target_pos[2], digits=4))
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
    final_error = abs(total_work - final_energy) / abs(total_work) * 100
    
    println("Final Summary:")
    println("=" ^ 50)
    println("Force Configuration:")
    println("  Magnitude: $(F_MAG) N")
    println("  Angle: $(FORCE_ANGLE_DEGREES)°")
    println("  Target mass: ($(FORCE_TARGET_ROW), $(FORCE_TARGET_COL))")
    println("  Components: Fx = $(round(fx, digits=3)) N, Fy = $(round(fy, digits=3)) N")
    println()
    println("Spring Network Statistics:")
    println("  Nearest neighbor springs: $(length(nearest_neighbor_connections))")
    println("  Diagonal springs: $(length(diagonal_connections))")
    println("  Total springs: $(length(nearest_neighbor_connections) + length(diagonal_connections))")
    println("  Damping coefficient: $(C_DAMPING) N·s/m (nearest neighbors only)")
    println()
    println("Energy Conservation:")
    println("  Total work input: $(total_work)")
    println("  Final total energy: $(final_energy)")
    println("  Energy error: $(abs(total_work - final_energy))")
    println("  Final % error: $(final_error)%")
    
    return fig, sol, total_work, final_energy
end


########################################################################
#  OPTIONAL: SAVE ANIMATION TO FILE WITH FORCE CONFIGURATION
########################################################################
function save_animation_to_file(filename = "lattice_anim.mp4")
    """
    Save the animation to a video file with force configuration display.
    """
    println("Creating animation file: $filename")
    
    # Solve the system
    u0 = zeros(2 * TOTAL_DOF)
    tspan = (0.0, T_END)
    prob = ODEProblem(lattice_2d_rhs_with_diagonals!, u0, tspan)
    sol = solve(prob, Vern9();
                reltol = REL_TOL, 
                abstol = ABS_TOL,
                saveat = OUTPUT_INTERVAL,
                dense = false)
    
    # Create equilibrium grid and spring connections
    equilibrium_grid = create_equilibrium_grid()
    nearest_neighbor_connections, diagonal_connections = create_spring_connections_with_diagonals()
    
    # Extract position data
    all_positions = [reshape(view(sol.u[i], 1:TOTAL_DOF), 2, TOTAL_MASSES) for i in 1:length(sol.t)]
    
    # Create figure for recording
    fig = Figure(size = (1400, 900), fontsize = 14)
    ax = Axis(fig[1, 1], 
              title = "2D Mass-Spring Lattice with Configurable Force",
              xlabel = "X Position",
              ylabel = "Y Position",
              aspect = DataAspect())
    
    # Set axis limits
    grid_min = minimum(equilibrium_grid) - 0.5
    grid_max = maximum(equilibrium_grid) + 0.5
    limits!(ax, grid_min, grid_max, grid_min, grid_max)
    
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
        
        # Plot masses
        mass_points = Point2f[current_positions[:, k] for k in 1:TOTAL_MASSES]
        colors = fill(:black, TOTAL_MASSES)
        target_idx = lattice_idx(FORCE_TARGET_ROW, FORCE_TARGET_COL)
        colors[target_idx] = :orange
        scatter!(ax, mass_points, markersize = NODE_SIZE, color = colors, strokewidth = 2, strokecolor = :white)
        
        # Add force arrow if active
        if sol.t[frame] <= F_ACTIVE_TIME
            force_arrow_points = create_force_arrow_points(equilibrium_grid, target_idx, F_MAG, FORCE_ANGLE_DEGREES)
            lines!(ax, force_arrow_points, color = :green, linewidth = 4)
            if length(force_arrow_points) >= 2
                scatter!(ax, [force_arrow_points[2]], markersize = 12, color = :green, marker = :circle)
            end
        end
        
        # Add title with current configuration
        fx, fy = calculate_force_components(F_MAG, FORCE_ANGLE_DEGREES)
        ax.title = "Force: $(F_MAG)N @ $(FORCE_ANGLE_DEGREES)° on ($(FORCE_TARGET_ROW),$(FORCE_TARGET_COL)) - Time: $(round(sol.t[frame], digits=3))s"
    end
    
    println("Animation saved to: $filename")
end


########################################################################
#  RUN SIMULATION
########################################################################
# Run the interactive simulation
fig, sol, work_total, final_energy = run_2d_simulation_with_configurable_force()


# Uncomment the following line to save animation to file
save_animation_to_file("animations/lattice_anim.mp4")
