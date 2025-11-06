using DifferentialEquations
using LinearAlgebra
using Printf

########################################################################
#  PHYSICAL PARAMETERS
########################################################################
const MASS        = 1.0          # kg for every mass
const K_COUPLING  = 100.0        # N m⁻¹ (identical for all springs)
const F_ACTIVE_TIME = 0.05       # duration of the driving pulse (s)
const F_MAG       = 10.0         # N, applied diagonally to corner node

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
#  UTILITY: 2-D ↔ 1-D index conversion and state vector indexing
########################################################################
@inline lattice_idx(i,j) = (i-1)*N + j                    # (i,j) → mass index 1…25
@inline lattice_i(k) = 1 + div(k-1,N)                     # mass index → row
@inline lattice_j(k) = 1 + mod(k-1,N)                     # mass index → column

# State vector indexing for 2D motion
# State vector: [x1, y1, x2, y2, ..., x25, y25, vx1, vy1, vx2, vy2, ..., vx25, vy25]
@inline pos_x_idx(mass_idx) = 2 * mass_idx - 1            # x position index
@inline pos_y_idx(mass_idx) = 2 * mass_idx                # y position index
@inline vel_x_idx(mass_idx) = TOTAL_DOF + 2 * mass_idx - 1 # x velocity index
@inline vel_y_idx(mass_idx) = TOTAL_DOF + 2 * mass_idx     # y velocity index

########################################################################
#  2D SPRING FORCE CALCULATION
########################################################################
function spring_force_2d(pos1, pos2, k)
    """
    Calculate 2D spring force between two masses.
    Returns force on mass 1 due to mass 2.
    """
    displacement = pos2 - pos1  # vector from mass1 to mass2
    return k * displacement     # Hooke's law in 2D
end

########################################################################
#  ODE RIGHT-HAND SIDE FOR 2D MOTION
########################################################################
function lattice_2d_rhs!(du, u, p, t)
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
        
        # Horizontal springs (left-right neighbors)
        if j > 1  # left neighbor exists
            neighbor_idx = lattice_idx(i, j-1)
            pos_neighbor = pos[:, neighbor_idx]
            force = spring_force_2d(pos_k, pos_neighbor, K_COUPLING)
            dvel[:, k] += force
        end
        
        if j < N  # right neighbor exists
            neighbor_idx = lattice_idx(i, j+1)
            pos_neighbor = pos[:, neighbor_idx]
            force = spring_force_2d(pos_k, pos_neighbor, K_COUPLING)
            dvel[:, k] += force
        end
        
        # Vertical springs (up-down neighbors)
        if i > 1  # up neighbor exists
            neighbor_idx = lattice_idx(i-1, j)
            pos_neighbor = pos[:, neighbor_idx]
            force = spring_force_2d(pos_k, pos_neighbor, K_COUPLING)
            dvel[:, k] += force
        end
        
        if i < N  # down neighbor exists
            neighbor_idx = lattice_idx(i+1, j)
            pos_neighbor = pos[:, neighbor_idx]
            force = spring_force_2d(pos_k, pos_neighbor, K_COUPLING)
            dvel[:, k] += force
        end
    end

    # --- external driving force on corner node (1,1) ---
    if t <= F_ACTIVE_TIME
        corner_idx = lattice_idx(1, 1)
        # Apply diagonal force (equal x and y components)
        diagonal_force = F_MAG / sqrt(2)  # normalize to maintain total magnitude
        dvel[1, corner_idx] += diagonal_force  # x component
        dvel[2, corner_idx] += diagonal_force  # y component
    end

    # --- divide by mass to obtain accelerations ---
    dvel ./= MASS
    
    return nothing
end

########################################################################
#  ENERGY & WORK CALCULATIONS FOR 2D MOTION
########################################################################
function kinetic_energy_2d(vel_matrix)
    """
    Calculate total kinetic energy for 2D motion.
    vel_matrix is 2×25 matrix where each column is [vx, vy] for one mass.
    """
    return 0.5 * MASS * sum(vel_matrix.^2)
end

function potential_energy_2d(pos_matrix)
    """
    Calculate total potential energy for 2D spring system.
    pos_matrix is 2×25 matrix where each column is [x, y] for one mass.
    """
    pe = 0.0
    
    # Horizontal springs
    for i in 1:N, j in 1:N-1
        k1 = lattice_idx(i, j)
        k2 = lattice_idx(i, j+1)
        displacement = pos_matrix[:, k2] - pos_matrix[:, k1]
        pe += 0.5 * K_COUPLING * sum(displacement.^2)
    end
    
    # Vertical springs
    for i in 1:N-1, j in 1:N
        k1 = lattice_idx(i, j)
        k2 = lattice_idx(i+1, j)
        displacement = pos_matrix[:, k2] - pos_matrix[:, k1]
        pe += 0.5 * K_COUPLING * sum(displacement.^2)
    end
    
    return pe
end

function work_done_2d(sol)
    """
    Calculate total work done by diagonal external force.
    Work = F · dr for the driven corner mass.
    """
    w = 0.0
    corner_idx = lattice_idx(1, 1)
    diagonal_force = F_MAG / sqrt(2)
    
    for n in 2:length(sol.t)
        t₁, t₂ = sol.t[n-1], sol.t[n]
        if t₁ ≥ F_ACTIVE_TIME
            break
        end
        
        # Get positions at current and previous time
        pos_prev = reshape(view(sol.u[n-1], 1:TOTAL_DOF), 2, TOTAL_MASSES)
        pos_curr = reshape(view(sol.u[n], 1:TOTAL_DOF), 2, TOTAL_MASSES)
        
        # Calculate displacement of corner mass
        displacement = pos_curr[:, corner_idx] - pos_prev[:, corner_idx]
        
        # Work = force dot displacement
        force_vector = [diagonal_force, diagonal_force]
        w += dot(force_vector, displacement)
    end
    
    return w
end

########################################################################
#  MAIN SIMULATION DRIVER
########################################################################
function run_2d_simulation()
    println("5×5 Lattice of 2D Linear Mass–Spring Oscillators")
    println("=" ^ 65)
    println("Physical Parameters:")
    println("  Total masses: $(TOTAL_MASSES)")
    println("  DOF per mass: $(DOF_PER_MASS) (x, y)")
    println("  Total DOF: $(TOTAL_DOF)")
    println("  Total state vector size: $(2 * TOTAL_DOF)")
    println("  Spring constant: $(K_COUPLING) N/m")
    println("  Mass: $(MASS) kg")
    println("")
    println("External Force:")
    println("  Magnitude: $(F_MAG) N")
    println("  Direction: Diagonal (45°)")
    println("  Applied to: Corner mass (1,1)")
    println("  Duration: $(F_ACTIVE_TIME) s")
    println("")
    println("Simulation Parameters:")
    println("  Total time: $(T_END) s")
    println("  Output interval: $(OUTPUT_INTERVAL) s")
    println("  Relative tolerance: $(REL_TOL)")
    println("  Absolute tolerance: $(ABS_TOL)")
    println("")

    # Initial state: all masses at rest at origin
    u0 = zeros(2 * TOTAL_DOF)
    tspan = (0.0, T_END)

    # Create and solve ODE problem
    prob = ODEProblem(lattice_2d_rhs!, u0, tspan)
    sol = solve(prob, Vern9();
                reltol = REL_TOL, 
                abstol = ABS_TOL,
                saveat = OUTPUT_INTERVAL)

    if sol.retcode != :Success
        error("Solver failure: ", sol.retcode)
    end

    println("Solver completed successfully!")
    println("Final time: $(sol.t[end]) s")
    println("Time steps: $(length(sol.t))")
    println("")

    # Calculate total work done by external force
    total_work = work_done_2d(sol)

    # Display results
    println("    t (s)   |       KE           PE          E_tot      |  Work    |  %Error")
    println("-" ^ 80)

    max_error = 0.0
    for (i, t) in enumerate(sol.t)
        # Extract positions and velocities
        pos = reshape(view(sol.u[i], 1:TOTAL_DOF), 2, TOTAL_MASSES)
        vel = reshape(view(sol.u[i], TOTAL_DOF+1:2*TOTAL_DOF), 2, TOTAL_MASSES)
        
        # Calculate energies
        ke = kinetic_energy_2d(vel)
        pe = potential_energy_2d(pos)
        total_energy = ke + pe
        
        # Energy conservation error
        if abs(total_work) > 1e-16
            percent_error = abs(total_work - total_energy) / abs(total_work) * 100
            max_error = max(max_error, percent_error)
        else
            percent_error = 0.0
        end
        
        if i % 20 == 1  # Print every 0.2 s
            @printf("%9.2f  | %12.6e %12.6e %12.6e | %8.6f | %7.3f\n",
                    t, ke, pe, total_energy, total_work, percent_error)
        end
    end

    println("\nFinal Summary:")
    println("=" ^ 50)
    
    # Final state analysis
    pos_final = reshape(view(sol.u[end], 1:TOTAL_DOF), 2, TOTAL_MASSES)
    vel_final = reshape(view(sol.u[end], TOTAL_DOF+1:2*TOTAL_DOF), 2, TOTAL_MASSES)
    
    ke_final = kinetic_energy_2d(vel_final)
    pe_final = potential_energy_2d(pos_final)
    energy_final = ke_final + pe_final
    
    final_error = abs(total_work - energy_final) / abs(total_work) * 100
    
    println("Total work input:    $(total_work)")
    println("Final kinetic energy: $(ke_final)")
    println("Final potential energy: $(pe_final)")
    println("Final total energy:   $(energy_final)")
    println("Energy error:         $(abs(total_work - energy_final))")
    println("Final % error:        $(final_error)%")
    println("Maximum % error:      $(max_error)%")
    
    # Display final positions of corner masses for verification
    println("\nFinal positions of corner masses:")
    for (desc, i, j) in [("Top-left", 1, 1), ("Top-right", 1, N), 
                         ("Bottom-left", N, 1), ("Bottom-right", N, N)]
        k = lattice_idx(i, j)
        pos = pos_final[:, k]
        @printf("  %s (%d,%d): x=%.6f, y=%.6f\n", desc, i, j, pos[1], pos[2])
    end
    
    return sol, total_work, energy_final
end

########################################################################
#  RUN SIMULATION
########################################################################
sol, work_total, final_energy = run_2d_simulation()
