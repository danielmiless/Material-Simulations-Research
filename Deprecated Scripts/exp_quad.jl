using DifferentialEquations
using DiffEqCallbacks
using SparseArrays
using LinearAlgebra
using Printf

# Physical parameters for exponential springs
const MASS = 1.0
const K_VAL = 100.0         # Reference spring constant for ground spring
const K_COUPLING = 100.0    # Reference spring constant for coupling springs
const ALPHA_GROUND = 2.0    # Exponential parameter for ground spring
const ALPHA_COUPLING = 2.0  # Exponential parameter for coupling springs
const F_ACTIVE_TIME = 0.05

# Simulation parameters with maximum precision
const T_END = 5.0
const OUTPUT_INTERVAL = 0.01

# Ultra-tight tolerances for energy conservation
const REL_TOL = 1e-14
const ABS_TOL = 1e-16

function exponential_spring_force(displacement, k0, alpha)
    """
    Calculate force for exponential spring: F = k0 * (exp(alpha * x) - 1)
    This ensures F = 0 at equilibrium (x = 0) and exponential growth away from equilibrium.
    """
    return k0 * (exp(alpha * displacement) - 1.0)
end

function exponential_spring_potential(displacement, k0, alpha)
    """
    Calculate potential energy for exponential spring.
    U = k0 * (exp(alpha * x)/alpha - x - 1/alpha)
    Setting U = 0 at x = 0
    """
    if abs(alpha) < 1e-12
        # Linear limit as alpha → 0: U ≈ k0 * x²/2
        return 0.5 * k0 * displacement^2
    else
        return k0 * ((exp(alpha * displacement) / alpha) - displacement - (1.0 / alpha))
    end
end

function exponential_spring_work(x1, x2, k0, alpha)
    """
    Calculate exact work done by exponential spring moving from x1 to x2.
    W = ∫F dx = ∫k0 * (exp(alpha * x) - 1) dx = k0 * (exp(alpha * x)/alpha - x) + C
    """
    if abs(alpha) < 1e-12
        # Linear limit: W = k0 * (x2² - x1²)/2
        return 0.5 * k0 * (x2^2 - x1^2)
    else
        work_x2 = k0 * (exp(alpha * x2) / alpha - x2)
        work_x1 = k0 * (exp(alpha * x1) / alpha - x1)
        return work_x2 - work_x1
    end
end

function mass_spring_system_exponential!(du, u, p, t)
    """
    Defines the 4-mass exponential spring system dynamics with exact force calculations.
    """
    
    # Extract positions and velocities
    x1, v1, x2, v2, x3, v3, x4, v4 = u
    
    # External force (applied only to first mass)
    external_force = (t <= F_ACTIVE_TIME) ? 10.0 : 0.0
    
    # Calculate exponential spring forces
    # Ground spring force (on mass 1)
    f_ground = -exponential_spring_force(x1, K_VAL, ALPHA_GROUND)
    
    # Coupling spring forces
    f_12 = exponential_spring_force(x2 - x1, K_COUPLING, ALPHA_COUPLING)  # Spring 1-2
    f_23 = exponential_spring_force(x3 - x2, K_COUPLING, ALPHA_COUPLING)  # Spring 2-3  
    f_34 = exponential_spring_force(x4 - x3, K_COUPLING, ALPHA_COUPLING)  # Spring 3-4
    
    # Equations of motion
    # Mass 1: ground spring + coupling to mass 2 + external force
    du[1] = v1  # dx1/dt = v1
    du[2] = (f_ground + f_12 + external_force) / MASS
    
    # Mass 2: coupling forces from masses 1 and 3
    du[3] = v2  # dx2/dt = v2
    du[4] = (-f_12 + f_23) / MASS
    
    # Mass 3: coupling forces from masses 2 and 4
    du[5] = v3  # dx3/dt = v3
    du[6] = (-f_23 + f_34) / MASS
    
    # Mass 4: coupling force from mass 3 only
    du[7] = v4  # dx4/dt = v4
    du[8] = (-f_34) / MASS
    
    return nothing
end

function calculate_energies_exponential(u, t)
    """
    Calculate kinetic energy, potential energy, and total energy for exponential springs.
    """
    x1, v1, x2, v2, x3, v3, x4, v4 = u
    
    # Kinetic energy
    KE = 0.5 * MASS * (v1^2 + v2^2 + v3^2 + v4^2)
    
    # Potential energy from exponential springs
    PE_ground = exponential_spring_potential(x1, K_VAL, ALPHA_GROUND)
    PE_12 = exponential_spring_potential(x2 - x1, K_COUPLING, ALPHA_COUPLING)
    PE_23 = exponential_spring_potential(x3 - x2, K_COUPLING, ALPHA_COUPLING)  
    PE_34 = exponential_spring_potential(x4 - x3, K_COUPLING, ALPHA_COUPLING)
    
    PE = PE_ground + PE_12 + PE_23 + PE_34
    
    return KE, PE, KE + PE
end

function calculate_total_work_exact(solution_times, solution_states)
    """
    Calculate total work done by external force using exact integration.
    This is the key fix for energy conservation!
    """
    work_total = 0.0
    
    for i in 2:length(solution_times)
        t_prev, t_curr = solution_times[i-1], solution_times[i]
        x1_prev, x1_curr = solution_states[i-1][1], solution_states[i][1]
        
        # Only calculate work when force is active
        if t_prev < F_ACTIVE_TIME
            if t_curr <= F_ACTIVE_TIME
                # Both times in active region
                work_increment = 10.0 * (x1_curr - x1_prev)
                
            else
                # Transition across cutoff - use linear interpolation
                dt_total = t_curr - t_prev
                dt_active = F_ACTIVE_TIME - t_prev
                
                # Linear interpolation to find position at cutoff
                x1_cutoff = x1_prev + (x1_curr - x1_prev) * (dt_active / dt_total)
                
                # Work only during active period
                work_increment = 10.0 * (x1_cutoff - x1_prev)
            end
            
            work_total += work_increment
        end
    end
    
    return work_total
end

function run_exponential_spring_simulation()
    """
    Main simulation function for exponential spring system with exact energy conservation.
    """
    println("Ultra-High Precision Exponential Spring System Simulation")
    println("=" ^ 60)
    println("Physical Parameters:")
    println("  Mass: $MASS kg")
    println("  Ground spring k₀: $K_VAL N/m, α: $ALPHA_GROUND")
    println("  Coupling spring k₀: $K_COUPLING N/m, α: $ALPHA_COUPLING")
    println("  External force: 10 N (applied for $(F_ACTIVE_TIME) s)")
    println()
    println("Exponential Force Law: F = k₀ * (exp(α * x) - 1)")
    println("Key Fix: Using exact work integration instead of linear approximation")
    println()
    println("Simulation Parameters:")
    println("  Total time: $T_END s")
    println("  Output interval: $OUTPUT_INTERVAL s")
    println("  Relative tolerance: $REL_TOL")
    println("  Absolute tolerance: $ABS_TOL")
    println()
    
    # Initial conditions (all masses at rest)
    u0 = zeros(8)
    tspan = (0.0, T_END)
    
    # Create ODE problem
    prob = ODEProblem(mass_spring_system_exponential!, u0, tspan)
    
    # Use highest-order solver with ultra-tight tolerances
    sol = solve(prob, Vern9();
                reltol = REL_TOL,
                abstol = ABS_TOL,
                saveat = OUTPUT_INTERVAL,
                dense = false,
                maxiters = 10000000,
                force_dtmin = true)
    
    if sol.retcode != :Success
        error("Solver failed: $(sol.retcode)")
    end
    
    println("Solver completed successfully!")
    println("Final time: $(sol.t[end])")
    println("Time steps: $(length(sol.t))")
    println()
    
    # Calculate exact total work using proper integration
    work_total = calculate_total_work_exact(sol.t, sol.u)
    
    println("Time    X1      V1        X2      V2        Energy    Work      %Error")
    println("-" ^ 75)
    
    # Process and display results
    max_error = 0.0
    for i in 1:length(sol.t)
        t = sol.t[i]
        u = sol.u[i]
        
        # Calculate energies
        KE, PE, total_energy = calculate_energies_exponential(u, t)
        
        # Energy conservation error
        if abs(work_total) > 1e-16
            percent_error = abs((work_total - total_energy) / work_total) * 100.0
            max_error = max(max_error, percent_error)
        else
            percent_error = 0.0
        end
        
        if i % 50 == 1  # Print every 50th point
            @printf("%6.3f %7.4f %9.6f %7.4f %9.6f %9.6f %9.6f %8.4f\n",
                    t, u[1], u[2], u[3], u[4], total_energy, work_total, percent_error)
        end
    end
    
    # Final energy analysis
    final_u = sol.u[end]
    final_KE, final_PE, final_energy = calculate_energies_exponential(final_u, sol.t[end])
    
    println("\nFinal Results:")
    println("=" ^ 40)
    println("Final kinetic energy: $(final_KE)")
    println("Final potential energy: $(final_PE)")
    println("Final total energy: $(final_energy)")
    println("Final work done: $(work_total)")
    println("Energy error: $(abs(work_total - final_energy))")
    println("Maximum error during simulation: $(max_error)%")
    
    if abs(work_total) > 1e-16
        final_percent_error = abs((work_total - final_energy) / work_total) * 100.0
        println("Final percentage error: $(final_percent_error)%")
    else
        println("Final percentage error: 0.0%")
    end
    
    # Display force characteristics at final positions
    println("\nFinal Spring Forces:")
    println("Ground spring force: $(exponential_spring_force(final_u[1], K_VAL, ALPHA_GROUND)) N")
    println("Spring 1-2 force: $(exponential_spring_force(final_u[3] - final_u[1], K_COUPLING, ALPHA_COUPLING)) N")
    println("Spring 2-3 force: $(exponential_spring_force(final_u[5] - final_u[3], K_COUPLING, ALPHA_COUPLING)) N")  
    println("Spring 3-4 force: $(exponential_spring_force(final_u[7] - final_u[5], K_COUPLING, ALPHA_COUPLING)) N")
    
    return sol, work_total, final_energy
end

# Run the simulation
sol, work_total, final_energy = run_exponential_spring_simulation()
