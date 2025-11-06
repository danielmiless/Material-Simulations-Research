using DifferentialEquations
using DiffEqCallbacks
using SparseArrays
using LinearAlgebra
using Printf

# Physical parameters
const MASS = 1.0
const K_VAL = 100.0
const K_COUPLING = 100.0
const F_ACTIVE_TIME = 0.05

# Simulation parameters with maximum precision
const T_END = 5.0
const OUTPUT_INTERVAL = 0.01
const ULTRA_HIGH_PRECISION = true

# Ultra-tight tolerances for energy conservation
const REL_TOL = 1e-14
const ABS_TOL = 1e-16

function mass_spring_system!(du, u, p, t)
    x1, v1, x2, v2, x3, v3, x4, v4 = u
    
    external_force = (t <= F_ACTIVE_TIME) ? 10.0 : 0.0
    
    du[1] = v1
    du[2] = (-K_VAL * x1 + K_COUPLING * (x2 - x1) + external_force) / MASS
    du[3] = v2
    du[4] = (K_COUPLING * (x1 - x2) + K_COUPLING * (x3 - x2)) / MASS
    du[5] = v3
    du[6] = (K_COUPLING * (x2 - x3) + K_COUPLING * (x4 - x3)) / MASS
    du[7] = v4
    du[8] = (K_COUPLING * (x3 - x4)) / MASS
    
    return nothing
end

function run_ultra_precision_simulation()
    println("Ultra-High Precision Mass-Spring Simulation")
    println("=" ^ 50)
    
    u0 = zeros(8)
    tspan = (0.0, T_END)
    
    prob = ODEProblem(mass_spring_system!, u0, tspan)
    
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
    
    # Improved work and energy calculations
    work_total = 0.0
    prev_x1, prev_t = 0.0, 0.0
    
    println("\nTime    X1      V1        Energy    Work      %Error")
    println("-" ^ 60)
    
    for i in 1:length(sol.t)
        t = sol.t[i]
        u = sol.u[i]
        
        # Calculate work with improved method
        if i > 1
            x1_curr, x1_prev = u[1], prev_x1
            
            if prev_t < F_ACTIVE_TIME && t <= F_ACTIVE_TIME
                work_increment = 10.0 * (x1_curr - x1_prev)
            elseif prev_t < F_ACTIVE_TIME && t > F_ACTIVE_TIME
                dt_total = t - prev_t
                dt_active = F_ACTIVE_TIME - prev_t
                x1_cutoff = x1_prev + (x1_curr - x1_prev) * (dt_active / dt_total)
                work_increment = 10.0 * (x1_cutoff - x1_prev)
            else
                work_increment = 0.0
            end
            
            work_total += work_increment
        end
        
        # Calculate energy
        KE = 0.5 * MASS * (u[2]^2 + u[4]^2 + u[6]^2 + u[8]^2)
        PE = (0.5 * K_VAL * u[1]^2 + 
              0.5 * K_COUPLING * (u[3] - u[1])^2 +
              0.5 * K_COUPLING * (u[5] - u[3])^2 +
              0.5 * K_COUPLING * (u[7] - u[5])^2)
        total_energy = KE + PE
        
        # Energy error calculation
        if abs(work_total) > 1e-16
            percent_error = abs((work_total - total_energy) / work_total) * 100.0
        else
            percent_error = 0.0
        end
        
        if i % 50 == 1  # Print every 50th point
            @printf("%6.3f %8.5f %9.6f %9.6f %9.6f %8.4f\n",
                    t, u[1], u[2], total_energy, work_total, percent_error)
        end
        
        prev_x1, prev_t = u[1], t
    end
    
    # Final energy analysis
    final_u = sol.u[end]
    final_KE = 0.5 * MASS * (final_u[2]^2 + final_u[4]^2 + final_u[6]^2 + final_u[8]^2)
    final_PE = (0.5 * K_VAL * final_u[1]^2 + 
                0.5 * K_COUPLING * (final_u[3] - final_u[1])^2 +
                0.5 * K_COUPLING * (final_u[5] - final_u[3])^2 +
                0.5 * K_COUPLING * (final_u[7] - final_u[5])^2)
    final_energy = final_KE + final_PE
    
    println("\nFinal Results:")
    println("Final energy: $(final_energy)")
    println("Final work: $(work_total)")
    println("Energy error: $(abs(work_total - final_energy))")
    
    if abs(work_total) > 1e-16
        final_percent_error = abs((work_total - final_energy) / work_total) * 100.0
        println("Final percent error: $(final_percent_error)%")
    else
        println("Final percent error: 0.0%")
    end
    
    return sol, work_total, final_energy
end

# Run the simulation
sol, work_total, final_energy = run_ultra_precision_simulation()
