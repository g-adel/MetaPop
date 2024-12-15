# using DifferentialEquations
using Plots

function f(du, u, p, t)
    μ = 0.02
    β = 0.25
    λ = 1e10

    # Initialize du to match the structure of u
    du .= 0.0

    du[1] = β * u[1] * (1 - u[1]) + 2μ * u[10] * (u[2] - u[1])
    for i in 2:4
        du[i] = β * u[i] * (1 - u[i]) + μ * u[4 + i] * (u[i-1] - u[i]) + μ * u[9 + i] * (u[i+1] - u[i])
    end
    du[5] = β * u[5] * (1 - u[5]) + μ * u[9] * (u[4] - u[5])

    for i in 1:4
        du[5 + i] = -λ * μ * u[5 + i] * u[i]
    end

    for i in 1:4
        du[9 + i] = -λ * μ * u[9 + i] * u[i + 1]
    end
end

function integrate_f()
    # Initial condition as a one-dimensional array
    u0 = [1e-5, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1]

    # Time span
    tspan = (0.0, 135.0)

    # Define the ODE problem
    prob = ODEProblem(f, u0, tspan)

    # Solve the ODE problem
    sol = solve(prob, Rosenbrock23(), reltol = 1e-8, abstol = 1e-8)

    # Plot the solution
    p=plot(size=(1000,450))
    for i in 1:5
        plot!(p,sol.t,[x[i] for x in sol.u[:]], linewidth =2, title = "Solution to the ODE with one-dimensional array",
            xaxis = "Time (t)", yaxis = "u(t)", label = "Solution")
    end
    plot(p)
end

# Define the SIR model differential equations
function sir_euler(S, I, R, β, γ, dt)
    dS = -β * S * I
    dI = β * S * I - γ * I
    dR = γ * I
    S_new = S + dS * dt
    I_new = I + dI * dt
    R_new = R + dR * dt
    return S_new, I_new, R_new
end

# Function to solve the SIR model using Euler integration
function solve_sir_euler(β, γ, S0, I0, R0, tspan, dt)
    t = tspan[1]:dt:tspan[2]
    S, I, R = [S0], [I0], [R0]
    for _ in t[2:end]
        S_new, I_new, R_new = sir_euler(S[end], I[end], R[end], β, γ, dt)
        push!(S, S_new)
        push!(I, I_new)
        push!(R, R_new)
    end
    return t, S, I, R
end

# Generate the phase space plot
function plot_phase_space(β, γ, tspan, S_range, I_range, dt)
    # Create a grid for the quiver plot
    S_vals = range(S_range[1], S_range[2], length=20)
    I_vals = range(I_range[1], I_range[2], length=20)
    S_grid, I_grid = [S for S in S_vals, I in I_vals], [I for S in S_vals, I in I_vals]
    dS = -β .* S_grid .* I_grid
    dI = β .* S_grid .* I_grid .- γ .* I_grid

    # Plot the quiver plot
    quiver(S_vals, I_vals, quiver=(dS, dI), xlabel="S", ylabel="I", title="Phase Space of SIR Model")

    # Plot solutions for different initial conditions
    # for S0 in 0.05:0.1:0.95
    #     plot!(S, I, label=false)
    # end
    S0=0.95
    t, S, I, R = solve_sir_euler(β, γ, S0, 1-S0, 0.0, tspan, dt)
    plot(t,[S I R],label=["Susceptible" "Infected" "Recovered"],size=(500,500),lw=3)
    # # Draw the diagonal line R = 0
    # plot!([0, 1], [1, 0], linestyle=:dash, color=:red, label="R=0",size=(500,500))

end

# Parameters
β = 0.25
γ = 0.07
tspan = (0.0, 100.0)
S_range = (0.0, 1.0)
I_range = (0.0, 1.0)
dt = 0.1

# Plot the phase space
plot_phase_space(β, γ, tspan, S_range, I_range, dt)