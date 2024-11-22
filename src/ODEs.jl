using DifferentialEquations
using Plots

# Define the function f
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