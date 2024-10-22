import sympy as sp

def solve_unadaptive_spread_rate():
    # Define the symbols
    kappa, mu, delta_t = sp.symbols('kappa mu delta_t')
    
    # Define the equation
    equation = sp.Eq(kappa * delta_t + sp.ln(mu * delta_t) + 1, 0)
    
    # Solve the equation for delta_t
    solution = sp.solve(equation, delta_t)
    
    
    return solution

def solve_adaptive_pop():
    # Define the symbols with assumptions
    t, lambda_, mu, i0 = sp.symbols('t lambda mu i0')
    kappa = sp.Symbol('kappa1', positive=True)
    I_1 = i0 * sp.exp(kappa * t)
    r = sp.Function('r')(t)
    ics = {r.subs(t, 0): 1}
    
    # Define the ODE
    ode = sp.Eq(r.diff(t), -lambda_ * mu * r * i0 * sp.exp(kappa * t))
    
    # Solve the ODE
    solution = sp.dsolve(ode, r,ics=ics)
    print("ODE: $$", sp.latex(ode), "$$")
    print("Solution : $$", sp.latex(solution), "$$")
    
    return solution

def solve_adaptive_system():
    # Define the symbols
    t, lambda_, mu, i0 = sp.symbols('t lambda mu i0')
    kappa = sp.Symbol('kappa', positive=True)

    rho_bar_2 = sp.Function('rho_b2')(t)
    rho_bar_3 = sp.Function('rho_b3')(t)
    I_1 = i0 * sp.exp(kappa * t)
    I_2 = sp.Function('I_2')(t)
    I_3 = sp.Function('I_3')(t)
    
    # Define the ODEs
    ode1 = sp.Eq(rho_bar_2.diff(t), -lambda_ * mu * I_1)
    ode2 = sp.Eq(I_2.diff(t), kappa * I_2 + mu * rho_bar_2 * I_1)
    
    # Define initial conditions
    ics = {rho_bar_2.subs(t, 0): 1, I_2.subs(t, 0): 0, rho_bar_3.subs(t, 0): 1, I_3.subs(t, 0): 0}
    
    # Solve the system of ODEs with initial conditions
    solutions = sp.dsolve([ode1, ode2], ics=ics)
    
    print("ODE1 : $$", sp.latex(ode1), "$$")
    print("ODE2 : $$", sp.latex(ode2), "$$")
    for sol in solutions:
        print("Solution : $$", sp.latex(sol), "$$")

    return solutions, ode1, ode2

def solve_mock_model():
    # Define the symbols
    t, mu, i0 = sp.symbols('t mu i0')
    lambda_ = sp.Symbol('lambda', real=True, nonzero=True)
    kappa = sp.Symbol('kappa', positive=True)

    rho_bar_j1 = sp.Function('rho_b_j1')(t)
    I_j = sp.Function('I_j')(t)
    
    # Define the ODEs
    ode1 = sp.Eq(I_j.diff(t), kappa * I_j)
    ode2 = sp.Eq(rho_bar_j1.diff(t), -lambda_ * mu * rho_bar_j1 * I_j)
    f_j1 = rho_bar_j1 * mu * I_j

    # Define initial conditions
    ics = {rho_bar_j1.subs(t, 0): 1, I_j.subs(t, 0): i0}
    
    # Solve the system of ODEs with initial conditions
    solutions = sp.dsolve([ode1, ode2], ics=ics)
    
    # Extract the solved functions
    rho_bar_j1_sol = solutions[0].rhs
    I_j_sol = solutions[1].rhs
    
    # Substitute the solutions into f_j1
    f_j1_sol = f_j1.subs({rho_bar_j1: rho_bar_j1_sol, I_j: I_j_sol})
    
    # Solve for t_f such that f_j1 = 0
    t_f = sp.solve(f_j1_sol, t)
    
    # Integrate f_j1 from t=0 to t=t_f
    integral_f_j1 = sp.integrate(f_j1_sol, (t, 0, t_f))

    print("f_j1 : $$", sp.latex(f_j1), "$$")
    print("Integral of f_j1 : $$", sp.latex(integral_f_j1), "$$","$ t_f =",sp.latex(t_f),"$")
    for sol in solutions:
        print("Solution : $$", sp.latex(sol), "$$")

if __name__ == "__main__":
    # unadaptive_spread_rate = solve_unadaptive_spread_rate()
    # print("Delta t", unadaptive_spread_rate)
    
    # # Print solutions in LaTeX format
    # print("Delta t LaTeX format:", sp.latex(unadaptive_spread_rate))

    # solve_mock_model()
    # solve_adaptive_system()
    solve_adaptive_pop()