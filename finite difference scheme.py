import numpy as np

# Finite difference step function
def finite_difference_step(u, delta_x, delta_t):
    """
    Perform one finite difference time step for the inviscid Burgers' equation:
        u_t + (u^2 / 2)_x = 0
    using a first-order conservative scheme with periodic boundary conditions.
    """
    Nx = len(u) - 1
    u_next = np.zeros_like(u)

    for j in range(1, Nx + 1):
        u_next[j] = u[j] - delta_t / (2 * delta_x) * (u[j]**2 - u[j - 1]**2)

    # Periodic boundary condition
    u_next[0] = u_next[Nx]
    return u_next


# Function and initial condition
def f(x):
    return 1 * np.exp(-(x - 0.5)**2 / (2 * 0.1**2))


# Setup
a = 0
b = 1
t_end = 2
Nx = 50
Nt = 100
delta_x = (b - a) / Nx
delta_t = t_end / Nt
x = np.linspace(a, b, Nx + 1)  # correct indexing

