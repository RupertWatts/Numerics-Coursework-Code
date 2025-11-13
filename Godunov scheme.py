import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize_scalar

# === Function and initial condition ===
def f(x):
    return np.exp(-(x - 0.5) ** 2 / (2 * 0.1 ** 2))

# === Flux function ===
def F(x):
    return 0.5 * x**2


# === Godunov scheme (one time step) ===
def godunov_step(u_G, delta_x, delta_t, F):
    """
    Perform one time step of the Godunov finite volume scheme for
    the inviscid Burgers' equation using the given flux function F(u).

    Parameters:
        u_G       : cell-averaged values (1D array)
        delta_x   : spatial step size
        delta_t   : time step size
        F         : flux function (callable)

    Returns:
        u_new     : updated array after one time step
    """
    Nx = len(u_G)
    Flux = np.zeros(Nx + 1)

    # Compute numerical fluxes at cell interfaces
    for i in range(Nx + 1):
        if i == 0:
            uL = u_G[-1]   # periodic BC
            uR = u_G[0]
        elif i == Nx:
            uL = u_G[-1]
            uR = u_G[0]
        else:
            uL = u_G[i - 1]
            uR = u_G[i]

        # Godunov flux selection
        if uL <= uR:
            # monotone interval → flux = min F(u) over [uL, uR]
            res = minimize_scalar(F, bounds=(uL, uR), method='bounded')
            Flux[i] = res.fun
        else:
            # discontinuous interval → flux = max F(u) over [uR, uL]
            res = minimize_scalar(lambda x: -F(x), bounds=(uR, uL), method='bounded')
            Flux[i] = -res.fun

    # Update step
    u_new = np.zeros_like(u_G)
    for i in range(Nx):
        u_new[i] = u_G[i] - (delta_t / delta_x) * (Flux[i + 1] - Flux[i])

    return u_new


# === Setup ===
a = 0
b = 1
t_end = 2
Nx = 50
Nt = 100
delta_x = (b - a) / Nx
delta_t = t_end / Nt
x_FD = np.linspace(a, b, Nx + 1)

# Boundaries and cell centers
boundaries = np.linspace(a, b, Nx + 1)
x_1 = a + delta_x / 2
x_Nx = b - delta_x / 2
u0 = np.linspace(x_1, x_Nx, Nx)

# Mean value over each cell
Mean_vals = np.zeros(Nx)
for i in range(Nx):
    result, _ = quad(f, boundaries[i], boundaries[i + 1])
    Mean_vals[i] = result / delta_x

# Initial conditions
u_FD = f(x_FD)
u1_FD = np.zeros_like(u_FD)
u_G = Mean_vals.copy()
