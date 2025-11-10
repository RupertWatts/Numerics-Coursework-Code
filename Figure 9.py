import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.optimize import minimize_scalar

# Function and IC
def f(x):
    return np.exp(-(x - 0.5) ** 2 / (2 * 0.1 ** 2))

# Setup
a = 0
b = 1
t_end = 2
Nx = 50
Nt = 100
delta_x = (b - a) / Nx
delta_t = t_end / Nt
x_FD = np.linspace(a, b, Nx + 1)

# Boundaries and midpoints
boundaries = np.linspace(a, b, Nx + 1)
x_1 = a + delta_x / 2
x_Nx = b - delta_x / 2
u0 = np.linspace(x_1, x_Nx, Nx)

# Mean value over boundaries 
Mean_vals = np.zeros(Nx)
for i in range(Nx):
    result, _ = quad(f, boundaries[i], boundaries[i + 1])
    Mean_vals[i] = result / delta_x

# Flux function
def F(x):
    return 0.5 * x**2

# Set-up vectors
u_FD = f(x_FD)
u1_FD = np.zeros(Nx + 1)
u_G = Mean_vals.copy()

# Times to plot
times_to_plot = [0.2, 1.2]
solutions_FD = {}
solutions_G = {}

for n in range(1, Nt + 1):
    # Finite difference
    for j in range(1, len(x_FD)):
        u1_FD[j] = u_FD[j] - delta_t / (2 * delta_x) * (u_FD[j] ** 2 - u_FD[j - 1] ** 2)
    u1_FD[0] = u1_FD[Nx]  # periodic BCs
    u_FD[:] = u1_FD
    t = n * delta_t

    # Fluxes for Godunov
    Flux = np.zeros(Nx + 1)
    for i in range(Nx + 1):
        if i == 0:
            uL = u_G[Nx - 1]
            uR = u_G[0]
            if uL <= uR:
                temp = minimize_scalar(F, bounds=(uL, uR), method='bounded')
                Flux[i] = np.minimum(F(uL), F(uR))
            else:
                temp = minimize_scalar(lambda x: -F(x), bounds=(uR, uL), method='bounded')
                Flux[i] = -temp.fun

        Flux[Nx] = Flux[0]

        if 1 <= i <= Nx - 1:
            uL = u_G[i - 1]
            uR = u_G[i]
            if uL <= uR:
                temp = minimize_scalar(F, bounds=(uL, uR), method='bounded')
                Flux[i] = temp.fun
            else:
                temp = minimize_scalar(lambda x: -F(x), bounds=(uR, uL), method='bounded')
                Flux[i] = -temp.fun

    u_new = np.zeros(Nx)
    for i in range(Nx):
        u_new[i] = u_G[i] - delta_t / delta_x * (Flux[i + 1] - Flux[i])

    u_G = u_new.copy()

    # Save snapshots if close to target times
    for target in times_to_plot:
        if abs(t - target) < delta_t / 2:
            solutions_FD[target] = u_FD.copy()
            solutions_G[target] = u_G.copy()

# === Plot both times side by side ===
fig, axes = plt.subplots(1, 2, figsize=(12, 4))

for ax, t in zip(axes, times_to_plot):
    if t in solutions_FD:
        # Finite difference
        ax.plot(x_FD, solutions_FD[t], 'r', label='Finite Difference')
        # Godunov (piecewise constant)
        x = np.linspace(0, 1, 1000)
        y = np.zeros_like(x)
        for i in range(Nx):
            mask = (x >= boundaries[i]) & (x < boundaries[i + 1])
            y[mask] = solutions_G[t][i]
        ax.plot(x, y, 'b', label='Godunov')
        ax.set_title(f't = {t}')
        ax.set_xlabel('x')
        ax.set_ylabel('u(x,t)')
        ax.legend()
    else:
        ax.text(0.5, 0.5, f"No data for t={t}", ha='center', va='center')

plt.tight_layout()
plt.show()
