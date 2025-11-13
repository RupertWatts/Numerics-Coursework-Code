import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.optimize import minimize_scalar


# === Setup ===
a = 0
b = 1
t_end = 0.2
Nx = 50
Nt = 10
delta_x = (b - a) / Nx
delta_t = t_end / Nt

# === Spatial grid ===
boundaries = np.linspace(a, b, Nx + 1)
x_centers = (boundaries[:-1] + boundaries[1:]) / 2

# === Compute cell-averaged initial condition ===
u_G = np.zeros(Nx)
for i in range(Nx):
    result, _ = quad(f, boundaries[i], boundaries[i + 1])
    u_G[i] = result / delta_x

# === Time integration ===
time = []
error_L2 = []

for n in range(1, Nt + 1):
    # Advance one step using the Godunov function
    u_G = godunov_step(u_G, delta_x, delta_t, F)

    # Analytical (approximate) reference solution via characteristics
    t = n * delta_t
    x_char = x_centers + f(x_centers) * t
    u_exact = f(x_centers)
    u_exact_interp = np.interp(x_centers, x_char % 1, u_exact)  # periodic mapping

    # Compute L2 error
    E_L2 = np.sqrt(np.mean((u_G - u_exact_interp)**2))
    time.append(t)
    error_L2.append(E_L2)

# === Plot log-log error curve ===
plt.figure(figsize=(6, 5))
plt.plot(np.log(time), np.log(error_L2), 'o-', color='crimson')
plt.xlabel('log(Time t)')
plt.ylabel('log(L2 error)')
plt.show()
