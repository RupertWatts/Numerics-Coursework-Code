import numpy as np
import matplotlib.pyplot as plt

# Function and initial condition
def f(x):
    return np.exp(-(x - 0.5)**2 / (2 * 0.1**2))

# Setup
a = 0
b = 1
t_end = 2
Nx = 50
Nt = 100
delta_x = (b - a) / Nx
delta_t = t_end / Nt
x = np.linspace(a, b, Nx + 1)

# Initial condition
u = f(x)
u1 = np.zeros(Nx + 1)

#set up
mass = []
ke = []
time = []

# Time stepping
for i in range(1, Nt + 1):
    # Finite difference
    for j in range(1, len(x)):
        u1[j] = u[j] - delta_t / (2 * delta_x) * (u[j]**2 - u[j - 1]**2)
    u1[0] = u1[Nx]  # periodic BCs
    u[:] = u1
    t = i * delta_t

    # mass and ke
    total_mass = np.trapezoid(u, x)
    total_ke = np.trapezoid(0.5 * u**2, x)

    mass.append(total_mass)
    ke.append(total_ke)
    time.append(t)

#plotting
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

ax1.plot(time, mass, color='blue', linewidth=2)
ax1.set_xlabel('Time', fontsize=12)
ax1.set_ylabel('Mass', fontsize=12)

ax2.plot(time, ke, color='red', linewidth=2)
ax2.set_xlabel('Time', fontsize=12)
ax2.set_ylabel('Kinetic Energy', fontsize=12)

plt.tight_layout()
plt.show()
