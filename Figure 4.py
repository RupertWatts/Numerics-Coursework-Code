#Author: Rupert Watts
import numpy as np
import matplotlib.pyplot as plt


# === Initial condition ===
u = f(x)

# === Storage for mass, kinetic energy, time ===
mass = []
ke = []
time = []

# === Time stepping ===
for i in range(1, Nt + 1):
    # Use finite difference function
    u = finite_difference_step(u, delta_x, delta_t)
    t = i * delta_t

    # Compute total mass and kinetic energy
    total_mass = np.trapezoid(u, x)
    total_ke = np.trapezoid(0.5 * u**2, x)

    mass.append(total_mass)
    ke.append(total_ke)
    time.append(t)

# === Plotting ===
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

ax1.plot(time, mass, color='blue', linewidth=2)
ax1.set_xlabel('Time', fontsize=12)
ax1.set_ylabel('Mass', fontsize=12)


ax2.plot(time, ke, color='red', linewidth=2)
ax2.set_xlabel('Time', fontsize=12)
ax2.set_ylabel('Kinetic Energy', fontsize=12)

plt.tight_layout()
plt.show()
