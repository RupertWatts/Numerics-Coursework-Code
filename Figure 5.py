import numpy as np
import matplotlib.pyplot as plt




# change Setup 
a = 0
b = 1
t_end = 2
Nx = 5000
Nt = 10000
delta_x = (b - a) / Nx
delta_t = t_end / Nt
x = np.linspace(a, b, Nx + 1)


# === Inviscid Burgers Equation ===
u = f(x)
ke_burgers, time = [], []

for i in range(1, Nt + 1):
    # Use modular finite difference step
    u = finite_difference_step(u, delta_x, delta_t)
    t = i * delta_t

    total_ke = np.trapezoid(0.5 * u**2, x)
    ke_burgers.append(total_ke)
    time.append(t)


# === Linear Advection Equation (c = 1) ===
c = 1.0
u = f(x)
u1 = np.zeros_like(u)
ke_adv = []

for i in range(1, Nt + 1):
    for j in range(1, len(x)):
        u1[j] = u[j] - c * delta_t / (2 * delta_x) * (u[j] - u[j - 1])
    u1[0] = u1[-1]  # periodic BC
    u[:] = u1

    total_ke = np.trapezoid(0.5 * u**2, x)
    ke_adv.append(total_ke)


# === Plot ===
plt.figure(figsize=(8, 5))
plt.plot(time, ke_burgers, color='red', linewidth=2, label='Inviscid Burgers')
plt.plot(time, ke_adv, color='blue', linewidth=2, label='Linear Advection (c=1)')
plt.xlabel('Time', fontsize=12)
plt.ylabel('Kinetic Energy', fontsize=12)
plt.legend()
plt.title('Kinetic Energy Evolution: Burgers vs Linear Advection')
plt.tight_layout()
plt.show()
