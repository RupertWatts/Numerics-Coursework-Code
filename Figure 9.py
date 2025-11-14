#Author: Rupert Watts
import matplotlib.pyplot as plt

# Times to plot
times_to_plot = [0.2, 1.2]
solutions_FD = {}
solutions_G = {}

# === Time loop ===
for n in range(1, Nt + 1):
    # Finite difference update
    for j in range(1, len(x_FD)):
        u1_FD[j] = u_FD[j] - delta_t / (2 * delta_x) * (u_FD[j] ** 2 - u_FD[j - 1] ** 2)
    u1_FD[0] = u1_FD[-1]  # periodic BC
    u_FD[:] = u1_FD
    t = n * delta_t

    # Godunov update
    u_G = godunov_step(u_G, delta_x, delta_t, F)

    # Save snapshots
    for target in times_to_plot:
        if abs(t - target) < delta_t / 2:
            solutions_FD[target] = u_FD.copy()
            solutions_G[target] = u_G.copy()

# === Plot results ===
fig, axes = plt.subplots(1, 2, figsize=(12, 4))

for ax, t in zip(axes, times_to_plot):
    if t in solutions_FD:
        # Finite difference (red)
        ax.plot(x_FD, solutions_FD[t], 'r', label='Finite Difference')

        # Godunov (blue, piecewise constant)
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
