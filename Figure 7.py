import numpy as np
import matplotlib.pyplot as plt




t_final = 0.2  # until shock
delta_t_values = [0.01, 0.02, 0.03, 0.05]

# === Run for different Δt values ===
plt.figure(figsize=(8, 5))

for delta_t in delta_t_values:
    u = f(x)
    steps = int(t_final / delta_t)

    for _ in range(steps):
        u = finite_difference_step(u, delta_x, delta_t)

    plt.plot(x, u, label=f'Δt = {delta_t}')

plt.legend()
plt.xlabel('x', fontsize=12)
plt.ylabel('u(x, t)', fontsize=12)
plt.ylim(0, 1)
plt.tight_layout()
plt.show()

