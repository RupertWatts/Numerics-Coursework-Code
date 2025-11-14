#Author: Rupert Watts
import numpy as np
import matplotlib.pyplot as plt

# Parameters
a = 0
b = 1
Nx = 50
x = np.linspace(a, b, Nx, endpoint=False)  
delta_x = (b - a) / Nx

# Initial condition
def f(x):
    return np.exp(-(x - 0.5)**2 / (2 * 0.1**2))

f = f(x)

# Inviscid upwind scheme with periodic BC (minimal edits)
def inviscid_upwind(f, delta_t, delta_x, Nx, steps):
    u = f.copy()
    u_hist = np.zeros((Nx, steps+1))
    u_hist[:, 0] = u
    u1 = np.zeros(Nx)
    C = delta_t / delta_x

    for i in range(1, steps+1):
        for j in range(Nx):
            u1[j] = u[j] -  C* (u[j] - u[j-1])
        u1[0] = u1[-1]
        u[:] = u1
        u_hist[:, i] = u.copy()
    return u_hist

# Amplification factor of kth mode over time (normalize each mode by its own initial amplitude)
def amp_factor_over_time(u_hist, Nmodes, k):
    steps = u_hist.shape[1]  # number of time steps + 1
    FMs = np.fft.fft(u_hist, axis=0)        # FFT over spatial axis for each time
    FMs_filtered = FMs[:Nmodes, :]          # keep first Nmodes
    init = np.abs(FMs_filtered[:, 0]).copy()
    init[init == 0] = 1e-14                 # avoid divby-zero
    return np.abs(FMs_filtered[k, :]) / init[k]

# Time parameters
t_final = 0.2
delta_t_vals = np.arange(0.001, 0.1 + 0.001, 0.001)
Nmodes = 6

# Store maximum amplification factors
Max_amp_factor = np.zeros((Nmodes, len(delta_t_vals)))
Var_amp_factor = np.zeros((Nmodes, len(delta_t_vals)))

# Loop over delta_t
for idx, delta_t in enumerate(delta_t_vals):
    steps = int(np.floor(t_final / delta_t))
    if steps < 1:
        continue

    # Solve PDE
    u_hist = inviscid_upwind(f, delta_t, delta_x, Nx, steps)

    # Maximum amplification of each Fourier mode
    for k in range(Nmodes):
        amps_over_time = amp_factor_over_time(u_hist, Nmodes, k)
        Max_amp_factor[k, idx] = np.max(amps_over_time)
        Var_amp_factor[k, idx] = np.var(amps_over_time)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Plot maximum amplification vs Δt
for k in range(Nmodes):
    ax1.plot(delta_t_vals, Max_amp_factor[k, :], label=f'Mode {k}')
    ax2.plot(delta_t_vals, Var_amp_factor[k, :], label=f'Mode {k}')

ax1.set_xlabel(r'$\Delta t$')
ax1.set_ylabel('Maximum amplification factor')
ax1.set_ylim(0, 2.0)
ax1.legend()

ax2.set_xlabel(r'$\Delta t$')
ax2.set_ylabel('Variance in amplification factor')
ax2.legend()

plt.tight_layout()
plt.show()
