import numpy as np
import matplotlib.pyplot as plt


# Function and initial condition (IC)
def f(x):
    return np.exp(-(x-0.5)**2/(2*0.1**2))

# Setup
a = 0
b = 1
t_end = 0.2
Nx = 50
Nt = 10
delta_x = (b - a) / Nx
delta_t = t_end / Nt
x = np.linspace(a, b, Nx + 1)  # correct indices

# Set-up vectors
u = f(x)                
u1 = np.zeros(Nx + 1)    
x_char = np.zeros(Nx + 1)  

# Storage for discretization error
error_L2 = []
time = []


fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

ax1.plot(x, f(x), color='blue', alpha=0.3)


# Time-stepping loop
for i in range(1, Nt + 1):
    # Finite difference method
    for j in range(1, len(x)):
        u1[j] = u[j] - delta_t/(2*delta_x) * (u[j]**2 - u[j-1]**2)
    u1[0] = u1[Nx]  # periodic BC to conserve mass
    u[:] = u1       # update for next time step

    t = i * delta_t
    time.append(t)

    # Analytical solution along characteristics
    for j in range(len(x)):
        x_char[j] = x[j] + f(x[j]) * t

    u_exact = f(x)  # IC values
    # Interpolate analytic solution onto uniform FD grid
    u_exact_on_grid = np.interp(x, x_char, u_exact)

    # Compute L2 error
    E_L2 = np.sqrt(np.mean((u - u_exact_on_grid)**2))
    error_L2.append(E_L2)

    # Plot evolution on left subplot
    if i == 1:
        ax1.plot(x, u, color='r', alpha=0.3, label='Finite difference solution')
        ax1.plot(x_char, f(x), color='b', alpha=0.3, label='Characteristic solution')
    else:
        ax1.plot(x, u, color='r', alpha=0.3)
        ax1.plot(x_char, f(x), color='b', alpha=0.3)


# left subplot

ax1.legend(loc='upper left', fontsize=9)
ax1.set_xlabel('x')
ax1.set_ylabel('u(x,t)')
ax1.set_title(f't $\in$ [0,{t_end}]')

# Right subplot: L2 error over time
ax2.plot(np.log(time), np.log(error_L2), marker='o', color='crimson')
ax2.set_xlabel('log(Time t)')
ax2.set_ylabel('log(L2 error)')

plt.show()
