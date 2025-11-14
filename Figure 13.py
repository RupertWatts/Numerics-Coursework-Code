# Ewan Strathdee
import numpy as np
import matplotlib.pyplot as plt


def one_timestep(u, Nx, Nt):
    # Calculates the function at the next timestep
    funct = []
    for i in range(Nx+1):
        funct.append(LWM(u, i, Nx, Nt))
    return funct


# Setup
Nx = 50
Nt = 50
dx = 1/Nx
dt = 1/Nt
x = np.linspace(0, 1, Nx + 1)  # correct indices
max_t = 0.2

# Set-up vectors
u = I(x)                
u1 = np.zeros(Nx + 1)    
x_char = np.zeros(Nx + 1)  

# Storage for discretization error
error_L2 = []
time = []

# Calculate the initial condition function
function = IC(Nx)

# Time-stepping loop
for i in range(1, int(Nt*max_t) + 1):
    # Lax-Wendroff method

    times = []
    times.append(function)
    function = one_timestep(function, Nx, Nt)
    
    u[:] = function
    t = i * dt
    time.append(t)
    # Analytical solution along characteristics (Rupert Watts)
    for j in range(len(x)):
        x_char[j] = x[j] + I(x[j]) * t

    u_exact = I(x)  # IC values
    # Interpolate analytic solution onto uniform FD grid (Rupert Watts)
    u_exact_on_grid = np.interp(x, x_char, u_exact)

    # Compute L2 error (Rupert Watts)
    E_L2 = np.sqrt(np.mean((u - u_exact_on_grid)**2))
    error_L2.append(E_L2)


# Right subplot: L2 error over time
plt.plot(np.log(time), np.log(error_L2), marker='o', color='crimson')
plt.xlabel('log(Time t)')
plt.ylabel('log(L2 error)')
plt.title('Lax-Wendroff Accuracy')

plt.show()


