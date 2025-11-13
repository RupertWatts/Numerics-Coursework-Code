from scipy.integrate import quad
from scipy.optimize import minimize_scalar
import numpy as np
import matplotlib.pyplot as plt

# IC
def f(x):
    return np.exp(-(x-0.5)**2/(2*0.1**2))

# Interval
a = 0
b = 1
Nx = 50
delta_x = (b - a) / Nx

# Boundaries of volumes and midpoints
boundaries = np.linspace(a, b, Nx + 1)
x_1 = a + delta_x / 2
x_Nx = b - delta_x / 2
u0 = np.linspace(x_1, x_Nx, Nx)

# Mean value over boundaries
Mean_vals = np.zeros(Nx)
for i in range(Nx):
    result, _ = quad(f, boundaries[i], boundaries[i + 1])
    Mean_vals[i] = result / delta_x

# Godunov IC - as a function
def f_G(x):
    if a <= x <= b:
        for i in range(Nx):
            if boundaries[i] <= x < boundaries[i + 1]:
                return Mean_vals[i]
    else:
        return 'NA'

f_G = np.vectorize(f_G)
x = np.linspace(0, 1, 10000)
y = f_G(x)
plt.plot(x, y, label='initial condition', alpha=0.3)

# Time 
t_end = 2
Nt = 100
delta_t = t_end / Nt

# Flux
def F(x):
    return 0.5 * x**2

# Initialize
u = Mean_vals.copy()
mass = []         # store mass values
time = []         # store time steps

# Compute initial mass
mass.append(np.sum(u**2) * delta_x)
time.append(0)

# Time loop
for n in range(Nt + 1):
    # Fluxes
    Flux = np.zeros(Nx + 1)
    for i in range(Nx + 1):
        if i == 0:
            uL = u[Nx - 1]
            uR = u[0]
            if uL <= uR:
                temp = minimize_scalar(F, bounds=(uL, uR), method='bounded')
                Flux[i] = np.minimum(F(uL), F(uR))
            else:
                temp = minimize_scalar(lambda x: -F(x), bounds=(uR, uL), method='bounded')
                Flux[i] = -temp.fun

        if 1 <= i <= Nx - 1:
            uL = u[i - 1]
            uR = u[i]
            if uL <= uR:
                temp = minimize_scalar(F, bounds=(uL, uR), method='bounded')
                Flux[i] = temp.fun
            else:
                temp = minimize_scalar(lambda x: -F(x), bounds=(uR, uL), method='bounded')
                Flux[i] = -temp.fun  

    Flux[Nx] = Flux[0]  # periodic BCs

    # Update
    u_new = np.zeros(Nx)
    for i in range(Nx):
        u_new[i] = u[i] - delta_t / delta_x * (Flux[i + 1] - Flux[i])

    # Plot solution over time
    def u_next(x):
        if a <= x <= b:
            for i in range(Nx):
                if boundaries[i] <= x < boundaries[i + 1]:
                    return u_new[i]
        else:
            return 'NA'

    t = n * delta_t
    u_next = np.vectorize(u_next)
    x = np.linspace(0, 1, 1000)
    y = u_next(x)
    if n == 1:
        plt.plot(x, y, color='r', alpha=0.3, label='Godunov difference solution')
    else:
        plt.plot(x, y, color='r', alpha=0.3)
    plt.title('t= ' + str(round(t, 3)))
    plt.xlabel('x')
    plt.ylabel('u(x,t)')

    # Store mass at this time step
    mass.append(np.sum(u_new**2*delta_x))
    time.append(t)

    u = u_new

plt.legend(loc='upper left', fontsize=9)
plt.show()

# ====== MASS PLOT ======
plt.figure()
plt.plot(time, mass,alpha=0.7)
plt.xlabel('time')
plt.ylabel('Total Mass')
plt.show()