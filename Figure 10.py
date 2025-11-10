from scipy.integrate import quad
import numpy as np
import matplotlib.pyplot as plt


def f(x):
    return np.exp(-(x - 0.5)**2 / (2 * 0.1**2))

# Flux function
def F(u):
    return 0.5 * u**2

a= 0
b=1
Nx = 50
delta_x = (b - a) / Nx
boundaries = np.linspace(a, b, Nx + 1)
x = (boundaries[:-1] + boundaries[1:]) / 2  # cell centers

# Compute cell-averaged initial condition
u = np.zeros(Nx)
for i in range(Nx):
    result, _ = quad(f, boundaries[i], boundaries[i + 1])
    u[i] = result / delta_x

# Time setup
t_end = 0.2
Nt = 10
delta_t = t_end / Nt

# Storage for L2 error
time = []
error_L2 = []


for n in range(1, Nt + 1):
    Flux = np.zeros(Nx + 1)

    # Compute fluxes (Godunov for inviscid Burgers)
    for i in range(Nx + 1):
        if i == 0:
            uL, uR = u[-1], u[0]  # periodic BC
        elif i == Nx:
            uL, uR = u[-1], u[0]
        else:
            uL, uR = u[i - 1], u[i]

        if uL <= uR:
            if uL >= 0:
                Flux[i] = F(uL)
            elif uR <= 0:
                Flux[i] = F(uR)
            else:
                Flux[i] = 0.0
        else:
            if (uL + uR) / 2 >= 0:
                Flux[i] = F(uL)
            else:
                Flux[i] = F(uR)

    # Update solution
    u_new = u - (delta_t / delta_x) * (Flux[1:] - Flux[:-1])
    u = u_new.copy()

    # Analytical approximation (shift along characteristics)
    t = n * delta_t
    x_char = x + f(x) * t
    u_exact = f(x)
    u_exact_interp = np.interp(x, x_char % 1, u_exact)  # periodic

    # Compute L2 error
    E_L2 = np.sqrt(np.mean((u - u_exact_interp)**2))
    time.append(t)
    error_L2.append(E_L2)


plt.figure(figsize=(6, 5))
plt.plot(np.log(time), np.log(error_L2), 'o-', color='crimson')
plt.xlabel('log(Time t)')
plt.ylabel('log(L2 error)')
plt.title('Godunov Scheme Accuracy')

plt.show()
