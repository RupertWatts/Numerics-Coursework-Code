import numpy as np
import matplotlib.pyplot as plt


# Function and initial condition (IC)
def f(x):
    return np.exp(-(x-0.5)**2/(2*0.1**2))


def F(u):
    return u**2/2


def LWM(u, i):
    # Lax-Wendroff Method, returns u value at time-step n+1

    # Setting up periodic boundry conditions
    if i == 0:
        il = Nx
    else:
        il = i - 1

    if i == Nx:
        ih = 0
    else:
        ih = i + 1
    
    # Calculating the Lax Wendroff at the next timestep for a given position i
    c = delta_t/delta_x
    Ap = (u[i] + u[ih])/2
    Am = (u[i] + u[il])/2

    T1 = -(c/2)*(F(u[ih]) - F(u[il]))
    T2 = ((c**2)/2)*(Ap*(F(u[ih]) - F(u[i])) - Am*(F(u[i]) - F(u[il])))

    return u[i] + T1 + T2



def timestep(u):
    # Calculates the function at the next timestep
    funct = []
    for i in range(Nx+1):
        funct.append(LWM(u, i))
    return funct


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

# Calculate the initial condition function
function = []
for j in range(Nx+1):
    function.append(f(j/Nx))

# Time-stepping loop
for i in range(1, Nt + 1):
    # Lax-Wendroff method

    times = []
    times.append(function)
    function = timestep(function)
    
    u[:] = function
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


# Right subplot: L2 error over time
plt.plot(np.log(time), np.log(error_L2), marker='o', color='crimson')
plt.xlabel('log(Time t)')
plt.ylabel('log(L2 error)')
plt.title('Lax-Wendroff Accuracy')

plt.show()
