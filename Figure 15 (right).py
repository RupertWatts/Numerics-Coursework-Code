import numpy as np
import matplotlib.pyplot as plt


def I(x):
    # Initial condition (gaussian)
    return np.exp((-(x - 0.5)**2)/(2*0.1**2))


def f(u):
    return (u**2)/2


def LWM(u, i, N, Nt):
    # Lax-Wendroff Method, returns u_i value at time-step n+1
    dx = 1/N
    dt = 1/Nt
    # adds periodic boundry conditions
    if i == 0:
        il = N
    else:
        il = i - 1

    if i == N:
        ih = 0
    else:
        ih = i + 1
    
    # implement the Lax-Wendroff finite vomume formula
    c = dt/dx
    Ap = (u[i] + u[ih])/2
    Am = (u[i] + u[il])/2
    As.append(Ap)

    T1 = -(c/2)*(f(u[ih]) - f(u[il]))
    T2 = ((c**2)/2)*(Ap*(f(u[ih]) - f(u[i])) - Am*(f(u[i]) - f(u[il])))

    return u[i] + T1 + T2


def timestep(u, N, Nt):
    # Calculates the function at the next timestep
    funct = []
    for i in range(N+1):
        funct.append(LWM(u, i, N, Nt))
    return funct


def Speeds(t, N, Nt, c):
    # Calculate the initial condition function
    u = []
    for j in range(N+1):
        u.append(I(j/N))
    
    # Use Lax-Wendroff to calculate the function after each timestep
    j = 0
    times = []
    while j < max_t*Nt + 1:
        times.append(u)
        u = timestep(u, N, Nt)
        j += 1



# set up
max_t = 0.15

t = max_t*100
N = 50
dx = 1/N

As = []   # List of characteristic speeds at 1+1/2
# record the maximum characteristic speeds for each temporal resolution
maxes = []
for i in range(10, 100):
    As = []
    Speeds(t, N, i, 'r')
    maxes.append(max(As))

fig, ax = plt.subplots()

x = np.linspace(10, 100, 90)
ax.plot(x, maxes)

plt.xlabel('1/Δt')
plt.ylabel('Max(A)')
