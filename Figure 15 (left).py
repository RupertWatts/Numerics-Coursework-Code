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
    T1 = -(c/2)*(f(u[ih]) - f(u[il]))
    T2 = ((c**2)/2)*(Ap*(f(u[ih]) - f(u[i])) - Am*(f(u[i]) - f(u[il])))

    return u[i] + T1 + T2


def timestep(u, N, Nt):
    # Calculates the function at the next timestep
    funct = []
    for i in range(N+1):
        funct.append(LWM(u, i, N, Nt))
    return funct


def stability(t, N, Nt, col):
    # Find the solutuion at a given t for a given N, Nt
    dt = round(100*(1/Nt))/100
    
    #Calculate the initial conmdition function
    u = []
    for j in range(N+1):
        u.append(I(j/N))
    
    # Find the solutions for each timestep
    j = 0
    times = []
    while j < max_t*Nt + 1:
        times.append(u)
        u = timestep(u, N, Nt)
        j += 1
    
    # plot the solution at time t
    t_A = int(t*(Nt/100))
    x = np.linspace(0, 1, N+1)
    ax.plot(x, times[t_A], color=col, alpha = 0.5, label=f'Δt = {dt}')

    plt.xlabel('x')
    plt.ylabel('u(x, t)')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.3)
    ax.legend(fontsize = 10)

# set-up
max_t = 0.15

t = max_t*100
N = 50
dx = 1/N

fig, ax = plt.subplots()

# Plot solutions for given temporal resolutions
stability(t, N, 100, 'b')
stability(t, N, 50, 'y')
stability(t, N, 33, 'g')
stability(t, N, 25, 'r')
