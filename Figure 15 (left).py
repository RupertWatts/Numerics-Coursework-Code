##### Ewan Strathdee #####

import numpy as np
import matplotlib.pyplot as plt


def stability(t, Nx, Nt, col):
    # Find the solutuion at a given time (t) for a given N, Nt
    # The plot the solution with a line of colour, col
    dt = round(100*(1/Nt))/100
    
    #Calculate the initial conmdition function
    u = IC(Nx)
    
    # Find the solutions for each timestep
    times = timesteps(u, Nx, Nt)
    
    # plot the solution at time t in colour col
    t_A = int(t*(Nt/100))
    x = np.linspace(0, 1, Nx+1)
    ax.plot(x, times[t_A], color=col, alpha = 0.5, label=f'Δt = {dt}')

    plt.xlabel('x')
    plt.ylabel('u(x, t)')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.3)
    ax.legend(fontsize = 10)

# setup
max_t = 0.15
t = max_t*100
Nx = 50

fig, ax = plt.subplots()

# Plot solutions for given temporal resolutions
stability(t, Nx, 100, 'b')
stability(t, Nx, 50, 'y')
stability(t, Nx, 33, 'g')
stability(t, Nx, 25, 'r')


