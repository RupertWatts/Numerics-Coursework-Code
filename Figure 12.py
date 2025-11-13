import numpy as np
import matplotlib.pyplot as plt


def multiplot(lines):
    # Plots the solution at many times on the same graph
    fig, ax = plt.subplots()
    
    spacing = 50   # keeps the line spacing the same for all time resolutions
    # Plotting the solutions at many times
    ax.plot(x, times[1], color='r', label='Lax-Wendroff solution', alpha = 0.5)
    for j in range(int((spacing/Nt)*(lines - 2))):
        ax.plot(x, times[int(((Nt/spacing)*j) + 2)], color='r', alpha = 0.5)
    
    # Plotting the initial condition
    ax.plot(x, u, color='b', label='Initial condition', alpha = 0.8)
    plt.xlabel('x')
    plt.ylabel('u(x, t)')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.3)
    
    ax.legend()


def plot_t(t):
    # Plot a graph of the solution at the timestep number t along with the IC
    fig, ax = plt.subplots()
    ax.plot(x, u, color='b', label='Initial condition', alpha = 0.8)
    ax.plot(x, times[t], color='r', alpha = 0.5, label=f'Solution at t={1/Nt * t}')
    plt.xlabel('x')
    plt.ylabel('u(x, t)')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.3)
        
    ax.legend()
    

# setup
Nx = 50
Nt = 50
x = np.linspace(0, 1, Nx+1)
max_t = 2
u = IC(Nx)

# calculating the function at each timestep
times = timesteps(u, Nx, Nt)

multiplot(Nt*max_t)
plot_t(7)
