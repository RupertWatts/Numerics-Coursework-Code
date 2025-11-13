import numpy as np
import matplotlib.pyplot as plt


mass = []         # store mass values


def I(x):
    # Initial condition (gaussian)
    return np.exp((-(x - 0.5)**2)/(2*0.1**2))


def f(u):
    return (u**2)/2


def LWM(u, i):
    # Lax-Wendroff Method, returns u_i value at time-step n+1
   
    # adds periodic boundry conditions
    if i == 0:
        il = N
    else:
        il = i - 1

    if i == N:
        ih = 0
    else:
        ih = i + 1
   
    # implement the Lax-Wendroff finite difference formula
    c = dt/dx
    Ap = (u[i] + u[ih])/2
    Am = (u[i] + u[il])/2

    T1 = -(c/2)*(f(u[ih]) - f(u[il]))
    T2 = ((c**2)/2)*(Ap*(f(u[ih]) - f(u[i])) - Am*(f(u[i]) - f(u[il])))

    return u[i] + T1 + T2



def timestep(u):
    # Calculate the all x values at a given for the (n+1)th time
    funct = []
    for i in range(N+1):
        funct.append(LWM(u, i))
    return funct


def animation(g, e):
    # updates the graph live at each new time step
    for j in range(int(max_t*Nt + 1)):
        '''
        fig, ax = plt.subplots()
       
        ax.plot(x, times[j], color='r', label='Lax-Wendroff solution', alpha = 0.5)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1.3)
        ax.set_title(f't = {j*dt}')
        ax.legend()
        '''

    return g, e

def multiplot(lines):
    # Plots the solution at many times on the same graph
    fig, ax = plt.subplots()

    ax.plot(x, times[1], color='r', label='Lax-Wendroff solution', alpha = 0.5)
    spacing = 50   # keeps the graph line spacing the same for all time resolutions
    for j in range(int((spacing/Nt)*(lines - 2))):
        ax.plot(x, times[int(((Nt/spacing)*j) + 2)], color='r', alpha = 0.5)
        #print(f't = {(j+2)*0.02}')

    ax.plot(x, u, color='b', label='Initial condition', alpha = 0.8)
    plt.xlabel('x')
    plt.ylabel('u(x, t)')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.3)
   
    ax.legend()


def plot_t(t):
    # Plots a graph of the solution at a fixed time, along with the IC
    # Only the timestep number t is plotted
    fig, ax = plt.subplots()
    ax.plot(x, u, color='b', label='Initial condition', alpha = 0.8)
    ax.plot(x, times[t], color='r', alpha = 0.5, label=f'Solution at t={1/Nt * t}')
    plt.xlabel('x')
    plt.ylabel('u(x, t)')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.3)
       
    ax.legend()


N = 50   # The number of points x is calculated at
Nt = 50   # The number of points times calculated at
dx = 1/N   # spatial resolution
dt = 1/Nt   # time resolution

max_t = 2   # end time

u = []
# generating the initial condition u in descrete form
for j in range(N+1):
    u.append(I(j/N))

x = np.linspace(0, 1, N+1)

mass.append(np.sum(np.asarray(u)**2)*dx)

# calculating the function at each timestep and updates u = function for the next step
function = u
j = 0
times = []
while j < max_t*Nt + 1:
    times.append(function)
    function = timestep(function)
    #Integrating mass and energy to show conservation
    mass.append(np.sum(np.asarray(function)**2)*dx)
    j += 1



#multiplot(Nt*max_t)
#plot_t(6)

plt.figure()
plt.plot(np.linspace(0,2,102), mass,alpha=0.7,color='r')
plt.xlabel('time')
plt.ylabel('Total Kinetic Energy')
plt.show()
