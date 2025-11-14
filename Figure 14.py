##### Joseph Winter #####

import numpy as np
import matplotlib.pyplot as plt

mass = []         # store mass values
ke = []           # store kinetic energy values

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

# Compute initial mass and ke
mass.append(np.sum(np.asarray(u))*dx)
ke.append(np.sum(np.asarray(u)**2)*dx)

# calculating the function at each timestep and updates u = function for the next step
function = u
j = 0
times = []
while j < max_t*Nt + 1:
    times.append(function)
    function = timestep(function)
    #Integrating mass and energy to show conservation
    mass.append(np.sum(np.asarray(function))*dx)
    ke.append(np.sum(np.asarray(function)**2)*dx)
    j += 1

# ====== MASS PLOT ======

plt.figure()
plt.plot(np.linspace(0,2,102), mass,alpha=0.7)
plt.xlabel('time')
plt.ylabel('Total Mass')
plt.show()

# ====== KE PLOT ======
plt.figure()
plt.plot(np.linspace(0,2,102), ke,'r',alpha=0.7)
plt.xlabel('time')
plt.ylabel('Total Kinetic Energy')
plt.show()

