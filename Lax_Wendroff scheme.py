import numpy as np


def I(x):
    # Initial condition (gaussian)
    return np.exp((-(x - 0.5)**2)/(2*0.1**2))


def IC(Nx):
    # generating the initial condition u
    ic = []
    for j in range(Nx+1):
        ic.append(I(j/Nx))
    
    return ic


def f(u):
    return (u**2)/2


As = []   # List of characteristic speeds at i+1/2
def LWM(u, i, Nx, Nt):
    """
    Perform one Lax-Wendroff timestep at a given position, using periodic
    boundry conditions:
    
        Parameters:
        u         : input function
        i         : position number on the function
        Nx        : number of position points the function is calculated for in one unit position
        Nt        : the number of timesteps the function is calculated for in 1 unit time

    Returns:
        u_new     : updated function value at the next timestep for the given position
    """
    dx = 1/Nx
    dt = 1/Nt
    # adds periodic boundry conditions
    if i == 0:
        il = Nx
    else:
        il = i - 1

    if i == Nx:
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

    u_new =  u[i] + T1 + T2
    
    return u_new


max_t = 2
def timesteps(u, Nx, Nt):
    # Calculates the function at each timestep
    times = []
    # Start with the initial condition
    times.append(u)
    # Calculate the function for all timesteps
    j = 0
    while j < max_t*Nt + 1:
        # Find the function value at each position
        funct = []
        for i in range(Nx+1):
            funct.append(LWM(u, i, Nx, Nt))
        # Find the solutions for each timestep
        times.append(funct)
        u = funct
        j += 1

    return times
