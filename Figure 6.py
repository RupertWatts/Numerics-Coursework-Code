import numpy as np
import matplotlib.pyplot as plt

# Parameters
a=0
b=1
Nx=50
delta_x=(b-a)/Nx #fixed
x=np.linspace(a, b, Nx+1)

t_final = 0.2#until shock
delta_t_values = [0.01,0.02,0.03,0.05]

#initial condition
def f(x):
    return np.exp(-(x-0.5)**2/(2*0.1**2))

f=f(x)

#upwind scheme up to a certain time step
def inviscid_upwind(f, delta_t, Nx, delta_x, steps):
    #initialising
    u=f.copy()
    u1=np.zeros(Nx+1)
    for i in range (steps):
        #finite difference method
        for j in range(1,Nx+1):
            u1[j]=u[j]-delta_t/(2*delta_x)*(u[j]**2-u[j-1]**2) 
        u1[0] = u1[Nx]#periodic BCs to conserve mass
        u[:] = u1 #for next time step
    return u

# Experiment with different dt

for delta_t in delta_t_values:
    steps = int(t_final/delta_t)
    u = inviscid_upwind(f, delta_t, Nx, delta_x, steps)
    plt.plot(x, u, label=f'$\Delta t$={delta_t}')

plt.legend()
plt.xlabel('x')
plt.ylim(0,1)
plt.ylabel('u')
plt.show()
