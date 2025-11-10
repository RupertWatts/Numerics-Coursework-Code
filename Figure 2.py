import numpy as np
import matplotlib.pyplot as plt

#function and IC
def f(x):
    return 1*np.exp(-(x-0.5)**2/(2*0.1**2))


#setup
a=0
b=1
t_end=2
Nx=50
Nt=100
delta_x=(b-a)/Nx
delta_t=t_end/Nt
x=np.linspace(a, b, Nx+1)#correct ind

#IC
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Plot IC on left
ax1.plot(x, f(x), color='b', label='Initial condition', alpha=0.3)
ax2.plot(x, f(x), color='b', label='Initial condition', alpha=0.3)
ax1.set_xlabel('x')
ax1.set_ylabel('u(x,t)')
ax1.set_xlim(0, 1)

#set-up vectors
u=f(x)
u1=np.zeros(Nx+1)
t_target=0.2

for i in range (1,Nt+1):
    #finite difference method
    for j in range(1,Nx+1):
        u1[j]=u[j]-delta_t/(2*delta_x)*(u[j]**2-u[j-1]**2) 
    u1[0] = u1[Nx]#periodic BCs to conserve mass
    u[:] = u1 #for next time step
    t=0+i*delta_t
    
    if abs(t - t_target) < delta_t/2:
      u_at_target = u.copy()
      t_actual=t

    if i == 1:
        ax1.plot(x, u, color='r', alpha=0.3, label='Finite difference solution')
    else:
        ax1.plot(x, u, color='r', alpha=0.3)

ax2.plot(x, u_at_target, color='r', alpha=0.7, label='Finite difference solution')
ax2.set_xlabel('x')
ax2.set_ylabel('u(x,t)')
ax2.set_xlim(0, 1)
ax2.set_title(f't = {t_actual}')
ax1.set_title(f't $\in$ [0,{t_end}]')
ax1.legend(loc='upper left')
ax2.legend(loc='upper left')

plt.tight_layout()
plt.show()
    
    