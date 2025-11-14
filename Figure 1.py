#Author: Rupert Watts
import numpy as np
import matplotlib.pyplot as plt

#initial condition
def f(x):
    return np.exp(-(x - 0.5)**2 / (2 * 0.1**2))


# Discretisation
N = 200
x = np.linspace(0, 1, N)
u = f(x)

# vectors for velocity visualisation
x_vec = np.linspace(0.5, 0.8, 10)
y_vec = f(x_vec)
u_vec = y_vec
v_vec = np.zeros_like(x_vec)

# Shock formation (first instance)
t_candidates = (x[1:] - x[:-1]) / (u[:-1] - u[1:])#all intersections
valid = (t_candidates > 0) & np.isfinite(t_candidates)#choose only positive and finite ones
i_first = np.argmin(np.where(valid, t_candidates, np.inf)) 
t_shock = t_candidates[i_first]
x_shock = x[i_first] + u[i_first] * t_shock

#initialisation for further shocks
shock_t = [t_shock] 
shock_x = [x_shock]
iL, iR = i_first, i_first + 1 #indicies of first characteristics that intersect

#shock propagation
for i in range(int(N/2)): #try to find N/2 points of intersection (i not used)
    if iL == 0 or iR == N - 1:
        break #stop if we start to go out of bounds
    iL_next = iL - 1
    iR_next = iR + 1
    xL, uL = x[iL_next], u[iL_next]
    xR, uR = x[iR_next], u[iR_next]
    t_next = (xR - xL) / (uL - uR)
    print(t_next)
    if t_next <= shock_t[-1]: #stops when weve reached
        break
    x_next = xL + uL * t_next
    shock_t.append(t_next)
    shock_x.append(x_next)
    iL, iR = iL_next, iR_next #for next iteration

#Set-up plots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12,5))

# Left: function and velocities
ax1.plot(x, u, color='blue', lw=1.5)
ax1.quiver(x_vec, y_vec, u_vec, v_vec, scale=1.5, color='red', width=0.003)
ax1.set_xlabel('x')
ax1.set_ylabel('f(x)')
ax1.set_xlim(0, 1.5)
ax1.set_ylim(0, 1.1)

# characteristics
t_vals = np.linspace(0,1,200)
for i in range(len(x)):
    xi = x[i]
    ui = u[i]
    ax2.plot(xi + ui * t_vals, t_vals, color='steelblue', lw=0.4)

#shock and trajectory
ax2.scatter(shock_x[0], shock_t[0], color='crimson', s=100,
            label=f'Shock forms at (x={shock_x[0]:.3f}, t={shock_t[0]:.3f})',zorder=5)
ax2.plot(shock_x, shock_t, color='crimson', lw=3, label='Shock trajectory')
ax2.set_xlabel('x'); ax2.set_ylabel('t')
ax2.set_xlim(0, 1); ax2.set_ylim(0, 1)
ax2.legend()

# Show legend in right plot
ax2.legend()

plt.show()



