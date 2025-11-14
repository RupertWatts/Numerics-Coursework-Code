##### Ewan Strathdee #####

import numpy as np
import matplotlib.pyplot as plt

# setup
max_t = 0.15
Nx = 50
u = IC(Nx)

# record the maximum characteristic speeds for each temporal resolution
maxes = []
for i in range(10, 100):
    Nt = i
    dt = 1/i
    As = []
    timesteps(u, Nx, i)
    maxes.append(max(As))

fig, ax = plt.subplots()

x = np.linspace(10, 100, 90)
ax.plot(x, maxes)

plt.xlabel('1/Δt')
plt.ylabel('Max(A)')


