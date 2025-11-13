import matplotlib.pyplot as plt

# IC
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Plot IC on left
ax1.plot(x, f(x), color='b', label='Initial condition', alpha=0.3)
ax2.plot(x, f(x), color='b', label='Initial condition', alpha=0.3)
ax1.set_xlabel('x')
ax1.set_ylabel('u(x,t)')
ax1.set_xlim(0, 1)

# Set up vectors
u = f(x)
t_target = 0.2

for i in range(1, Nt + 1):
    # Finite difference step
    u = finite_difference_step(u, delta_x, delta_t)
    t = i * delta_t

    # Save solution at target time
    if abs(t - t_target) < delta_t / 2:
        u_at_target = u.copy()
        t_actual = t

    # Plot evolution
    if i == 1:
        ax1.plot(x, u, color='r', alpha=0.3, label='Finite difference solution')
    else:
        ax1.plot(x, u, color='r', alpha=0.3)

# Final plot at t_target
ax2.plot(x, u_at_target, color='r', alpha=0.7, label='Finite difference solution')
ax2.set_xlabel('x')
ax2.set_ylabel('u(x,t)')
ax2.set_xlim(0, 1)
ax2.set_title(f't = {t_actual:.2f}')
ax1.set_title(f't ∈ [0,{t_end}]')
ax1.legend(loc='upper left')
ax2.legend(loc='upper left')

plt.tight_layout()
plt.show()

plt.tight_layout()
plt.show()
    
    
