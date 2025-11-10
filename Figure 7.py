import numpy as np
import matplotlib.pyplot as plt

# Domain
a=0 
b=1
Nx = 50
x = np.linspace(a, b, Nx+1, endpoint=False)

# Original function
def f(x):
    return np.exp(-(x - 0.5)**2 / (2 * 0.1**2))

f_true = f(x)

# Fourier transform
f_FT = np.fft.fft(f_true) #coef.s of modes
k = np.arange(Nx)# mode numbers

# first 5 modes
N_modes = 6
f_FT_filtered = np.zeros_like(f_FT, dtype=complex)
f_FT_filtered[:N_modes+1] = f_FT[:N_modes+1] #up to 5th mode
f_FT_filtered[-N_modes:] = f_FT[-N_modes:] # last 5 modes

# Reconstruct truncated signal
f_approx = np.real(np.fft.ifft(f_FT))

# Plot function and spectrum
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Function vs Fourier approximiation
ax1.plot(x, f_true, 'k', lw=2, label='$f(x)$')
ax1.plot(x, f_approx, 'r--', lw=2, label=f'First {N_modes} modes')
ax1.set_xlabel('$x$')
ax1.legend()

# power spectrum
ax2.plot(np.abs(f_FT[:int(Nx/2)]), 'b-', lw=2)
ax2.set_xlabel('Mode number $k$')
ax2.set_xlim(0,10)
ax2.set_ylabel('Power')

plt.show()
