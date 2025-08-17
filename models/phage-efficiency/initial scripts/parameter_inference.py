import numpy as np
import matplotlib.pyplot as plt

# Parameters (common)
r = 0.5          # Bacterial growth rate
K = 1.0          # Carrying capacity
kmax = 1.0       # Max infection rate
Ks = 0.1         # Half-saturation constant
beta = 5         # Burst size (reduced to avoid explosions)
m = 0.1          # Phage decay rate

# Spatial domain
L = 1.0          # Length of biofilm (cm)
Nx = 100         # Number of spatial grid points
dx = L / (Nx-1)
x = np.linspace(0, L, Nx)

# Time domain
dt = 0.0001       # Small time step for stability
Tmax = 1.0       # Total simulation time
Nt = int(Tmax/dt)

# Biofilm density scenarios (different diffusion)
diffusion_scenarios = {
    'Low Diffusion': 0.001,
    'Medium Diffusion': 0.01,
    'High Diffusion': 0.1
}

# Colors for plotting
colors = ['tab:blue', 'tab:green', 'tab:red']

# Initialize figure
plt.figure(figsize=(15,5))

# Simulate each scenario
for i, (label, Dp) in enumerate(diffusion_scenarios.items()):
    # Initialize
    N = np.ones(Nx)*0.2
    P = np.zeros(Nx)
    P[0:5] = 0.5

    # Time loop
    for t in range(Nt):
        # Infection term
        infection = kmax * N * P / (Ks + N)

        # Bacteria update
        dNdt = r*N*(1 - N/K) - infection

        # Phage diffusion
        d2Pdx2 = np.zeros(Nx)
        d2Pdx2[1:-1] = (P[2:] - 2*P[1:-1] + P[:-2]) / dx**2

        # Phage update
        dPdt = Dp*d2Pdx2 + beta*infection - m*P

        # Euler update
        N += dt*dNdt
        P += dt*dPdt

        # Enforce nonnegative
        N = np.maximum(N,0)
        P = np.maximum(P,0)

    # Plot final timepoint
    plt.plot(x, N, color=colors[i], linestyle='-', label=f'{label} Bacteria')
    plt.plot(x, P, color=colors[i], linestyle='--', label=f'{label} Phage')

# Formatting
plt.xlabel('Distance into biofilm')
plt.ylabel('Density')
plt.title('Bacteria and Phage Densities for Different Biofilm Diffusion')
plt.legend()
plt.grid()
plt.tight_layout()
plt.show()
