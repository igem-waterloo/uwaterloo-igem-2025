import numpy as np
from scipy.optimize import least_squares
import matplotlib.pyplot as plt

# ==========================================================
# 1. PDE simulation function
# ==========================================================
def simulate_phage_bacteria_pde(params, Nx, dx, dt, Nt, N0, P0):
    """
    Simulates 1D reaction-diffusion PDE model.
    
    Args:
        params: dict with keys 'k', 'beta', 'm'
        Nx: number of spatial grid points
        dx: spatial step size
        dt: time step size
        Nt: number of time steps
        N0, P0: initial profiles (arrays of length Nx)
    
    Returns:
        N_history, P_history: arrays shape (Nt, Nx)
    """
    k = params['k']
    beta = params['beta']
    m = params['m']
    
    # Fixed parameters
    r = 0.5
    K = 1.0
    D_p = 0.01
    
    # Initialize
    N = N0.copy()
    P = P0.copy()
    
    N_history = [N.copy()]
    P_history = [P.copy()]
    
    for t in range(Nt):
        infection = k * N * P
        dNdt = r * N * (1 - N / K) - infection
        dPdt_reaction = beta * infection - m * P
        
        # Diffusion term
        d2Pdx2 = np.zeros(Nx)
        d2Pdx2[1:-1] = (P[2:] - 2*P[1:-1] + P[:-2]) / dx**2
        
        # Euler update
        N += dt * dNdt
        P += dt * (D_p * d2Pdx2 + dPdt_reaction)
        
        N = np.clip(N, 0, None)
        P = np.clip(P, 0, None)
        
        N_history.append(N.copy())
        P_history.append(P.copy())
    
    return np.array(N_history), np.array(P_history)

# ==========================================================
# 2. Generate synthetic data
# ==========================================================
# Spatial grid
L = 1.0
Nx = 50
dx = L / (Nx - 1)
x = np.linspace(0, L, Nx)

# Time grid
dt = 0.01
Tmax = 2.0
Nt = int(Tmax/dt)
time_points = np.linspace(0, Tmax, 20)
time_indices = np.round(time_points/dt).astype(int)

# Initial conditions
N0 = np.ones(Nx)*0.5
P0 = np.zeros(Nx)
P0[0:5] = 0.2

# True parameters
true_params = {'k':2.0, 'beta':20.0, 'm':0.3}

# Simulate true solution
N_hist_true, P_hist_true = simulate_phage_bacteria_pde(
    true_params, Nx, dx, dt, Nt, N0, P0
)

# Add noise to integrated data
N_integrated_true = N_hist_true[time_indices].sum(axis=1)
P_integrated_true = P_hist_true[time_indices].sum(axis=1)
noise_std = 0.5
N_data = N_integrated_true + np.random.normal(0, noise_std, N_integrated_true.shape)
P_data = P_integrated_true + np.random.normal(0, noise_std, P_integrated_true.shape)

# ==========================================================
# 3. Residuals function
# ==========================================================
def residuals(param_array, Nx, dx, dt, Nt, N0, P0, N_data, P_data, time_indices):
    k, beta, m = param_array
    params = {'k':k, 'beta':beta, 'm':m}
    N_hist, P_hist = simulate_phage_bacteria_pde(params, Nx, dx, dt, Nt, N0, P0)
    N_model = N_hist[time_indices].sum(axis=1)
    P_model = P_hist[time_indices].sum(axis=1)
    return np.concatenate([N_model - N_data, P_model - P_data])

# ==========================================================
# 4. Estimate parameters
# ==========================================================
param_guess = [1.0, 10.0, 0.1]
bounds_lower = [0,0,0]
bounds_upper = [5,50,2]

result = least_squares(
    residuals,
    param_guess,
    args=(Nx, dx, dt, Nt, N0, P0, N_data, P_data, time_indices),
    bounds=(bounds_lower, bounds_upper)
)

k_est, beta_est, m_est = result.x
print("Estimated parameters:")
print(f"k = {k_est:.3f}")
print(f"beta = {beta_est:.3f}")
print(f"m = {m_est:.3f}")

# ==========================================================
# 5. Simulate estimated model
# ==========================================================
params_est = {'k':k_est, 'beta':beta_est, 'm':m_est}
N_hist_est, P_hist_est = simulate_phage_bacteria_pde(
    params_est, Nx, dx, dt, Nt, N0, P0
)

# ==========================================================
# 6. Plot results
# ==========================================================
plt.figure(figsize=(15,5))

# (1) Total bacteria over time
plt.subplot(1,3,1)
plt.scatter(time_points, N_data, color='k', label='Noisy data')
plt.plot(time_points, N_integrated_true, 'b--', label='True')
plt.plot(time_points, N_hist_est[time_indices].sum(axis=1), 'r-', label='Fit')
plt.xlabel('Time')
plt.ylabel('Total Bacteria')
plt.legend()
plt.grid()

# (2) Total phage over time
plt.subplot(1,3,2)
plt.scatter(time_points, P_data, color='k', label='Noisy data')
plt.plot(time_points, P_integrated_true, 'b--', label='True')
plt.plot(time_points, P_hist_est[time_indices].sum(axis=1), 'r-', label='Fit')
plt.xlabel('Time')
plt.ylabel('Total Phage')
plt.legend()
plt.grid()

# (3) Spatial profile at final time
plt.subplot(1,3,3)
plt.plot(x, N_hist_true[-1], 'b--', label='Bacteria True')
plt.plot(x, N_hist_est[-1], 'b-', label='Bacteria Fit')
plt.plot(x, P_hist_true[-1], 'r--', label='Phage True')
plt.plot(x, P_hist_est[-1], 'r-', label='Phage Fit')
plt.xlabel('Space')
plt.ylabel('Density')
plt.legend()
plt.grid()

plt.tight_layout()
plt.show()
