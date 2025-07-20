import numpy as np
from scipy.integrate import odeint
from scipy.optimize import least_squares
import matplotlib.pyplot as plt

# ==========================================================
# 1. Reusable ODE simulation function
# ==========================================================
def simulate_phage_bacteria_odes(params, y0, t):
    """
    Simulates the phage-bacteria ODE model.

    Args:
        params: dict with keys 'k', 's', 'd', 'a', 'r'
        y0: list or array of initial conditions [P0, C0, N0]
        t: array of time points

    Returns:
        solution: array of shape (len(t), 3) with columns [P, C, N]
    """
    k = params['k']
    s = params['s']
    d = params['d']
    a = params['a']
    r = params['r']

    def model(y, t):
        P, C, N = y
        dPdt = k * s * C - d * P
        dCdt = a * N - r * P - s * C
        dNdt = -a * N - r * P
        return [dPdt, dCdt, dNdt]

    solution = odeint(model, y0, t)
    return solution

# ==========================================================
# 2. Generate synthetic data (simulate "observations")
# ==========================================================

# "True" parameters for generating data
true_params = {'k': 50, 's': 1.0, 'd': 0.1, 'a': 0.01, 'r': 0.01}

# Initial conditions
y0 = [1.0, 0.0, 100.0]

# Time points
t_obs = np.linspace(0, 10, 20)

# Simulate true solution
sol_true = simulate_phage_bacteria_odes(true_params, y0, t_obs)

# Add noise to simulate measurement error
noise_std = 5.0
data_noisy = sol_true + np.random.normal(0, noise_std, sol_true.shape)

# ==========================================================
# 3. Define residuals function for parameter estimation
# ==========================================================
def residuals(param_array, t, data_obs, y0):
    # Convert array to dict
    params = {
        'k': param_array[0],
        's': param_array[1],
        'd': param_array[2],
        'a': param_array[3],
        'r': param_array[4]
    }
    # Simulate model with current parameters
    sol = simulate_phage_bacteria_odes(params, y0, t)
    # Flatten residuals into 1D array
    return (sol - data_obs).ravel()

# ==========================================================
# 4. Estimate parameters via least squares
# ==========================================================
# Initial guess
param_guess = [40, 0.8, 0.05, 0.02, 0.02]

# Bounds for parameters
bounds_lower = [0, 0, 0, 0, 0]
bounds_upper = [100, 10, 10, 1, 1]

# Optimization
result = least_squares(
    residuals,
    param_guess,
    args=(t_obs, data_noisy, y0),
    bounds=(bounds_lower, bounds_upper)
)

# Extract estimated parameters
k_est, s_est, d_est, a_est, r_est = result.x

# Print estimated values
print("Estimated parameters:")
print(f"k = {k_est:.3f}")
print(f"s = {s_est:.3f}")
print(f"d = {d_est:.3f}")
print(f"a = {a_est:.4f}")
print(f"r = {r_est:.4f}")

# ==========================================================
# 5. Simulate model with estimated parameters for comparison
# ==========================================================
params_estimated = {
    'k': k_est,
    's': s_est,
    'd': d_est,
    'a': a_est,
    'r': r_est
}
sol_est = simulate_phage_bacteria_odes(params_estimated, y0, t_obs)

# ==========================================================
# 6. Plot results
# ==========================================================
plt.figure(figsize=(15,5))

labels = ['Phage (P)', 'Infected Complex (C)', 'Bacteria (N)']
colors = ['tab:red', 'tab:green', 'tab:blue']

for i in range(3):
    plt.subplot(1,3,i+1)
    # Noisy data points
    plt.scatter(t_obs, data_noisy[:,i], color='k', label='Noisy data')
    # True solution
    plt.plot(t_obs, sol_true[:,i], 'b--', label='True')
    # Estimated solution
    plt.plot(t_obs, sol_est[:,i], color=colors[i], label='Fit')
    plt.xlabel('Time')
    plt.ylabel(labels[i])
    plt.legend()
    plt.grid()
    plt.tight_layout()

plt.show()