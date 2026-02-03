import numpy as np
import matplotlib.pyplot as plt

def calculate_derivatives(state_matrix, current_time, initial_human_pop, phi):
    h_naive                 = state_matrix[:, 0]
    h_i1                    = state_matrix[:, 1]
    h_i2                    = state_matrix[:, 2]
    h_r1                    = state_matrix[:, 3]
    h_r2                    = state_matrix[:, 4]
    h_s2                    = state_matrix[:, 5]
    h_s1                    = state_matrix[:, 6]
    h_i12                   = state_matrix[:, 7]
    h_i21                   = state_matrix[:, 8]
    h_r_total               = state_matrix[:, 9]
    vec_s                   = state_matrix[:, 10]
    vec_i1                  = state_matrix[:, 11]
    vec_i2                  = state_matrix[:, 12]

    total_vectors = vec_s + vec_i1 + vec_i2
    
    # Constants
    mu_h, gamma, sigma, mu_v = 1/65, 52, 2, 36.5
    biting_rate = 2 * mu_v
    beta = 2 * gamma
    
    # Seasonality
    season_freq = 2 * np.pi * 1 
    v_birth = mu_v * (1 + 0.4 * np.cos(season_freq * current_time))

    # Force of infection (Vectors to Humans)
    lambda_h1 = (beta / total_vectors) * vec_i1
    lambda_h2 = (beta / total_vectors) * vec_i2

    # Derivatives
    d_h_naive = -h_naive * (lambda_h1 + lambda_h2) + mu_h * (initial_human_pop - h_naive)
    d_h_i1 = h_naive * lambda_h1 - (gamma + mu_h) * h_i1
    d_h_i2 = h_naive * lambda_h2 - (gamma + mu_h) * h_i2
    d_h_r1 = gamma * h_i1 - (sigma + mu_h) * h_r1
    d_h_r2 = gamma * h_i2 - (sigma + mu_h) * h_r2
    d_h_s2 = sigma * h_r1 - h_s2 * lambda_h2 - mu_h * h_s2
    d_h_s1 = sigma * h_r2 - h_s1 * lambda_h1 - mu_h * h_s1
    d_h_i12 = h_s2 * lambda_h2 - (gamma + mu_h) * h_i12
    d_h_i21 = h_s1 * lambda_h1 - (gamma + mu_h) * h_i21
    d_h_r_total = gamma * (h_i12 + h_i21) - mu_h * h_r_total

    # Vector dynamics
    inf_h1 = h_i1 + phi * h_i21
    inf_h2 = h_i2 + phi * h_i12
    
    d_vec_s = v_birth * total_vectors - (biting_rate / initial_human_pop) * vec_s * (inf_h1 + inf_h2) - mu_v * vec_s
    d_vec_i1 = (biting_rate / initial_human_pop) * vec_s * inf_h1 - mu_v * vec_i1
    d_vec_i2 = (biting_rate / initial_human_pop) * vec_s * inf_h2 - mu_v * vec_i2

    return np.array([
        d_h_naive, 
        d_h_i1,
        d_h_i2,
        d_h_r1,
        d_h_r2,
        d_h_s2,
        d_h_s1,
        d_h_i12,
        d_h_i21,
        d_h_r_total,
        d_vec_s,
        d_vec_i1,
        d_vec_i2]).T

def runge_kutta_step(state, t, dt, total_pop, phi):
    k1 = calculate_derivatives(state, t, total_pop, phi)
    k2 = calculate_derivatives(state + 0.5 * dt * k1, t + 0.5 * dt, total_pop, phi)
    k3 = calculate_derivatives(state + 0.5 * dt * k2, t + 0.5 * dt, total_pop, phi)
    k4 = calculate_derivatives(state + dt * k3, t + dt, total_pop, phi)
    return np.maximum(state + (dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4), 0)

def run_simulation(initial_state, years, phi, samples=1, seed=None):
    if seed is not None: np.random.seed(seed)
    dt = 1/365
    steps = int(years * 365)
    state = np.tile(initial_state, (samples, 1)).astype(float)
    
    if samples > 1:
        state = state * np.random.uniform(0.95, 1.05, state.shape)
    
    initial_human_pop = np.sum(state[:, :10], axis=1)
    history = np.zeros((steps, samples, 13))
    curr_t = 0
    for s in range(steps):
        state = runge_kutta_step(state, curr_t, dt, initial_human_pop, phi)
        history[s] = state
        curr_t += dt
    return history, np.linspace(0, years, steps)

def plot_all_variables(history, time, phi_val):
    var_names = ["S (Naive)", "I1 (Pri)", "I2 (Pri)", "R1 (Cross)", "R2 (Cross)", 
                 "S2 (Susc to 2)", "S1 (Susc to 1)", "I12 (Sec)", "I21 (Sec)", 
                 "R (Final)", "Sv (Vec)", "V1 (Vec)", "V2 (Vec)"]
    
    fig, axes = plt.subplots(4, 4, figsize=(20, 15))
    axes_flat = axes.flatten()
    
    for i in range(13):
        ax = axes_flat[i]
        # Plotting the 5 samples
        ax.plot(time, history[:, :, i], linewidth=0.8, alpha=0.7)
        ax.set_title(var_names[i], fontsize=10, fontweight='bold')
        ax.grid(True, alpha=0.3)
        if i >= 9: ax.set_xlabel("Years")

    for i in range(13, 16): fig.delaxes(axes_flat[i])
    
    plt.suptitle(f"Multi-Strain Dynamics with ADE (phi = {phi_val})", fontsize=18)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    
    filename = f"results_phi_{str(phi_val).replace('.', '_')}.png"
    plt.savefig(filename, dpi=150)
    print(f"Saved plot: {filename}")
    plt.close()

# --- EXPERIMENT START ---
base_init = [700, 200, 100, 0, 0, 0, 0, 0, 0, 0, 9000, 500, 500]

print("Phase 1: Reaching equilibrium (300 years)...")
hist_eq, _ = run_simulation(base_init, years=300, phi=2.5, samples=1)
final_state_eq = hist_eq[-1, 0, :]

print("Phase 2: Running phi=0.8 (100 years, 5 samples)...")
hist_08, time_08 = run_simulation(final_state_eq, years=100, phi=0.8, samples=5, seed=123)
plot_all_variables(hist_08, time_08, 0.8)

print("Phase 3: Running phi=2.5 (100 years, 5 samples)...")
hist_25, time_25 = run_simulation(final_state_eq, years=100, phi=2.5, samples=5, seed=123)
plot_all_variables(hist_25, time_25, 2.5)

print("Experiment complete")
