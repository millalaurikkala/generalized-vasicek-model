import pickle
import numpy as np
import matplotlib.pyplot as plt

from vasicek import Vasicek
from simulate import simulate_interest_rates
import vasicek_params

# arrays for saving estimates
thetas = []
b_hats = []
sigmas = []

for randseed in range(0, vasicek_params.dim*vasicek_params.num_it, vasicek_params.dim):
    # initialize model and fit parameters
    r = simulate_interest_rates(vasicek_params.b, vasicek_params.theta, vasicek_params.sigma, vasicek_params.n,
                                vasicek_params.T, vasicek_params.H, vasicek_params.dim, randseed)
    solver = Vasicek(r, vasicek_params.H, vasicek_params.T)
    solver.estimate_covariance(vasicek_params.s)
    solver.estimate_matrices(vasicek_params.s)
    
    # save estimates
    thetas.append(solver.theta)
    b_hats.append(solver.b)
    sigmas.append(solver.sigma_product)
    print(randseed)
    

# save parameter estimates to files
with open(f'theta_hats_{vasicek_params.filename}.pkl', 'wb') as f:
    pickle.dump(thetas, f)

with open(f'b_hats_{vasicek_params.filename}.pkl', 'wb') as f:
    pickle.dump(b_hats, f)

with open(f'sigma_hats_{vasicek_params.filename}.pkl', 'wb') as f:
    pickle.dump(sigmas, f)

# calculate componentwise differences to true value of theta    
theta_diffs = []
for i, t in enumerate(thetas):
    theta_diffs.append(vasicek_params.theta - t)
theta_diffs = np.array(theta_diffs)

# plot differences as histogram
fig, ax = plt.subplots(vasicek_params.dim, vasicek_params.dim, figsize = (12,8), tight_layout=True)
for j in range(vasicek_params.dim):
    for k in range(vasicek_params.dim):
        ax[j, k].hist(theta_diffs[:, j, k], bins=50, facecolor='#f792e8')
        ax[j, k].set_title(f"Component {j},{k}")
        ax[j, k].set_ylabel("Frequency")
        ax[j,k].grid(True, color ='grey', linestyle ='-.', linewidth = 0.5, alpha = 0.6) 
        ax[j,k].tick_params(bottom=False, left=False)
        
        ax[j,k].spines['top'].set_visible(False)
        ax[j,k].spines['right'].set_visible(False)
        ax[j,k].spines['left'].set_visible(False)
        ax[j,k].spines['bottom'].set_color('grey')
        
        ax[j,k].axvline(x = 0, color = '#be49fc', linestyle='--')

fig.suptitle("Componentwise difference from true value", fontsize=16)

plt.savefig(f"theta_difference_{vasicek_params.filename}.png", dpi=300)