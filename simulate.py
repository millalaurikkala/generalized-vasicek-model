import numpy as np
from stochastic.processes.continuous import FractionalBrownianMotion

def simulate_interest_rates(b, theta, sigma, n: int, T: int, H, dim: int, randseed: int):
    """
    Function for simulating interest rate data from the Vasicek model with given parameters and Bm/fBm as noise.
    Inputs:
    b - long term average, array of size dim
    theta - mean reversion rate, dim x dim positive definite matrix
    n - number of datapoints
    T - length of time series
    H - Hurst indices for fBm, array of size dim
    dim - dimension of interest rate vector, should correspond to other parameters
    randseed - random seed for simulating fBm
    
    Returns:
    r - interest rates, dim x n array
    """
    
    X = []
    for d in range(dim): 
        f = FractionalBrownianMotion(H[d], T, np.random.default_rng(randseed+d))
        X.append(f.sample(n-1))
    X = np.matrix(X)

    dt = 1/n*T
    t = np.linspace(0., T, n)  # Vector of times.
    r = np.matrix(np.zeros((dim,n)))
    r[:, 0] = b

    for i in range(n-1):
        # dr = theta*(b-r)dt + sigma*dX
        r[:, i+1] = r[:, i] + theta @ (b-r[:, i])*dt + sigma @ X[:, i+1] - sigma @ X[:, i]
        
    return r