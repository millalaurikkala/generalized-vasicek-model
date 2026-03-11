import numpy as np
from stochastic.processes.continuous import FractionalBrownianMotion, PoissonProcess

def simulate_interest_rates(b, theta, sigma, n: int, T: int, H, dim: int, randseed: int, process="fbm"):
    """
    Function for simulating interest rate data from the Vasicek model with given parameters and Bm/fBm as noise.
    Inputs:
    b - long term average, array of size dim
    theta - mean reversion rate, dim x dim positive definite matrix
    n - number of datapoints
    T - length of time series
    H - Hurst indices for fBm, array of size dim, or lambdas for Poisson
    dim - dimension of interest rate vector, should correspond to other parameters
    randseed - random seed for simulating fBm
    
    Returns:
    r - interest rates, dim x n array
    """
    
    X = []
    for d in range(dim):
        if process == "fbm":
            f = FractionalBrownianMotion(H[d], T, np.random.default_rng(randseed+d))
            X.append(f.sample(n-1))
        elif process == "poisson":
            rng = np.random.default_rng(seed=randseed+d)
            p = rng.poisson(H[d], size=n)
            s = [np.sum(p[0:i]) for i in range(n)]
            t = np.linspace(0, n, n)
            X.append(s - H[d]*t)
        elif process == "chi2":
            f = FractionalBrownianMotion(H[d], T, np.random.default_rng(randseed+d))
            n_chi2 = 5
            X.append(sum([f.sample(n-1)**2 for i in range(n_chi2)]) - 1)
        elif process == "t":
            f = FractionalBrownianMotion(H[d], T, np.random.default_rng(randseed+d))
            rng = np.random.default_rng(seed=randseed+d)
            X.append(f.sample(n-1)/ (np.linalg.norm(rng.normal(size=1000))))
            

    X = np.matrix(X)

    dt = 1/n*T
    t = np.linspace(0., T, n)  # Vector of times.
    r = np.matrix(np.zeros((dim,n)))
    r[:, 0] = b

    for i in range(n-1):
        # dr = theta*(b-r)dt + sigma*dX
        r[:, i+1] = r[:, i] + theta @ (b-r[:, i])*dt + sigma @ X[:, i+1] - sigma @ X[:, i]
        
    return r