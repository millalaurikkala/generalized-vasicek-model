import numpy as np
from stochastic.processes.continuous import FractionalBrownianMotion

def simulate_OU(l, T, n, randseed=0):
    # simulate an Ornstein-Uhlenbeck process with n observations in the time interval [0, T]
    f = FractionalBrownianMotion(0.5, T, np.random.default_rng(randseed))
    X = f.sample(n-1)
    dt = T/n
    U = np.zeros(n)
    U[0] = 0.1 # idk what this should be but 0 does not work (division by 0)
    for i in range(n-1):
        U[i+1] = U[i] - l * U[i]*dt + (X[i+1] - X[i])
    return U

def simulate_t(l, m, T, n, randseed=0):
    # first simulate m OU-processes
    # randseed info: we need m+1 independent brownian motions so m+1 randseeds
    # it is not possible to take just randseed + i because simulate.py uses that
    # therefore here use randseed*(m+1) + i - it should work as long as dim < m
    U_col = np.zeros((m, n))
    for i in range(m):
        U = simulate_OU(l, T, n, randseed=randseed*(m+1) + i)
        U_col[i, :] = U
    # compute denominator
    S = np.sqrt(1/m * np.sum(np.power(U_col, 2), axis=0))
    # generate brownian motion
    f = FractionalBrownianMotion(0.5, T, np.random.default_rng(randseed*(m+1) + m))
    W = f.sample(n-1)
    # simulate t process starting from 0
    X = np.zeros(n)
    for i in range(n-1):
        X[i+1] = X[i] + (W[i+1] - W[i]) / S[i]

    return X
