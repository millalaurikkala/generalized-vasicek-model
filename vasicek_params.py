import numpy as np

n = 100000 # number of observations
T = 4000 # length of time series

dim = 2 # dimension - number of interest rates to simulate
num_it = 5 # number of iterations - how many paths are simulated
# increasing num_it will increase runtime of demo_simulation.py
# starting with a smaller number to get an idea of your computer's runtime is recommended

H = np.array([0.8, 0.8]) # Hurst indices for components

b = np.matrix([[0],[0]]) # long term mean
theta = np.matrix([[0.1, 0], [0, 0.5]]) # mean reversion + interactions
sigma = np.diag([1,1]) # volatility

s = int(5/T * n) # upper bound for integrals

filename = f"BM_{H[0]}_{num_it}" # path for saving results