import numpy as np
import pandas as pd
from scipy import linalg
import matplotlib.pyplot as plt

from vasicek import Vasicek

# Reading data
data = pd.read_csv("data.csv")
data['date'] = pd.to_datetime(data['date'])
data = data.set_index('date')

# initialize variables
r = np.matrix(np.transpose(data.to_numpy())) # data as dim x T matrix
H = np.array([0.7, 0.7]) # Hurst indices for interest rates
T = len(data) # length of time series
s = 10 # upper limit used in estimating coefficient matrices

# initialize and fit model
solver = Vasicek(r, H, T)
solver.estimate_covariance(s)
solver.estimate_matrices(s)

# printing parameters
print("Theta")
print(pd.DataFrame(solver.theta))
print("\nb")
print(pd.DataFrame(solver.b))
print("\nsigma")
print(pd.DataFrame(linalg.sqrtm(solver.sigma_product)))

# predicting last 20 % of the time series with 1-step predictions
start = int(len(data)*0.8)
p = solver.predict(start)


# plot original interest rates and predictions
fig, ax = plt.subplots(1,1, figsize = (10,5), tight_layout=True)
ax.plot(data.index, np.transpose(r[0,:]), color='#f792e8', label="Euribor")
ax.plot(data.index, np.transpose(r[1,:]), color='#a4bff5', label="Federal Funds Effective Rate")

ax.plot(data.index[start:], np.transpose(p[0,:]), color='magenta', label="Euribor prediction")
ax.plot(data.index[start:], np.transpose(p[1,:]), color='blue', label="Federal Funds Effective Rate prediction")

ax.axhline(y = 0, color = 'grey', linewidth=0.5, linestyle='-')
ax.legend(  fontsize=16)
ax.tick_params(axis='both', which='major', labelsize=16)

plt.savefig("pred.png", dpi=300)