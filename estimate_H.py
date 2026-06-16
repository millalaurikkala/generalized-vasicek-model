import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.linear_model import LinearRegression

# Read data
data = pd.read_csv("data.csv")
data['date'] = pd.to_datetime(data['date'])
data = data.set_index('date')
split = 2350 # use only smoother part of data
data = data[split:]
N = len(data)

for rate in data.columns:
    r = pd.Series(data[rate])

    r_qv = []
    thin_range = range(300, 0, -1)
    N_range = []
    for n in thin_range:
        r_thin = r[1::n]
        qv = 0
        N_range.append(len(r_thin))
        for t in range(0, len(r_thin)-2):
            qv += (r_thin.iloc[t+2] - 2*r_thin.iloc[t+1] + r_thin.iloc[t])**2
        r_qv.append(qv)


    # variables for linear regression
    x = np.array(N_range)
    y = r_qv

    # start of stable period, take logarithm
    start = 285
    end = len(N_range)
    x_log = np.log(x[start:end]).reshape(-1, 1)
    y_log = np.log(y[start:end])

    # fit model
    model = LinearRegression()
    model.fit(x_log,y_log)
    coef = model.coef_[0]
    H = (1-coef)/2
    print(f"Hurst index for {rate}: {H}")
