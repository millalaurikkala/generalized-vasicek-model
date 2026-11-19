import numpy as np
import pandas as pd
from scipy import linalg
from sklearn.metrics import mean_absolute_error, mean_squared_error
from statsmodels.tsa.arima.model import ARIMA
from statsmodels.tsa.api import VAR
from statsmodels.tsa.stattools import adfuller, kpss
import matplotlib.pyplot as plt

from vasicek import Vasicek, predict_vasicek

def main():
    plt.rcParams['mathtext.fontset'] = 'stix'
    plt.rcParams['font.family'] = 'STIXGeneral'
    predict = False
    plotting = False
    arma = False
    var = False
    no_change = False
    noise = True
    expanding = False
    # Reading data
    data = pd.read_csv("data.csv")
    data['date'] = pd.to_datetime(data['date'])
    data = data.set_index('date')
    name_1 = "euribor"
    name_2 = "DFF"

    # initialize variables
    r = np.matrix(np.transpose(data.to_numpy())) # data as dim x T matrix
    T = len(data) # length of time series
    s = 10 # upper limit used in estimating coefficient matrices KCPRU-DGS10: 60, or 130 for expanding
    start = int(len(data)*0.8) # start of last 20 % of the data

    # You can use estimate_H.py to estimate H
    H = np.array([0.73, 0.37]) # KCPRU-DGS10: np.array([0.64, 0.57]) 
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

    if predict:
        # predicting last 20 % of the time series with 1-step predictions
        # fit model to first 80 % of data
        if not expanding:
            solver_pred = Vasicek(r[:, :start], H, start)
            solver_pred.estimate_covariance(s)
            solver_pred.estimate_matrices(s)

            # printing parameters
            print("Theta")
            print(pd.DataFrame(solver_pred.theta))
            print("\nb")
            print(pd.DataFrame(solver_pred.b))
            print("\nsigma")
            print(pd.DataFrame(linalg.sqrtm(solver_pred.sigma_product)))

            # predict last 20 %
            p = solver_pred.predict(start, r)

        if expanding:
            p = predict_vasicek(r, H, T, start, s)

        r_1 = np.asarray(r[0,start:])
        r_2 = np.asarray(r[1,start:])
        p_1 = np.asarray(p[0,:])
        p_2 = np.asarray(p[1,:])

        # compare and calculate metrics
        rmse_1 = np.sqrt(mean_squared_error(r_1, p_1))
        rmse_2 = np.sqrt(mean_squared_error(r_2, p_2))
        mae_1 = mean_absolute_error(r_1, p_1)
        mae_2 = mean_absolute_error(r_2, p_2)
        print(f"\nRMSE of component 1: {rmse_1}")
        print(f"RMSE of component 2: {rmse_2}")
        print(f"MAE of component 1: {mae_1}")
        print(f"MAE of component 2: {mae_2}")

    if noise:
        increments, _ = solver.calculate_noise()
        inc_1 = pd.Series(np.squeeze(np.asarray(increments[0, :])))
        inc_2 = pd.Series(np.squeeze(np.asarray(increments[1, :])))
        # Perform ADF test
        adf_1 = adfuller(inc_1)
        adf_2 = adfuller(inc_2)

        #perform KPSS test
        kpss_1 = kpss(inc_1, regression='c')
        kpss_2 = kpss(inc_2, regression='c')

        print(f"p-value {name_1}\n ADF: {adf_1[1]}")
        print(f"p-value {name_2}\n ADF: {adf_2[1]}")

        print(f"p-value {name_1}\n KPSS: {kpss_1[1]}")
        print(f"p-value {name_2}\n KPSS: {kpss_2[1]}")


        fig, ax = plt.subplots(1,2, figsize = (10,4), tight_layout=True)
        ax[0].scatter(data.index, inc_1, color='#f792e8')
        ax[0].axhline(y = 0, color = 'grey', linewidth=0.5, linestyle='-')
        ax[0].set_title(f"Increments of {name_1}", fontsize=18)
        ax[1].scatter(data.index, inc_2, color='#a4bff5')
        ax[1].axhline(y = 0, color = 'grey', linewidth=0.5, linestyle='-')
        ax[1].set_title(f"Increments of {name_2}", fontsize=18)
        ax[0].tick_params(axis='both', which='major', labelsize=16)
        ax[1].tick_params(axis='both', which='major', labelsize=16)

        #plt.savefig("increments_kcpru_dgs10.png", dpi=300)

    if plotting:
        # plot original interest rates and predictions
        fig, ax = plt.subplots(1,1, figsize = (10,5), tight_layout=True)
        ax.plot(data.index, np.transpose(r[0,:]), color='#f792e8', label=f"Euribor")
        ax.plot(data.index, np.transpose(r[1,:]), color='#a4bff5', label=f"{name_2}")

        ax.plot(data.index[start:], np.transpose(p[0,:]), color='magenta', label=f"Euribor prediction")
        ax.plot(data.index[start:], np.transpose(p[1,:]), color='blue', label=f"{name_2} prediction")

        ax.axhline(y = 0, color = 'grey', linewidth=0.5, linestyle='-')
        ax.legend(bbox_to_anchor=(0.63, 0.65),  fontsize=16)
        ax.tick_params(axis='both', which='major', labelsize=16)

        plt.savefig(f"pred_{name_1}_{name_2}.png", dpi=300)

    if arma:
        ser_1 = pd.Series(data[f"{name_1}"])
        ser_2 = pd.Series(data[f"{name_2}"])
        p, q = 1, 1  # You can change these values based on ACF and PACF plots
        arma_ser_1 = ARIMA(ser_1[:T], order=(p, 0, q)).fit()
        arma_ser_2 = ARIMA(ser_2[:T], order=(p, 0, q)).fit()

        predictions_ser_1 = arma_ser_1.predict(start=start, end=T-1)
        predictions_ser_2 = arma_ser_2.predict(start=start, end=T-1)

        # compare and calculate metrics
        rmse_1_ts = np.sqrt(mean_squared_error(ser_1[start:], predictions_ser_1))
        rmse_2_ts = np.sqrt(mean_squared_error(ser_2[start:], predictions_ser_2))
        mae_1_ts = mean_absolute_error(ser_1[start:], predictions_ser_1)
        mae_2_ts = mean_absolute_error(ser_2[start:], predictions_ser_2)
        print("Time series")
        print(f"\nRMSE of component 1: {rmse_1_ts}")
        print(f"RMSE of component 2: {rmse_2_ts}")
        print(f"MAE of component 1: {mae_1_ts}")
        print(f"MAE of component 2: {mae_2_ts}")

    if var:
        ser_1 = pd.Series(data[f"{name_1}"])
        ser_2 = pd.Series(data[f"{name_2}"])
        model = VAR(data)
        results = model.fit(1)

        pred = []
        for i in range(start, T):
            pred.append(results.forecast(data[:i].values, steps=1))

        predictions_ser_1 = list(map(lambda x: x[0][0], pred))
        predictions_ser_2 = list(map(lambda x: x[0][1], pred))

        # compare and calculate metrics
        rmse_1_ts = np.sqrt(mean_squared_error(ser_1[start:], predictions_ser_1))
        rmse_2_ts = np.sqrt(mean_squared_error(ser_2[start:], predictions_ser_2))
        mae_1_ts = mean_absolute_error(ser_1[start:], predictions_ser_1)
        mae_2_ts = mean_absolute_error(ser_2[start:], predictions_ser_2)
        print("VAR")
        print(f"\nRMSE of component 1: {rmse_1_ts}")
        print(f"RMSE of component 2: {rmse_2_ts}")
        print(f"MAE of component 1: {mae_1_ts}")
        print(f"MAE of component 2: {mae_2_ts}")


    if no_change:
        p = np.zeros((2, T-start))
        p = r[:, start-1:-1]

        r_1 = np.asarray(r[0,start:])
        r_2 = np.asarray(r[1,start:])
        p_1 = np.asarray(p[0,:])
        p_2 = np.asarray(p[1,:])

        # compare and calculate metrics
        rmse_1 = np.sqrt(mean_squared_error(r_1, p_1))
        rmse_2 = np.sqrt(mean_squared_error(r_2, p_2))
        mae_1 = mean_absolute_error(r_1, p_1)
        mae_2 = mean_absolute_error(r_2, p_2)
        print("No-change")
        print(f"\nRMSE of component 1: {rmse_1}")
        print(f"RMSE of component 2: {rmse_2}")
        print(f"MAE of component 1: {mae_1}")
        print(f"MAE of component 2: {mae_2}")
if __name__=="__main__":
    main()