# Generalized multivariate Vasicek model

Library for the generalized multivariate Vasicek model studied in [1]:

$dr_t = \Theta(b - r_t)dt + \sigma dX_t$,

where $r$ is the interest rates as a vector, $X$ is a random process and $\Theta$, $b$ and $\sigma$ are the model parameters.

The library includes 
* model fitting
* 1-step predictions
* simulating interest rates from the Vasicek model with an fBM noise
* demos for simulations and real data
* simulation results from [1].


## Installation
If you want to do simulations, install the stochastic library.
Otherwise, standard Python libraries are sufficient.

## Usage
Use the class Vasicek in vasicek.py to initialize your model, fit parameters and do predictions for your data.
To simulate data, use the function simulate_interest_rates in simulate.py.

The file demo_simulation.py contains a demo for simulating data with the Vasicek model with fBM noise with parameters given in vasicek_parameters.py.

The file demo_realdata.py contains a demo for using the library with two real interest rates.
The rates are the 1-month Euribor [2] and the Federal Funds effective rate [3] with daily observations from 4.1.1999 to 6.5.2025.

## References

[1] Ilmonen Pauliina, Laurikkala Milla, Ralchenko Kostiantyn, Sottinen Tommi and Viitasaari Lauri. Data driven modeling of multiple interest rates with generalized Vasicek-type models. In preparation.

[2] Bank of Finland. Euribor* rates and eonia* rate, monthly average. https://www.suomenpankki.fi/en/statistics/data-and-charts/interest-rates/tables/korot_taulukot_en/euribor_korot_long_en/, 2025. Accessed: 6.5.2025.

[3] Federal Reserve Bank of St. Louis. Federal funds effective rate. https://fred.stlouisfed.org/series/DFF, 2025. Accessed: 6.5.2025.