# Simulation parameters

### Increasing T

* Time series length $T = {100, 1075, 2050, 3025, 4000}
* Number of observations is $n = 2.5 T$ so $n = {250, 2687, 5125, 7562, 10000}$
* Dimension d = 2
* Number of realizations N = 1000

* Three cases where both components have the same Hurst index $H = {0.5, 0.6, 0.8}$. Special case: $H = \begin{pmatrix} 0.6 & 0.8 \begin{pmatrix}$
* Upper bound in integral $s = 0.00125 n$

Other parameters correspond to non-diagonal case:
* Long term mean $b = \begin{pmatrix} 1 & 3 \end{pmatrix}$ 
* $\Theta = \begin{pmatrix} 0.5 & 0.1 \\ 0.1 & 0.3 \end{pmatrix}$
* Volatility $\sigma = \begin{pmatrix} 1 & 0 \\ 0 & 2 \end{pmatrix}$ 