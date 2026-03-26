### 3D t distribution simulation parameters

* Time series length T = 4000
* Number of observations n = 10 000 and n = 100 000
* Dimension d = 3
* Number of realizations N = 1000

Noise
* 3-dimensional, each component given by $X_t^{(k)} = \frac{B_t^{(k)}}{\sqrt{\sum_{i=1}^{1000}(Z_i^{(k)})^2}}$
* $B_t$ is standard Brownian motion
* $Z_i^{(k)}$ i.i.d. standard normal variables

Parameters
* Long term mean $b = \begin{pmatrix} 0 & 1 & 3 \end{pmatrix}$ 
* $\Theta = \begin{pmatrix} 0.5 & 0.1 & 0.2 \\ 0.1 & 0.5 & 0.1 \\ 0.2 & 0.1 & 0.5 \end{pmatrix}$
* Volatility $\sigma = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix}$ 
* Covariance used in estimation: $V(t) = t^{2H}I$ for $H=0.5.$
* Upper bound in integral s = 5


### 3D Poisson simulation parameters

* Time series length T = 4000
* Number of observations n = 10 000 and n = 100 000
* Dimension d = 3
* Number of realizations N = 1000

Noise
* Centered Poisson process: $N_k(t) - \lambda_k t$ (three independent components)
* Rates $\lambda = \begin{pmatrix} 0.5 & 1 & 5 \end{pmatrix}$

Parameters
* Long term mean $b = \begin{pmatrix} 0 & 1 & 3 \end{pmatrix}$ 
* $\Theta = \begin{pmatrix} 0.5 & 0.1 & 0.2 \\ 0.1 & 0.5 & 0.1 \\ 0.2 & 0.1 & 0.5 \end{pmatrix}$
* Volatility $\sigma = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix}$  
* Covariance used in estimation: $V(t) = t^{2H}I$ for $H=0.5.$
* Upper bound in integral s = 5