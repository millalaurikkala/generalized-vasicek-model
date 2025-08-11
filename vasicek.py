import numpy as np
from scipy import linalg
from itertools import repeat, starmap
import multiprocessing

class Vasicek:
    def __init__(self, r, H, T=1):
        """
        r - interest rate data
        H - Hurst index/indices for random noise
        T - length of time interval
        """
        self.r = r
        # number of variables, number of observations
        self.d, self.N = np.shape(r)
        self.T = T
        self.delta_t = 1/self.N * T
        self.H = H
        
        self.b = None
        self.gamma = None
        self.model_fitted = False
        
    def calculate_b(self, b=None):
        """
        Estimates b or sets it to given value.
        Inputs:
        b (optional) - value for mean parameter
        """
        if b:
            self.b = b
        else:
            self.b = 1/self.T * np.sum(self.r, axis=1)*self.delta_t
    
    def calculate_gamma(self, b, s):
        """
        Inner function for estimating gamma.
        Much slower than calculate_gamma_numpy, but allows for custom b.
        Inputs:
        b - mean value used in covariance estimation
        s - lag
        """
        cov = np.zeros((self.d, self.d))
        for t in range(0, self.N-s):
            cov += np.outer(self.r[:, t] - b, self.r[:, t+s] - b) * self.delta_t

        return 1/self.T * cov      
        
    def find_gammas_parallel(self, b, s_range):
        """
        Outer function for estimating covariance matrices with different lags.
        Possible to use parallel computations and a custom mean b.
        
        Inputs:
        b - mean value used in covariance estimation
        s_range - lags as indices for which covariance will be estimated
        """

        pool_obj = multiprocessing.Pool()
        gamma_pool = pool_obj.starmap(self.calculate_gamma, zip(repeat(b), s_range))
        pool_obj.close()
        self.gamma = np.array(gamma_pool)
        
    def calculate_gamma_numpy(self, s):
        """
        Inner function for covariance estimation.
        Calculates covariance for lag s using numpy.cov function.
        Inputs:
        s - lag
        Returns:
        dim x dim covariance matrix
        """
        return np.cov(self.r[:, s:], self.r[:, :self.N-s])[0:self.d, self.d:]
    
    def find_gammas_numpy(self, s_range):
        """
        Outer function for estimating covariance matrices with different lags using numpy.
        This is much faster than find_gammas_parallel.
        
        Inputs:
        s_range - lags as indices for which covariance will be estimated.
        """
        self.gamma = np.array([self.calculate_gamma_numpy(s) for s in s_range])
    
    def calculate_B(self, s):
        """
        Calculates the matrix B from [1].
        Approximates the integral with a Riemann sum.
        
        Inputs:
        s - upper limit in the integral (given as index)
        """
        
        B = np.zeros((self.d, self.d))
        """
        for s_k in range(0, s):
            B += (gamma[s_k, :, :] - np.transpose(gamma[s_k, :, :])) * delta_t

        return B
        """
        # list comprehension version (faster)
        # S+1
        self.B = B + sum([(self.gamma[s_k, :, :] - np.transpose(self.gamma[s_k, :, :])) * self.delta_t for s_k in range(0, s+1)])
        
    def calculate_C_inner(self, s, s_k):
        """
        Calculates the inner integral of C for a fixed s_k as an approximation.
        
        Inputs:
        s - upper limit of integral
        s_k - fixed value used in calculations, comes from outer integral
        
        Returns:
        dim x dim matrix, inner integral for s_k
        """
        delta_t_squared = self.delta_t * self.delta_t
        C = np.zeros((self.d, self.d))

        # S+1
        C_p = sum([self.gamma[u_k - s_k, :, :] * delta_t_squared for u_k in range(s_k, s+1)])
        C_n = sum([np.transpose(self.gamma[s_k - u_k, :, :]) * delta_t_squared for u_k in range(s_k)])

        return C_p + C_n + C
    
    def calculate_C_parallel(self, s):
        """
        Outer function for estimating the matrix C in parallel.
        When dimension and the upper limit of the integral are small, this provides little to no increase in computation speed.
        s = upper limit in the integral (given as index)
        """
        # S+1
        s_range = list(range(s+1))

        pool_obj = multiprocessing.Pool()
        self.C = np.array(sum(pool_obj.starmap(self.calculate_C_inner, zip(repeat(s), s_range))))
        pool_obj.close()

    def calculate_C_outer(self, s):
        """
        Outer function for estimating the matrix C.
        s - upper limit in integral.
        """
        s_range = list(range(s+1))
        self.C = np.array(sum(starmap(self.calculate_C_inner, zip(repeat(s), s_range))))
        
    def calculate_sigma(self):
        """
        Calculates the product of sigma*sigma^T.
        For assumptions about the random process, see [1].
        """
        A = np.diag(self.delta_t**(2*self.H))
        
        sAs = np.zeros((self.d, self.d))
        for t in range(0, self.N-1):
            sAs += np.outer(self.r[:, t+1]-self.r[:, t], self.r[:, t+1]-self.r[:, t])

        self.sigma_product = 1/self.N * np.linalg.inv(A) * sAs


    def calculate_cov_diff(self, s):
        """
        Calculates Cov(r_t - r_0).
        s - time difference
        """
        return 2*self.gamma[0, :, :] - self.gamma[s, :, :] - np.transpose(self.gamma[s, :, :])
    
    def solve_are(self, round_C=10):
        """
        Finds estimate for theta by solving a continuous time Riccati equation.
        
        Inputs:
        round_C - if given, the matrices B, C and D are rounded to this precicion to ensure that C is symmetric.
        """
        R = np.eye(self.d)
        if round_C:
            B = np.round(self.B, round_C)
            C = np.round(self.C, round_C)
            D = np.round(self.D, round_C)
        self.theta = linalg.solve_continuous_are(B, linalg.cholesky(C, lower=True), D, R)
              
    def estimate_covariance(self, s=10, b=None, gamma=None):
        """
        Estimates covariances of r (gamma).
        Uses either find_gammas_numpy (faster) or find_gammas_parallel.
        
        It is also possible to use existing covariance estimates and save them to self.gamma with this function.
        
        Inputs:
        s - maximum lag for which covariance is estimated
        b (optional) - mean vector for time series. If not given, b is estimated.
        gamma (optional) - external covariance matrices. Give this if you have existing covariance estimates.
        """
        self.calculate_b(b)
        if gamma is not None:
            self.gamma = gamma
        else:
            # huom s+1
            s_range = list(range(s+1))
            if b is not None:
                self.find_gammas_parallel(b, s_range)
            else:
                self.find_gammas_numpy(s_range)
                
                
    def estimate_matrices(self, s=10):
        """
        Estimates coefficient matrices B, C, D for given time s.
        Solves theta from the CARE equation from them.
        
        Inputs:
        s - upper bound for integrals when estimating B, C, D
        """
        if self.gamma is None:
            print("Initialize gamma first.")
            return
        
        self.calculate_B(s)
        self.calculate_C_outer(s)
        self.calculate_sigma()
        # For now, the variance of noise is assumed to be t^2H for all t
        var_noise = np.power(s*self.delta_t*np.eye(self.d), 2*self.H)
        self.D = var_noise @ self.sigma_product - self.calculate_cov_diff(s)
        
        try:
            self.solve_are()
            self.model_fitted = True
        except np.linalg.LinAlgError:
            print("Could not solve")
        except:
            print("Something else went wrong")
            
    def predict_1step(self, s):
        """
        Do 1-step prediction from time step s (given as index).
        
        Inputs:
        s - time index where prediction starts
        
        Returns:
        prediction as dim x 1 array
        """
        return self.r[:, s] + self.theta @ (self.b-self.r[:, s])*self.delta_t
            
            
    def predict(self, start):
        """
        Do 1-step predictions from time step start until end of time series.
        Does not extrapolate!
        
        Inputs:
        start - time index where prediction starts
        
        Returns:
        pred - dim x (N - start) array of predictions
        """
        if not self.model_fitted:
            print("Fit model first.")
            return
        
        if start > self.N:
            print("Incorrect starting timestep.")
            return
        
        pred = np.matrix(np.zeros((self.d,self.N - start)))
        for s in range(start, self.N):
            pred[:, s-start] = self.predict_1step(s)
            
        return pred
    
    
    def calculate_noise(self):
        """
        Calculate estimates for the noise process and its increments using the model equation and estimated parameters.
        X_0 is assumed to be 0.
        
        Returns:
        increments - dim x N matrix of increments, increments[:, i] = X_i - X_i-1
        noise - dim x N matrix of noise values
        """
        if not self.model_fitted:
            print("Fit model first.")
            return
        sigma_inv = linalg.inv(linalg.sqrtm(self.sigma_product))
        
        increments = np.matrix(np.zeros((self.d, self.N)))
        noise = np.matrix(np.zeros((self.d, self.N)))
        for i in range(0, self.N-1):
            increments[:, i+1] = sigma_inv @ (self.r[:, i] - self.r[:, i+1] + self.theta @ (self.b-self.r[:, i])*self.delta_t)
            noise[:, i+1] = sigma_inv @ (self.r[:, i+1] - self.r[:, i+1] - 
                                         self.theta @ (self.b-self.r[:, i])*self.delta_t) + noise[:, i]
            
        return (increments, noise)
                
                
