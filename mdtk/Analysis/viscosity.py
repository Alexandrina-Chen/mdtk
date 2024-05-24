


from .functionals import double_exp_2
from .correlation import autocorrelation
import numpy as np
import scipy

# Constants
KB = 1.38064852e-23  # J/K
NA = 6.02214076e23  # 1/mol
ATM_TO_PA = 1.01325e5  # Pa
A_TO_M = 1e-10  # m
PS_TO_S = 1e-12  # s

def a_m_doubexp(x, A, alpha, tau1, tau2):
    """
    viscosity = A * alpha * tau1 * (1 - exp(-x / tau1)) + A * (1 - alpha) * tau2 * (1 - exp(-x / tau2))
    """
    return A * double_exp_2(x, alpha, 1 - alpha, tau1, tau2)

def a_tau_singleexp(t, a, tau):
    """
    viscosity = a * exp( - (t/tau))
    """
    return a * np.exp(-(t/tau))

def a_b_tau_singleexp(t, a, b, tau):
    """
    viscosity = a * (exp( - (t/tau))) ** b
    """
    return a * (np.exp( - (t/tau))) ** b

class ViscosityAnalysis:
    def __init__(self, timeseries, num_trjs, pxx, pyy, pzz, pxy, pxz, pyz, volume, temperature, method='1'):
        """
        calculate the viscosity from the pressure tensor autocorrelation function, averaged over all trajectories
        The function assumes that the pressure tensor are given in the unit of atm, and the time series is given in the unit of ps.
        
        Parameters
        ----------
        timeseries : array-like of shape (n,)
            time series in the unit of ps
        num_trjs : int
            number of trjs (trajectories)
        pxx : array-like of shape (num_trjs, n)
            xx component of the pressure tensor of each trajectory
        pyy : array-like of shape (num_trjs, n)
            yy component of the pressure tensor of each trajectory
        pzz : array-like of shape (num_trjs, n)
            zz component of the pressure tensor of each trajectory
        pxy : array-like of shape (num_trjs, n)
            xy component of the pressure tensor of each trajectory
        pxz : array-like of shape (num_trjs, n)
            xz component of the pressure tensor of each trajectory
        pyz : array-like of shape (num_trjs, n)
            yz component of the pressure tensor of each trajectory
        volume : float
            volume of the system in the unit of A^3
        temperature : float
            temperature of the system in the unit of K
        method : str, optional, default='1'
            method to calculate the viscosity, 1 is fitting the running integral, 2 is fitting the tail of the autocorrelation function with $a*exp(-t^beta)$
            method 1: see https://doi.org/10.1021/acs.jctc.5b00351, J. Chem. Theory Comput. 2015, 11, 8, 3537–3546
            method 2: see J. Chem. Theory Comput. 2019, 15, 5858−5871
        """
        self.timeseries = np.array(timeseries, dtype=float)
        self.num_trjs = num_trjs
        self.pxx = np.array(pxx, dtype=float)
        self.pyy = np.array(pyy, dtype=float)
        self.pzz = np.array(pzz, dtype=float)
        self.pxy = np.array(pxy, dtype=float)
        self.pxz = np.array(pxz, dtype=float)
        self.pyz = np.array(pyz, dtype=float)
        self.volume = volume
        self.temperature = temperature
        self.method = method
        # check the shape of the input
        for p in [self.pxx, self.pyy, self.pzz, self.pxy, self.pxz, self.pyz]:
            assert p.shape == (self.num_trjs, len(self.timeseries))

    def calculate_viscosity(self, **kwargs):
        """
        
        """
        if self.method == '1':
            return self._calculate_viscosity_1()
        elif self.method == '2':
            return self._calculate_viscosity_2(**kwargs)
        else:
            raise ValueError("method must be 1 or 2")

    def _calculate_viscosity_1(self):
        """
        Returns
        -------
        viscosity : float
            viscosity in the unit of Pa*s
        viscosity_fit : array-like of shape (4,)
            the fitted parameters of the double exponential function, A, alpha, tau1, tau2
            viscosity = A * alpha * tau1 + A * (1 - alpha) * tau2
            viscosity(t) = A * alpha * tau1 * (1 - exp(-t / tau1)) + A * (1 - alpha) * tau2 * (1 - exp(-t / tau2))
        """

        prefix = self.volume / (KB * self.temperature) * ATM_TO_PA**2 * PS_TO_S * A_TO_M**3 # in Pa*s
        pcorr_list = []
        viscosity_list = []
        for i_trj in range(self.num_trjs):
            pxy_corr = autocorrelation(pxy[i_trj])
            pxz_corr = autocorrelation(pxz[i_trj])
            pyz_corr = autocorrelation(pyz[i_trj])
            pxxyy_corr = autocorrelation(pxx[i_trj] - pyy[i_trj])
            pyyzz_corr = autocorrelation(pyy[i_trj] - pzz[i_trj])
            pxxzz_corr = autocorrelation(pxx[i_trj] - pzz[i_trj])
            pcorr = (pxy_corr + pxz_corr + pyz_corr) / 6 + (pxxyy_corr + pyyzz_corr + pxxzz_corr) / 24
            pcorr_list.append(pcorr)
            viscosity = prefix * scipy.integrate.cumtrapz(pcorr, self.timeseries)
            viscosity_list.append(viscosity)


        self.pcorr_array = np.array(pcorr_list)
        self.viscosity_raw_array = np.array(viscosity_list)
        self.viscosity_raw_mean = self.viscosity_raw_array.mean(axis=0)

        popt, pcov = self._fit_viscosity_1(self.timeseries, self.viscosity_raw_mean)
        self.popt = popt
        self.viscosity_fit = a_m_doubexp(self.timeseries, *popt)
        self.viscosity = self.popt[0] * self.popt[1] * self.popt[2] + self.popt[0] * (1 - self.popt[1]) * self.popt[3]

        print(f"Viscosity: {self.viscosity:.6e} Pa*s")
        return self.viscosity, self.viscosity_fit

        

        

    def _fit_viscosity_1(self, times, viscosity):
        self.pcorr_std = self.pcorr_array.std(axis=0)
        # the start of fitting is the first time point where the time is greater than 2.0 ps, because the fluctuations are large at the beginning
        mask_1 = times > 2.0
        # find the first index where the time is greater than 2.0
        start_index = np.argmax(mask_1)

        # the end of fitting is the first time point where the standard deviation is greater than 0.4 times the viscosity
        mask_2 = (times < 1500) * (times > 1000)
        mask_3 = self.pcorr_std > 0.4 * np.mean(viscosity[mask_2])
        end_index = np.argmax(mask_3)

        # fit the viscosity to a double exponential function
        popt, pcov = scipy.optimize.curve_fit(a_m_doubexp, times[start_index:end_index], viscosity[start_index:end_index], p0=[1e-3, 0.5, 1, 10], maxfev=1000000, sigma=self.pcorr_std[start_index:end_index], bounds=(0, [np.inf, 1, np.inf, np.inf]))
        return popt, pcov

    def _calculate_viscosity_2(self, switch_time = 10.0, end_time = 2000.0, my_fit_func = "a_tau_singleexp"):
        """
        Before the switch time, the integral is directly calculated from the autocorrelation function; after the switch time upto end_time, the autocorrelation function is fitted with a single exponential function, and the integral is calculated from the fitted function.
        Typically, the switch time is set to 1.0 ps, and the end time is set to 1500.0 ps.
        """
        prefix = self.volume / (KB * self.temperature) * ATM_TO_PA**2 * PS_TO_S * A_TO_M**3 # in Pa*s
        if my_fit_func == "a_tau_singleexp":
            my_fit = a_tau_singleexp
        elif my_fit_func == "a_b_tau_singleexp":
            my_fit = a_b_tau_singleexp
        else:
            raise ValueError("my_fit_func must be a_tau_singleexp or a_b_tau_singleexp")
        # prepare the autocorrelation function
        pcorr_list = []
        for i_trj in range(self.num_trjs):
            pxy_corr = autocorrelation(self.pxy[i_trj])
            pxz_corr = autocorrelation(self.pxz[i_trj])
            pyz_corr = autocorrelation(self.pyz[i_trj])
            # pxxyy_corr = autocorrelation(self.pxx[i_trj] - self.pyy[i_trj])
            # pyyzz_corr = autocorrelation(self.pyy[i_trj] - self.pzz[i_trj])
            # pxxzz_corr = autocorrelation(self.pxx[i_trj] - self.pzz[i_trj])
            # pcorr = (pxy_corr + pxz_corr + pyz_corr) / 6 + (pxxyy_corr + pyyzz_corr + pxxzz_corr) / 24
            # pcorr = (pxy_corr + pxz_corr + pyz_corr) / 3
            pxx_corr = autocorrelation(self.pxx[i_trj] - (self.pxx[i_trj] + self.pyy[i_trj] + self.pzz[i_trj])/3)
            pyy_corr = autocorrelation(self.pyy[i_trj] - (self.pxx[i_trj] + self.pyy[i_trj] + self.pzz[i_trj])/3)
            pzz_corr = autocorrelation(self.pzz[i_trj] - (self.pxx[i_trj] + self.pyy[i_trj] + self.pzz[i_trj])/3)
            pcorr = (pxy_corr * 2 + pxz_corr * 2 + pyz_corr * 2 + pxx_corr * 4 / 3 + pyy_corr * 4 / 3 + pzz_corr * 4 / 3) / 10
            pcorr_list.append(pcorr)

        self.pcorr_array = np.array(pcorr_list)

        # fit the tail of the autocorrelation function with a single exponential function
        self.tail_parms_array = []
        for i in range(self.num_trjs):
            
            popt, pcov = self._fit_viscosity_2(self.timeseries, self.pcorr_array[i], start_time=switch_time, end_time=end_time, my_fit_func=my_fit_func)
            self.tail_parms_array.append(popt)

        # find the mask
        mask_1 = self.timeseries > switch_time
        switch_index = np.argmax(mask_1)
        mask_2 = self.timeseries > end_time
        end_index = np.argmax(mask_2)

        # store the results
        self.viscosity_array = []
        # integrate the autocorrelation function
        for i in range(self.num_trjs):
            to_fit = list(self.pcorr_array[i][:switch_index]) + list(my_fit(self.timeseries[switch_index:], *self.tail_parms_array[i]))
            self.viscosity_array.append (prefix * scipy.integrate.cumtrapz(to_fit, self.timeseries, initial=0.0))

        self.viscosity_array = np.array(self.viscosity_array)
        self.viscosity_mean = self.viscosity_array.mean(axis=0)
        self.viscosity = self.viscosity_mean[-1]
        self.viscosity_std = self.viscosity_array[:,-1].std()

        # print(f"Viscosity: {self.viscosity:.6e} +- {self.viscosity_std:.6e} Pa*s")

        return self.viscosity, self.viscosity_mean, self.viscosity_array


    def _fit_viscosity_2(self, times, acf, start_time = 10.0, end_time = 2000.0, my_fit_func = "a_tau_singleexp"):
        if my_fit_func == "a_tau_singleexp":
            my_fit = a_tau_singleexp
        elif my_fit_func == "a_b_tau_singleexp":
            my_fit = a_b_tau_singleexp
        else:
            raise ValueError("my_fit_func must be a_tau_singleexp or a_b_tau_singleexp")
        # the start of fitting is the first time point where the time is greater than 1.0 ps
        mask_1 = times > start_time
        # find the first index where the time is greater than start_time
        start_index = np.argmax(mask_1)

        # the end of fitting is the first time point where the time is greater than 1500.0 ps
        mask_2 = times > end_time
        end_index = np.argmax(mask_2)

        # fit the acf to a single exponential function
        print(f"start_index: {start_index} {times[start_index]}, end_index: {end_index} {times[end_index]}")
        # print(acf[start_index:end_index])
        if my_fit_func == "a_tau_singleexp":
            popt, pcov = scipy.optimize.curve_fit(my_fit, times[start_index:end_index], acf[start_index:end_index], p0=[1000.0, 1000.0], maxfev=1000000, bounds=(0, [np.inf, np.inf]))
        elif my_fit_func == "a_b_tau_singleexp":
            popt, pcov = scipy.optimize.curve_fit(my_fit, times[start_index:end_index], acf[start_index:end_index], p0=[10000.0, 1.0, 1000.0], maxfev=1000000, bounds=([100,0.05,100], [np.inf, 1.0, np.inf]))
        # print(popt)
        return popt, pcov
