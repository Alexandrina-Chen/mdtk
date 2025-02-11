# python 3.9
# last update: 06/22/2024
# Implementation of a discrete autocorrelation function

from copy import deepcopy
import numpy as np


"""
# copy from `https://gist.github.com/anyuzx/0888680ea61090102d706f35bd9b97ea`
# Since there is no real autocorrelation compute function in python and numpy
# I write my own one

import numpy as np

# Suppose data array represent a time series data
# each element represent data at given time
# Assume time are equally spaced
def acf(data):
	mean = np.mean(data)
	var = data.var()
	length = len(data)

	acf_array = []
	for t in np.arange(0,length):
		temp = np.mean((data[:(length-t)] - mean)*(data[t:] - mean))/var
		acf_array.append(temp)

	acf_array = np.array(acf_array)
	return acf_array
    
"""

def acf_bf(data):
    """
    brute force method to calculate autocorrelation of a 1D array of data.
    by definition, the autocorrelation of a time series is the correlation of the time series with itself at different time lags.
    acf(t) = \sum_{i=1}^{N-t} (x_{i} * x_{i+t}) / (N-t)
    
    Parameters
    ----------
    data : array-like
        1D array of data.

    Returns
    -------
    acf : array-like
        autocorrelation of the data.
    
    """
    length = len(data)
    acf = np.zeros(length)
    for t in range(length):
        acf[t] = np.mean((data[:length-t]) * (data[t:]))
    return acf

def autocorrelation(data):
    """
    Calculate the autocorrelation of a 1D array of data.
    """
    extended_data = np.concatenate([data, np.zeros_like(data)])
    fft_data = np.fft.fft(extended_data)
    ifft_data = np.fft.ifft(fft_data * np.conj(fft_data)).real
    autocorr = ifft_data[:len(data)] / np.arange(len(data), 0, -1)
    return autocorr

def correct_intermittency(set_list,intermittency=0):
    if intermittency == 0:
        return set_list

    set_list=deepcopy(set_list)

    for i,set in enumerate(set_list):
        # initially update each frame as seen 0 ago (now)
        seen_frame_ago={s:0 for s in set}
        for j in range(1,intermittency+2):
            for s in seen_frame_ago.keys():
                # no more frames:
                if i + j >= len(set_list):
                    continue

                # if the element is absent now
                if s not in set_list[i + j]:
                    seen_frame_ago[s] +=1
                    continue

                # the element is present
                if seen_frame_ago[s] == 0:
                    # the element was present in the last frame
                    continue

                # the element was absent more times than allowed
                if seen_frame_ago[s] > intermittency:
                    continue

                for k in range(seen_frame_ago[s],0,-1):
                    set_list[i + j - k ].add(s)

                seen_frame_ago[s] = 0
    
    return set_list

class ContinuousCorrelation:
    def __init__(self,set_list,tau_max,window_step=1,intermittency=0) -> None:
        self.set_list=set_list
        self.tau_max=tau_max
        self.window_step=window_step
        self.intermittency=intermittency

    @staticmethod
    def correct_intermittency(set_list,intermittency=0):
        if intermittency == 0:
            return set_list

        set_list=deepcopy(set_list)

        for i,set in enumerate(set_list):
            # initially update each frame as seen 0 ago (now)
            seen_frame_ago={s:0 for s in set}
            for j in range(1,intermittency+2):
                for s in seen_frame_ago.keys():
                    # no more frames:
                    if i + j >= len(set_list):
                        continue

                    # if the element is absent now
                    if s not in set_list[i + j]:
                        seen_frame_ago[s] +=1
                        continue

                    # the element is present
                    if seen_frame_ago[s] == 0:
                        # the element was present in the last frame
                        continue

                    # the element was absent more times than allowed
                    if seen_frame_ago[s] > intermittency:
                        continue

                    for k in range(seen_frame_ago[s],0,-1):
                        set_list[i + j - k ].add(s)

                    seen_frame_ago[s] = 0
        
        return set_list

    def continuous_correlation(self,tell_time=None):
        """Implementation of a discrete autocorrelation function.

        The autocorrelation of a property :math:`x` from a time :math:`t=t_0` to :math:`t=t_0 + \tau`
        is given by:

        .. math::
            C(\tau) = \langle \frac{ x(t_0)x(t_0 +\tau) }{ x(t_0)x(t_0) } \rangle

        where :math:`x` may represent any property of a particle, such as velocity or
        potential energy.

        This function is an implementation of a special case of the time
        autocorrelation function in which the property under consideration can
        be encoded with indicator variables, :math:`0` and :math:`1`, to represent the binary
        state of said property. This special case is often referred to as the
        survival probability (:math:`S(\tau)`). As an example, in calculating the survival
        probability of water molecules within 5 Å of a protein, each water
        molecule will either be within this cutoff range (:math:`1`) or not (:math:`0`). The
        total number of water molecules within the cutoff at time :math:`t_0` will be
        given by :math:`N(t_0)`. Other cases include the Hydrogen Bond Lifetime as
        well as the translocation rate of cholesterol across a bilayer.

        The survival probability of a property of a set of particles is
        given by:

        .. math::
            S(\tau) =  \langle \frac{ N(t_0, t_0 + \tau )} { N(t_0) }\rangle

        where :math:`N(t0)` is the number of particles at time :math:`t_0` for which the feature
        is observed, and :math:`N(t0, t_0 + \tau)` is the number of particles for which
        this feature is present at every frame from :math:`t_0` to :math:`N(t0, t_0 + \tau)`.
        The angular brackets represent an average over all time origins, :math:`t_0`.

        See [Araya-Secchi2014]_ for a description survival probability.

        Parameters
        ----------
        set_lists : list
        List of sets. Each set corresponds to data from a single frame. Each element in a set
        may be, for example, an atom id or a tuple of atoms ids. In the case of calculating the
        survival probability of water around a protein, these atom ids in a given set will be
        those of the atoms which are within a cutoff distance of the protein at a given frame.
        tau_max : int
        The last tau (lag time, inclusive) for which to calculate the autocorrelation. e.g if tau_max = 20,
        the survival probability will be calculated over 20 frames.
        window_step : int, optional
        The step size for t0 to perform autocorrelation. Ideally, window_step will be larger than
        tau_max to ensure independence of each window for which the calculation is performed.
        Default is 1.

        Returns
        --------
        tau_timeseries : list of int
            the values of tau for which the autocorrelation was calculated
        timeseries : list of int
            the autocorelation values for each of the tau values
        timeseries_data : list of list of int
            the raw data from which the autocorrelation is computed, i.e :math:`S(\tau)` at each window.
            This allows the time dependant evolution of :math:`S(\tau)` to be investigated.

        """

        tau_timeseries=list(range(1,self.tau_max+1))
        timeseries_data=[ [] for _ in range(self.tau_max)]

        self.set_list=self.correct_intermittency(self.set_list,intermittency=self.intermittency)

        # calculate autocorrelation
        for t in range(0,len(self.set_list),self.window_step):
            if (tell_time!= None) and t % tell_time == 0:
                print("Already calculate {} frames...".format(t))

            Nt=len(self.set_list[t])

            if Nt == 0:
                continue
            for tau in tau_timeseries:
                if tau + t >= len(self.set_list):
                    break

                Ntau = len(set.intersection(*self.set_list[t:t + tau + 1]))
                timeseries_data[tau-1].append(Ntau/float(Nt))

        timeseries = [np.mean(x) for x in timeseries_data]

        tau_timeseries.insert(0,0)
        timeseries.insert(0,1)

        return tau_timeseries,timeseries,timeseries_data


class IntermittentCorrelation:
    """
    Implementation of a discrete autocorrelation function called intermittent correlation.

    Define, for a given set of particles, the intermittent correlation function as:

    C(tau) = <\sum_{i,j} h_{ij}(t) h_{ij}(t+tau)> / <\sum_{i,j} h_{ij}(t) h_{ij}(t)>

    where h_{ij}(t) is 1 if the property (e.g., h-bond, covalent bond, ion pair, etc.) is present at time t, and 0 otherwise. 
    
    The summation is performed over all possible pairings, ij.
    
    Angular brackets represent an average over many different starting times in the trajectory.

    Intermittent correlation allows the property which were considered broken to be reformed and counted again at a future point in time; 
    
    therefore, measuring the time that a particular hydrogen bonded pair remains in the same vicinity, 
    
    yielding information on the structural relaxation time of the related property.

    """

    def __init__(self,set_list,tau_max,window_step=1) -> None:
        """
        initialize the IntermittentCorrelation object.

        Parameters
        ----------
        set_list : list
            List of sets. Each set corresponds to data from a single frame. Each element in a set
            may be, for example, an atom id or a tuple of atoms ids. In the case of calculating the
            survival probability of water around a protein, these atom ids in a given set will be
            those of the atoms which are within a cutoff distance of the protein at a given frame.
        tau_max : int
            The last tau (lag time, inclusive) for which to calculate the autocorrelation. e.g if tau_max = 20,
            the survival probability will be calculated over 20 frames.
        window_step : int, optional
            The step size for t0 to perform autocorrelation. Ideally, window_step will be larger than
            tau_max to ensure independence of each window for which the calculation is performed.
            Default is 1.
        """
        self.set_list = set_list
        self.tau_max = tau_max
        self.window_step = window_step

    def intermittent_correlation(self):
        """
        Calculate the intermittent correlation of a property of a set of particles.

        Returns
        -------
        tau_timeseries : list of int
            the values of tau for which the autocorrelation was calculated
        timeseries : list of int
            the autocorelation values for each of the tau values
        timeseries_data : list of list of int
            the raw data from which the autocorrelation is computed, i.e :math:`S(\tau)` at each window.
            This allows the time dependant evolution of :math:`S(\tau)` to be investigated.
        """
        tau_timeseries = list(range(1, self.tau_max + 1))
        timeseries_data = [[] for _ in range(self.tau_max)]

        # calculate autocorrelation
        for t in range(0, len(self.set_list), self.window_step):
            Nt = len(self.set_list[t])
            if Nt == 0:
                continue
            for tau in tau_timeseries:
                if tau + t >= len(self.set_list):
                    break

                # calculate the intersection of sets at time `t` and time `t + tau` only, not the entire range
                Ntau = len(set.intersection(self.set_list[t], self.set_list[t + tau]))
                timeseries_data[tau - 1].append(Ntau / float(Nt))

        timeseries = [np.mean(x) for x in timeseries_data]

        tau_timeseries.insert(0, 0)
        timeseries.insert(0, 1)

        return tau_timeseries, timeseries, timeseries_data

    def intermittent_correlation_v2(self):
        """
        Calculate the intermittent correlation of a property of a set of particles.

        Returns
        -------
        tau_timeseries : list of int
            the values of tau for which the autocorrelation was calculated
        timeseries : list of int
            the autocorelation values for each of the tau values
        timeseries_data : list of list of int
            the raw data from which the autocorrelation is computed, i.e :math:`S(\tau)` at each window.
            This allows the time dependant evolution of :math:`S(\tau)` to be investigated.
        """
        tau_timeseries = list(range(0, self.tau_max + 1))
        timeseries_data = [[] for _ in (range(self.tau_max+1)) ]

        # calculate autocorrelation
        for t in range(0, len(self.set_list), self.window_step):
            Nt = len(self.set_list[t])
            timeseries_data[0].append(Nt)
            if Nt == 0:
                continue
            for tau in tau_timeseries:
                if tau + t >= len(self.set_list):
                    break

                # calculate the intersection of sets at time `t` and time `t + tau` only, not the entire range
                Ntau = len(set.intersection(self.set_list[t], self.set_list[t + tau]))
                # timeseries_data[tau - 1].append(Ntau / float(Nt))
                timeseries_data[tau].append(Ntau)

        timeseries = [np.mean(x) for x in timeseries_data]
        timeseries = np.array(timeseries) / timeseries[0]

        

        return tau_timeseries, timeseries, timeseries_data