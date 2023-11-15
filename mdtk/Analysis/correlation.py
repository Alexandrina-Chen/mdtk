# python 3.9
# last update: 07/22/2022
# Implementation of a discrete autocorrelation function

from copy import deepcopy
import numpy as np

class Correlation:
    def __init__(self,list_of_set,tau_max,window_step=1,intermittency=0) -> None:
        self.SetList=list_of_set
        self.MaxTau=tau_max
        self.WindowStep=window_step
        self.Intermittency=intermittency

    @staticmethod
    def CorrectIntermittency(SetList,intermittency=0):
        if intermittency == 0:
            return SetList

        SetList=deepcopy(SetList)

        for i,set in enumerate(SetList):
            # initially update each frame as seen 0 ago (now)
            seen_frame_ago={s:0 for s in set}
            for j in range(1,intermittency+2):
                for s in seen_frame_ago.keys():
                    # no more frames:
                    if i + j >= len(SetList):
                        continue

                    # if the element is absent now
                    if s not in SetList[i + j]:
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
                        SetList[i + j - k ].add(s)

                    seen_frame_ago[s] = 0
        
        return SetList

    def AutoCorrelation(self,TellTime=None):
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
        list_of_sets : list
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

        TauTimeSeries=list(range(1,self.MaxTau+1))
        TimeSeriesData=[ [] for _ in range(self.MaxTau)]

        self.SetList=self.CorrectIntermittency(self.SetList,intermittency=self.Intermittency)

        # calculate autocorrelation
        for t in range(0,len(self.SetList),self.WindowStep):
            if (TellTime!= None) and t % TellTime == 0:
                print("Already calculate {} frames...".format(t))

            Nt=len(self.SetList[t])

            if Nt == 0:
                continue
            for tau in TauTimeSeries:
                if tau + t >= len(self.SetList):
                    break

                Ntau = len(set.intersection(*self.SetList[t:t + tau + 1]))
                TimeSeriesData[tau-1].append(Ntau/float(Nt))

        TimeSeries = [np.mean(x) for x in TimeSeriesData]

        TauTimeSeries.insert(0,0)
        TimeSeries.insert(0,1)

        return TauTimeSeries,TimeSeries,TimeSeriesData


