# Python 3.6.1

import sys
import numpy as np


class AnalysisBase(object):
    def __init__(self, coords,
                 box_volume=None, box_lengths=None, box_boundaries=None,
                 tell_time=False, err=sys.stderr, out=sys.stdout):
        """
        :param coords: numpy.array, xyz coordinates of ONLY necessary atoms/points.
        :param box_volume: float, volume of simulation box
        :param box_lengths: numpy.array, shape=(3,)
        :param box_boundaries: numpy.array, shape=(3,2)
        """
        self.coords = coords

        self._box_volume = None
        self._box_lengths = None
        self._box_boundaries = None
        if box_volume is not None:
            self.box_volume = box_volume
        if box_lengths is not None:
            self.box_lengths = box_lengths
        if box_boundaries is not None:
            self.box_boundaries = box_boundaries

        self._tell_time = tell_time
        self.err = err
        self.out = out
        self.tell_format = '{ClassName} Processing Frame {Time}'.format(ClassName=self.__class__.__name__,
                                                                        Time='{Time}')

    @property
    def box_volume(self):
        if self._box_volume is None:
            raise AttributeError("{ClassName}.box_volume is not set!".format(ClassName=self.__class__.__name__))
        return self._box_volume

    @box_volume.setter
    def box_volume(self, box_volume):
        self._box_volume = box_volume

    @property
    def box_lengths(self):
        if self._box_lengths is None:
            raise AttributeError("{ClassName}.box_lengths is not set!".format(ClassName=self.__class__.__name__))
        return self._box_lengths

    @box_lengths.setter
    def box_lengths(self, box_lengths):
        self._box_lengths = box_lengths
        self._box_volume = np.prod(self._box_lengths)

    @property
    def box_boundaries(self):
        if self._box_boundaries is None:
            raise AttributeError("{ClassName}.box_boundaries is not set!".format(ClassName=self.__class__.__name__))
        return self._box_boundaries

    @box_boundaries.setter
    def box_boundaries(self, box_boundaries):
        self._box_boundaries = box_boundaries
        self._box_lengths = np.diff(self._box_boundaries, axis=1).ravel()
        self._box_volume = np.prod(self._box_lengths)

    @classmethod
    def pbc_dist(cls, origin, coords, box_lengths):
        diff = np.abs(coords - origin)
        pdiff = box_lengths - diff
        mask = diff < box_lengths / 2
        diff = np.where(mask, diff, pdiff)  # mask'em
        dist = np.linalg.norm(diff, axis=1)
        return dist

