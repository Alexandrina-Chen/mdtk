# Python 3.6.1

import numpy as np

from .base import AnalysisBase


class DensityAnalysis(AnalysisBase):
    def __init__(self, coords, **kwargs):
        super().__init__(coords, **kwargs)

    def density(self, start, stop, bins, style='histogram', bandwidth=2.4):
        """
        :param start: numpy.array; One corner of the box
        :param stop: numpy.array; The other corner of the box
        :param bins: numpy.array; Number of the slices on x,y,z directions; use numpy.nan to indicate low-dimension bins
        :param style: str; Way of calculating the density. Supported: histogram, gaussian
        :param bandwidth: float; Bandwidth in Gaussian kernel, only valid with style "gaussian"
        :return: Histogram style: histogram densities and bin edges of the histogram;
                 Gaussian style: Gaussian densities and coords of the densities.
        """
        if np.sum(np.isnan(bins)) == 3:
            raise ValueError('At lease one dimension of bins should be given.')

        coords = self.coords[:, ~np.isnan(bins)]
        start_nec = start[~np.isnan(bins)]
        stop_nec = stop[~np.isnan(bins)]
        bins_nec = bins[~np.isnan(bins)].astype(np.int32)
        if style == 'histogram':
            hist, edges = np.histogramdd(coords, bins=bins_nec, range=np.vstack((start_nec, stop_nec)).T)
            bin_volume = np.prod(
                np.append(np.abs(stop - start)[np.isnan(bins)], (np.abs(stop - start)[~np.isnan(bins)] / bins[~np.isnan(bins)])))
            hist /= bin_volume
            return hist, edges
        elif style == 'gaussian':
            box_lengths = self.box_lengths[~np.isnan(bins)]
            edges = np.vstack((start, stop)).T
            xyz = [np.linspace(edge[0], edge[1], num_bin) for edge, num_bin in zip(edges, bins_nec)]
            origins = np.vstack(list(map(np.ravel, np.meshgrid(*xyz)))).T
            densities = np.array(list(map(lambda x: self.gaussian_kde_pbc(x, coords, box_lengths, bandwidth), origins)))
            return densities, origins
        else:
            raise ValueError('Unsupported style.')

    @classmethod
    def gaussian_kde_pbc(cls, origin, coords, box_lengths, bandwidth, dimension=3):
        diff = np.abs(coords - origin)  # if no pbc_dist
        pdiff = box_lengths - diff  # if all pbc_dist
        mask = diff < box_lengths / 2  # find which element need pbc_dist
        diff = np.where(mask, diff, pdiff)  # mask'em
        tdiff = diff / bandwidth
        if len(tdiff.shape) > 1:
            energy = np.sum(tdiff * tdiff, axis=1) / 2
        else:
            energy = tdiff * tdiff / 2
        return np.sum(np.exp(-energy)) * ((2 * np.pi * bandwidth)** -dimension/2)
