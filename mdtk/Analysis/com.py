# -*- encoding: utf-8 -*-
# python 3.10.12
'''
Filename         : com.py
Description      : calculate the com of a group of atoms given their coordinates, group, mass and box information
Time             : 2024/02/05
Author           : Sijia Chen
Version          : 1.0
Email            : sijiachen@uchicago.edu
'''

#-----------Import Modules-----------#
import numpy as np
from .base import AnalysisBase

class ComAnalysis(AnalysisBase):
    def __init__(self, coords, group, mass, unwrap = False, **kwargs):
        """
        Description
        -----------
        Calculate the com of a group of atoms given their coordinates, group index, mass and box information
        This class is only used to calculate the com of the same kind of molecules,

        Parameters
        ----------
        coords : numpy.array, shape=(n, 3)
            xyz coordinates of atoms
        group : numpy.array, shape=(# of molecules, # of atoms in each molecule)
            the group of atoms that form a molecule, each row is a molecule, which contains the indices of atoms (0-based) according to the input coords
        mass : numpy.array, shape=(# of atoms in each molecule, )
            the mass of each atom in the molecule
        unwrap : bool, default=False
            whether to unwrap the coordinates of atoms. 
            If they are wrapped in the box, set it to `True`. 
            If `True`, at least one of `box_lengths` and `box_boundaries` should be given

        Properties
        ----------
        com : numpy.array, shape=(# of molecules, 3)
            the center of mass of the group of atoms
        coords : numpy.array, shape=(n, 3)
            xyz coordinates of atoms
        group : numpy.array, shape=(# of molecules, # of atoms in each molecule)
            the group of atoms that form a molecule, each row is a molecule, which contains the indices of atoms (0-based) according to the input coords
        mass : numpy.array, shape=(# of atoms in each molecule, )
            the mass of each atom in the molecule
        unwrap : bool, default=False
            whether to unwrap the coordinates of atoms. 
            If they are wrapped in the box, set this to `True`. 
            If `True`, at least one of `box_lengths` and `box_boundaries` should be given
        box_lengths : None or numpy.array, shape=(3,)
            the lengths of the box
        box_boundaries : None or numpy.array, shape=(3,2)
            the boundaries of the box

        """
        super().__init__(coords, **kwargs)
        self.group = group
        self.mass = mass
        self.unwrap = unwrap
        if self.unwrap:
            if self._box_lengths is None and self._box_boundaries is None:
                raise ValueError('At least one of box_lengths and box_boundaries should be given if unwrap is True')
        # have not calculated the com yet
        self._com = None

    def _calc_com(self):
        """
        Description
        -----------
        Calculate the center of mass of the group of atoms

        Returns
        -------
        com : numpy.array, shape=(# of molecules, 3)
            the center of mass of the group of atoms
        """
        if self.unwrap:
            coords = self._unwrap_coords(self.coords[self.group], self.box_lengths, self.box_boundaries)
        else:
            coords = self.coords[self.group]
        self._com = np.average(coords, axis=1, weights=self.mass)
        return self._com

    def _unwrap_coords(self, coords, box_lengths, box_boundaries):
        """
        Description
        -----------
        Unwrap the coordinates of atoms in the box

        Parameters
        ----------
        coords : numpy.array, shape=(# of molecules, # of atoms in each molecule, 3)
            xyz coordinates of atoms
        box_lengths : None or numpy.array, shape=(3,)
            the lengths of the box
        box_boundaries : None or numpy.array, shape=(3,2)
            the boundaries of the box
        
        Returns
        -------
        unwrapped_coords : numpy.array, shape=(n, 3)
            unwrapped coordinates of atoms
        
        Notes
        -----
        At least one of `box_lengths` and `box_boundaries` should be given
        """
        if box_lengths is None:
            box_lengths = np.diff(box_boundaries, axis=1).ravel()
        for mol_coords in coords:
            for i, coord in enumerate(mol_coords[:-1]):
                diff = mol_coords[i+1:] - coord
                mask = np.abs(diff) > (box_lengths / 2)
                diff[mask] -= np.sign(diff[mask]) * box_lengths
                mol_coords[i+1:] = coord + diff
        return coords
    
    @property
    def group(self):
        return self._group

    @group.setter
    def group(self, group):
        self._group = group

    @property
    def mass(self):
        return self._mass

    @mass.setter
    def mass(self, mass):
        self._mass = mass

    @property
    def unwrap(self):
        return self._unwrap

    @unwrap.setter
    def unwrap(self, unwrap):
        self._unwrap = unwrap

    @property
    def com(self):
        if self._com is None:
            self._calc_com()
        return self._com
    
    @property
    def coords(self):
        return self._coords

    @coords.setter
    def coords(self, coords):
        self._coords = coords
        # after the coordinates are changed, the com should be recalculated
        self._com = None
            
    
def calc_com_from_traj(traj_reader, group, mass, unwrap = False, skip_beginning = 0, tell_time = 0, **kwargs):
    """
    Description
    -----------
    Calculate the center of mass of a group of atoms given their coordinates, group, mass and box information from a trajectory

    Parameters
    ----------
    traj_reader : TrajectoryReader (LammpsReader, XyzReader, etc.)
        the **Opened** trajectory reader
    group : numpy.array, shape=(# of molecules, # of atoms in each molecule)
        the group of atoms that form several molecules, 
        each row is a molecule, which contains the indices of atoms (0-based) according to the input coords. 
    mass : numpy.array, shape=(# of atoms in each molecule, )
        the mass of each atom in the molecule
    unwrap : bool, default=False
        whether to unwrap the coordinates of atoms. 
        If they are wrapped in the box, set it to `True`. 
        If `True`, at least one of `box_lengths` and `box_boundaries` should be given
    skip_beginning : int, default=0
        the number of frames to skip at the beginning
    tell_time : int, default=0
        frequency to print the time, if 0 or less, do not print the time; otherwise, print the time every `tell_time` frames
    kwargs : dict
        other parameters for the TrajectoryReader

    Returns
    -------
    com : numpy.array, shape=(# of frames, # of molecules, 3)
        the center of mass of the group of atoms for each frame
    timesteps : numpy.array, shape=(# of frames, )
        the timestep of each frame
    """
    com = []
    timesteps = []
    com_analysis = ComAnalysis(None, group, mass, unwrap, **kwargs)
    for i, frame in enumerate(traj_reader):
        if i < skip_beginning:
            continue
        if tell_time > 0 and i % tell_time == 0:
            print(f"Have calculated {i} frames...")
        com_analysis.coords = frame.unwrapped_positions
        com.append(com_analysis.com)
        timesteps.append(frame.timestep)
    return np.array(com), np.array(timesteps)
        

