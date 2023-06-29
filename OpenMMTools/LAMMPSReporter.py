"""
LAMMPSReporter.py: LAMMPS DUMP reporter for OpenMM, currently only support part of the functions
Authors: Sijia Chen
Last Update: 06/28/2023 (v0.0.2)
What's new (v0.0.2):
- add two new positional argument to handle type name and element names for atoms
- add `try` and `except` to ensure compatibility with old version OpenMM
"""
from __future__ import absolute_import
from __future__ import print_function

__author__ = "Sijia Chen"
__version__ = "0.0.2"

try:
    import openmm.unit as unit
    from openmm import *
    from openmm.app import *
# updated in v1.1: now work with old version OpenMM
except ModuleNotFoundError as mNFErr: 
    import simtk.unit as unit
    from simtk.openmm import *
    from simtk.openmm.app import *
import math
import numpy as np



class LAMMPSReporter(object):
    '''
    LAMMPSReporter outputs information about a simulation, such as atom id and positions, to a file.
    This reporter aims to create LAMMPS-like trajectories. 
    ----------
    file : string or file
        The file to write to, specified as a file name or file object
    reportInterval : int
        The interval (in time steps) at which to write frames
    id : bool = True
        whether to write the ids of atoms
    mass : bool = False
        whether to write the masses of atoms
    mol : bool = False
        whether to write the molecular ids of atoms
    type : bool = False
        whether to write the atom types of atoms
    charge : bool = False
        whether to write the atomic charges of atoms
    element : bool = False
        whether to write the element names of atoms
    elementNames: dict = None
        give the element name of each type, complusory if `element = True`
    wrappedPositions : bool = False
        whether to write the wrapped positions of atoms
    unwrappedPosition : bool = False
        whether to write the unwrapped positions of atoms
    boxImages : bool = False
        whether to write the box images of atoms
    velocity : bool = False
        whether to write the velocities of atoms
    force : bool = False
        whether to write the forces of atoms
    '''

    def __init__(self, file: str, reportInterval: int,  mass=False, id = True, mol = False, type = False, typeNames=None, element=False, elementNames=None, wrappedPosition = False, unwrappedPosition = False, boxImage = False, velocity= False, charge=False, force= False):
        """
        
        Arguments
        ---------
        - file: a string of output file name
        - reportInterval: an integer, output frequency (in timestep)
        - mass: bool, default = False, whether report mass
        - id: bool, default = True, whether report atom index (starting from 1)
        - mol: bool, default = False, whether report atom mol index (starting from 1)
        - type: bool, default = False, whether report user-defined atom type. If set to `True`, then user must provide `typeNames` to provide the type name of each "residue"-"atomname"
        - typeNames, dict, default = None, a dictionary contains information of atom type for each atom name in each residue.
        - element, bool, default = False, whether report element of each atom. If set to `True`, then user must provide `elementNames` to provide the element of each atom name in each residue
        - wrappedPosition: bool, default = False, whether report wrapped positions (in Angstrom) of each atoms, starting from 0 to box boundary
        - unwrappedPosition: bool, default = False, whether report unwrapped positions (in Angstrom) of each atoms
        - boxImage: bool, default = False, whether report box images of each atoms, the box stands in (0, box boundary)
        - velocity: bool, default = False, whether report velocities (in Angstrom/fs) of each atoms.
        - charge: bool, default = False, whether report the charges (in unit electron charge) on atoms
        - force: bool, default = False, whether report forces (in kcal/mol/Angstrom) acting on each atoms.
        
        
        Returns
        -------
        None
        
        """
        
        
        self._reportInterval = reportInterval
        self._openedFile = isinstance(file, str)
        if self._openedFile:
            self._out = open(file, 'w')
        else:
            self._out = file
        self._mass = mass
        self._id = id
        self._mol = mol
        self._type = type
        if self._type:
            if typeNames == None:
                raise ValueError("You must provide the names of each atom in a dictionay {'ResidueName-AtomName':'AtomTypeName'}")
            else:
                self._typeNames = typeNames
        self._element = element
        if self._element :
            if not isinstance(elementNames, dict):
                raise ValueError("You must provide the element name of each atom as a dictionary {'ResidueName-AtomName':'AtomElementName'}")
            else:
                self._elementNames=elementNames
        self._wrappedPosition = wrappedPosition
        self._unwrappedPosition = unwrappedPosition
        self._boxImage = boxImage
        self._velocity = velocity
        self._charge = charge
        self._force = force

        self._boxSizeList = [[], [], []]

        # check if position/velocity/force is needed at each report
        if self._wrappedPosition or self._unwrappedPosition or self._boxImage:
            self._needPositions = True
            if self._unwrappedPosition or self._boxImage:
                self._applyPBC = False
        else:
            self._needPositions = False
            self._applyPBC = None # let the openmm decide if use pbc or not
        if self._velocity:
            self._needVelocities = True
        else:
            self._needVelocities = False
        if self._force:
            self._needForces = True
        else:
            self._needForces = False

        self._needEnergies = False
        self._hasInitialized = False
        


    def describeNextReport(self, simulation):
        """Get information about the next report this object will generate.
        Parameters
        ----------
        simulation : Simulation
            The Simulation to generate a report for
        Returns
        -------
        tuple
            A five element tuple. The first element is the number of steps
            until the next report. The remaining elements specify whether
            that report will require positions, velocities, forces,
            energies, whether apply pbc on positions respectively.
        """
        steps = self._reportInterval - simulation.currentStep % self._reportInterval
        return (steps, self._needPositions, self._needVelocities, self._needForces, self._needEnergies, self._applyPBC)

    def report(self, simulation, state):
        """Generate a report.
        Parameters
        ----------
        simulation : Simulation
            The Simulation to generate a report for
        state : State
            The current state of the simulation
        """
        if not self._hasInitialized:
            self._initializeConstants(simulation)
            self._hasInitialized = True

        # Check for errors.
        self._checkForErrors(state)

        # construct head of each frame
        header = self._constructHeaders(simulation,state)
        # Query for the values
        values = self._constructReportValues(state)

        # Write the values.
        self._writeFrame(header,values)
        try:
            self._out.flush()
        except AttributeError:
            pass

    def _constructReportValues(self, state):
        """Query the simulation for the current state of our observables of interest.
        Parameters
        ----------
        simulation : Simulation
            The Simulation to generate a report for
        state : State
            The current state of the simulation
        Returns
        -------
        A 2D-list (# of atoms, # of properties) of properties for each atom in each row.
        The element in each row is what is required by the user.
        """
        values = []
        box = state.getPeriodicBoxVectors()
        boxSize = np.array([box[0][0].value_in_unit(unit.angstrom), box[1][1].value_in_unit(unit.angstrom), box[2][2].value_in_unit(unit.angstrom)])
        if self._id:
            values.append(self._atomIds)
        if self._mol:
            values.append(self._atomMols)
        if self._type:
            values.append(self._atomTypes)
        if self._charge:
            values.append(self._atomCharges)
        if self._mass:
            values.append(self._atomMasses)
        if self._element:
            values.append(self._atomElements)

        if self._needPositions:
            if self._applyPBC:
                wrappedPositions = state.getPositions(asNumpy = True)
                wrappedPositions = wrappedPositions.value_in_unit(unit.angstroms)
                values.append(list(wrappedPositions[:,0]))
                values.append(list(wrappedPositions[:,1]))
                values.append(list(wrappedPositions[:,2]))
            else:
                unwrappedPositions = state.getPositions(asNumpy = True)
                unwrappedPositions = unwrappedPositions.value_in_unit(unit.angstroms)
                if self._wrappedPosition:
                    wrappedPositions = np.divmod(unwrappedPositions,boxSize)
                    values.append(list(wrappedPositions[:,0]))
                    values.append(list(wrappedPositions[:,1]))
                    values.append(list(wrappedPositions[:,2]))
                if self._unwrappedPosition:
                    values.append(list(unwrappedPositions[:,0]))
                    values.append(list(unwrappedPositions[:,1]))
                    values.append(list(unwrappedPositions[:,2]))
                if self._boxImage:
                    boxImages = np.floor_divide(unwrappedPositions,boxSize).astype(np.int32)
                    values.append(list(boxImages[:,0]))
                    values.append(list(boxImages[:,1]))
                    values.append(list(boxImages[:,2]))

        if self._velocity:
            velocities = state.getVelocities(asNumpy = True)
            velocities = velocities.value_in_unit(unit.angstroms / unit.femtosecond)
            values.append(list(velocities[:,0]))
            values.append(list(velocities[:,1]))
            values.append(list(velocities[:,2]))
            
        if self._force:
            forces = state.getForces(asNumpy = True)
            forces = forces.value_in_unit(unit.kilocalories_per_mole / unit.angstroms)
            values.append(list(forces[:,0]))
            values.append(list(forces[:,1]))
            values.append(list(forces[:,2]))

        return np.array(values).T.tolist()

    def _initializeConstants(self, simulation):
        """Initialize a set of constants required for the reports
        Parameters
        - simulation (Simulation) The simulation to generate a report for
        - state
        """
        system = simulation.system
        self._numAtoms = system.getNumParticles()
        if self._id:
            self._atomIds = range(1, self._numAtoms+1)
        if self._mass:
            self._atomMasses = [system.getParticleMass(i).value_in_unit(unit.dalton) for i in range(system.getNumParticles())]
        if self._charge:
            nonbonded = [f for f in system.getForces() if isinstance(f, NonbondedForce)][0]
            self._atomCharges = []
            for i in range(self._numAtoms):
                charge,_,_ = nonbonded.getParticleParameters(i)
                self._atomCharges.append(charge.value_in_unit(unit.elementary_charge))
        if self._mol:
            self._atomMols = list()
            for i,res in enumerate(simulation.topology.residues()):
                for _ in res.atoms():
                    self._atomMols.append(i+1)
        if self._type:
            self._atomTypes = list()
            for residue in simulation.topology.residues():
                residueName=residue.name
                for atom in residue.atoms():
                    self._atomTypes.append(self._typeNames[f"{residueName}-{atom.name}"])
        if self._element:
            self._atomElements = list()
            for residue in simulation.topology.residues():
                residueName=residue.name
                for atom in residue.atoms():
                    self._atomElements.append(self._elementNames[f"{residueName}-{atom.name}"])

    def _constructHeaders(self,simulation,state):
        """Construct the headers for the CSV output
        Returns: a list of strings giving the title of each observable being reported on.
        """
        box = state.getPeriodicBoxVectors()
        boxSize = np.array([box[0][0].value_in_unit(unit.angstrom), box[1][1].value_in_unit(unit.angstrom), box[2][2].value_in_unit(unit.angstrom)])
        timeStep = simulation.currentStep
        headers = list()
        headers.append("ITEM: TIMESTEP")
        headers.append(f"{timeStep}")
        headers.append("ITEM: NUMBER OF ATOMS")
        headers.append(f"{self._numAtoms}")
        headers.append("ITEM: BOX BOUNDS pp pp pp")
        headers.append(f"0.000000 {boxSize[0]}")
        headers.append(f"0.000000 {boxSize[1]}")
        headers.append(f"0.000000 {boxSize[2]}")
        line="ITEM: ATOMS"
        if self._id:
            line += " id"
        if self._mol:
            line += " mol"
        if self._type:
            line += " type"
        if self._charge:
            line += " q"
        if self._mass:
            line += " mass"
        if self._element:
            line += " element"
        if self._wrappedPosition:
            line += " x y z"
        if self._unwrappedPosition:
            line += " xu yu zu"
        if self._boxImage:
            line += " ix iy iz"
        if self._velocity:
            line += " vx vy vz"
        if self._force:
            line += " fx fy fz"
        headers.append(line)
        return headers

    def _checkForErrors(self, state):
        """Check for errors in the current state of the simulation
        Parameters
         - simulation (Simulation) The Simulation to generate a report for
         - state (State) The current state of the simulation
        """
        if self._needEnergies:
            energy = (state.getKineticEnergy() + state.getPotentialEnergy()).value_in_unit(
                unit.kilojoules_per_mole)
            if math.isnan(energy):
                raise ValueError('Energy is NaN')
            if math.isinf(energy):
                raise ValueError('Energy is infinite')
        if self._needForces:
            forces = (state.getForces(asNumpy = True)).value_in_unit(unit.kilocalories_per_mole / unit.angstroms)
            if np.isnan(forces).any():
                raise ValueError("Force is NaN")
            if np.isinf(forces).any():
                raise ValueError("Force is infinite")

    def _writeFrame(self,header,values):
        for h in header:
            print(h,file=self._out)
        for value in values:
            print(" ".join(str(v) for v in value),file=self._out)

    def __del__(self):
        if self._openedFile:
            self._out.close()

