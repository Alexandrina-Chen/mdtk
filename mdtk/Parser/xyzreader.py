# -*- encoding: utf-8 -*-
# python 3.10.12
'''
Filename         : xyzreader.py
Description      : A class to read the xyz file and store the data
Time             : 06/22/2024
Last Modified    : 
Author           : Sijia Chen
Version          : 1.0
Email            : sijiachen@uchicago.edu
'''
import os
import numpy as np
import functools

class XyzFrame(object):
    """
    RaptorFrame stores the information of one frame in RAPTOR output file.

    Attributes
    ----------
    n_atoms : int
        The number of complexes in the simulation.
    index : int
        The index of current frame, starting from 0 for each trajectory file.
    comment : str
        the comment line of this frame, namely the second line of each frame in xyz file.
    elements : np.array of str, the element name of atoms in this frame with shape (n_atoms, )
        Get the charges of atoms in this frame.
    positions : np.array, the positions of atoms in this frame with shape (n_atoms, 3)
        Get the positions of atoms in this frame. 
        If the unwrapped positions are available, return the unwrapped positions; 
        otherwise, return the wrapped positions.

    Notes
    -----

    """
    def __init__(self, n_atoms: int, **kwargs) -> None:
        """
        Initialize the RaptorFrame object, which stores the information of one frame
        
        Parameters
        ----------
        n_atoms : int
        
        """
        self._n_atoms = n_atoms
        self._index = -1
        self._timestep = -1
        self._positions = np.zeros((n_atoms, 3), dtype=np.float64)
        self._comment = ""

    @property
    def n_atoms(self) -> int:
        """
        The number of atoms in the simulation.
        
        Returns
        -------
        - n_atoms : int
        
        """
        return self._n_atoms
        
    @n_atoms.setter
    def n_atoms(self, new_n_atoms):
        if new_n_atoms == self._n_atoms:
            pass
        else:
            raise ValueError("Cannot change the number of atoms after initilization!")

    @property
    def index(self) -> int:
        """
        The index of current frame, starting from 0.
        
        Returns
        -------
        - index : int
        
        """
        return self._index
        
    @index.setter
    def index(self, new_index):
        self._index = new_index

    @property
    def comment(self) -> str:
        """
        The comment line of this frame, namely the second line of each frame in xyz file.
        
        Returns
        -------
        - comment : str
        
        """
        return self._comment
        
    @comment.setter
    def comment(self, new_comment):
        self._comment = new_comment


    @property
    def elements(self):
        """
        Get the elements of atoms in this frame.
        
        Returns
        -------
        elements: np.array, the elements of atoms in this frame with shape (n_atoms, )

        """
        if not self._has_elements:
            raise ValueError("This frame does not have elements!")
        else:
            return self._elements
        
    @elements.setter
    def elements(self, new_elements):
        self._has_elements = True
        self._elements = new_elements

    @property
    def positions(self):
        """
        Get the positions of atoms in this frame. 
        If the unwrapped positions are available, return the unwrapped positions; 
        otherwise, return the wrapped positions.
        
        Returns
        -------
        positions: np.array, the positions of atoms in this frame with shape (n_atoms, 3)

        """
        return self._positions
    
    @positions.setter
    def positions(self, new_positions):
        self._positions = new_positions


class XyzReader(object):
    """
    XyzReader reads in the trajectory file and stores the information of all frames in the trajectory file.
    
    Parameters
    ----------
    filename : str
        Name of the RAPTOR output file
    verbose : bool, optional
        If True, print out more information. Mostly for debugging purpose.

    Attributes
    ----------
    n_frames : int
        The number of frames in the trajectory file.
    n_atoms : int
        The number of atoms in the simulation.
    timestep : int

    """
    # todo: add support for reading in multiple trajectory files; "_n_frames" should be the total number of frames in all trajectory files, therefore it should be a class property
    # maybe "_n_frames", "_offsets", "_cached" should be class properties? a experiment for now
    # a example of class property: https://stackoverflow.com/questions/5189699/how-to-make-a-class-property
    # _n_frames = 0
    # _frame = XyzFrame
    # _offsets = []

    def __init__(self, filename, verbose = False,**kwargs) -> None:

        # initialize filename and check sanity
        self.filename = filename
        self._check_sanity_filename()
        self._verbose = verbose # if True, print out more information. Mostly for debugging purpose

        # initialize the number of frames and number of atoms
        self._n_frames = 0
        self._frame = XyzFrame
        self._offsets = []
        self._n_atoms = 0

        # initialize frame class-related arguments
        self._frame_kwargs = self._parse_frame_kwargs(kwargs)
        # initialize cached dictionary
        self._cached = dict()

        self._reopen()
        self._has_initialized = False
        self._read_first_time()

    def _reopen(self):
        self.close()
        self._file = open(self.filename, 'r')
        self.frame= self._frame(self.n_atoms, **self._frame_kwargs)
        self.frame.index = -1
        
    def _check_sanity_filename(self):
        # check if the input file exists
        if isinstance(self.filename, str) and os.path.isfile(self.filename):
            pass
        elif not isinstance(self.filename, str):
            raise TypeError("Input file name must be a string!")
        elif not os.path.isfile(self.filename):
            raise FileNotFoundError("File {} does not exist!".format(self.filename))
        else:
            raise TypeError("Unknown error when parsing the argument {}!".format(self.filename))

        # check if the input file is empty
        if os.stat(self.filename).st_size == 0:            
            raise ValueError("Input file {} is empty!".format(self.filename))
        
    @staticmethod
    def _parse_frame_kwargs(kwargs) -> dict:
        frame_kwargs = dict()
        return frame_kwargs


    # @property
    # def time(self):
    #     return self.frame.time

    # @property
    # def frame(self) -> XyzFrame:
    #     return self._frame
    
    # @frame.setter
    # def frame(self, new_frame):
    #     self._frame = new_frame
    
    
    @property
    def index(self) -> int:
        """
        return the index of the current frame, 0-based

        Returns
        -------
        - index : int
        """
        return self.frame.index
    
    @functools.cached_property
    def n_atoms(self):
        with open(self.filename, 'r') as f:
            self._n_atoms = int(f.readline())
        self._cached["n_atoms"] = self._n_atoms
        return self._n_atoms

    @functools.cached_property
    def n_frames(self):
        lines_per_frame = 2 + self.n_atoms
        counter = 0
        offsets = []
        with open(self.filename, 'r') as f:
            line = True
            while line:
                if counter % lines_per_frame == 0 :
                    offsets.append(f.tell())
                counter += 1
                line = f.readline()
        self._offsets = offsets[:-1] # last position is EOF
        self._n_frames += len(self._offsets)
        self._cached["n_frames"] = self._n_frames
        return self._n_frames
    
    @property
    def elements(self):
        """
        return the elements of atoms in the current frame

        Returns
        -------
        - elements : np.array, the elements of atoms in this frame with shape (n_atoms, )
        """
        
        return self.frame.elements
    
    @property
    def positions(self):
        """
        return the positions of atoms in the current frame

        Returns
        -------
        - positions : np.array, the positions of atoms in this frame with shape (n_atoms, 3)
        """
        return self.frame.positions
    
    def __len__(self):
        return self.n_frames
    
    def next(self) -> XyzFrame:
        try:
            frame = self._read_next_frame()
        except EOFError:
            self.rewind()
            raise StopIteration from None
        else:
            pass
        
        return frame


    def __next__(self):
        return self.next()

    def close(self):
        if hasattr(self, '_file'):
            self._file.close()

    def _read_next_frame(self) -> XyzFrame:
        # initialize 
        if not self._has_initialized:
            raise RuntimeError("Please initialize the reader first!")
        
        f = self._file
        frame = self.frame
        frame.index += 1
        if frame.index >= len(self):
            raise EOFError("End of file reached!")
        

        # go to the TIMESTEP line of this frame
        f.seek(self._offsets[frame.index])

        # read number of atoms
        frame.n_atoms = int(f.readline().strip())
        
        # read the comment
        frame.comment = f.readline().strip()
        
        # read the remaing lines of this frame
        atom_info=[]
        for _ in range(self.n_atoms):
            atom_info.append(f.readline().strip().split())
        atom_info = np.array(atom_info)

        # parse the atom elements and positions
        frame.elements = atom_info[:,0]
        frame.positions = atom_info[:,1:4].astype(np.float64)

        return frame
        
        
    def rewind(self) -> XyzFrame:
        self._reopen()
        self.next()

    def _read_frame(self, index) -> XyzFrame:
        """
        Read a specific frame `index` from the trajectory and return it.
        
        Parameters
        ----------
        index : int
            Index of the frame to read. 0-based.

        Returns
        -------
        frame : RaptorFrame
            Frame object.
        """
        if index < 0 or index >= len(self):
            raise IndexError("Index out of range!")
        self.frame.index = index-1
        # self._file.seek(self._offsets[index])
        return self._read_next_frame()
    
    def __getitem__(self, index) -> XyzFrame:
        """Return the frame with index `index`.

        Parameters
        ----------
        index : int
            Index of the frame to read. 0-based.

        """
        if isinstance(index, int) and index >= 0:
            return self._read_frame(index)
        elif isinstance(index, slice):
            start = index.start or 0
            stop = index.stop or len(self)
            step = index.step or 1
            return [self._read_frame(i) for i in range(start, stop, step)]
        elif isinstance(index, (list, np.ndarray)):
            return [self._read_frame(i) for i in index]
        else:
            raise TypeError("Invalid index {} with type {}!".format(index,type(index)))

    def __del__(self):
        self.close()

    def __iter__(self):
        """ Iterate over trajectory frames. """
        self._reopen()
        return self

    def _read_first_time(self) -> None:
        """
        Check the information of the first frame and initialize the frame class
        initialized the frame class with the number of atoms and number of frames
        """
        _ = self.n_atoms
        _ = self.n_frames
        self._has_initialized = True

    def __repr__(self):
        return f"<XyzReader: {self.filename}>"

    def write_to_xyz(self, filename, start=0, stop=None, step=1):
        """
        write the trajectory to a new file according to the slice

        Parameters
        ----------
        filename : str
            The name of the new file
        start : int, optional
            The start index of the slice, default is 0
        stop : int, optional
            The stop index of the slice, default is None
        step : int, optional
            The step of the slice, default is 1
        """
        # check the sanity of the input filename
        if not isinstance(filename, str):
            raise TypeError("Input file name must be a string!")
        if os.path.isfile(filename):
            # warn the user that the file already exists
            print("Warning: file {} already exists and will be overwritten if you proceed!".format(filename))
            # let the user decide whether to overwrite the file
            user_input = input("Do you want to overwrite the file? (y/n)")
            if user_input.lower() not in ["y", "yes"]:
                print("User chose not to overwrite the file. Exiting...")
                return
            else:
                print("User chose to overwrite the file. Proceeding...")
        else:
            pass
        
        # check the sanity of the slice
        start, stop, step = self._slice_index_sanity_check(start, stop, step)
        print(f"Writing the trajectory to the new file {filename} with slice {start}:{stop}:{step}...")
        
        # write the trajectory to the new file
        with open(filename, 'w') as f:
            for i in range(start, stop, step):
                frame = self._read_frame(i)
                f.write(f"{frame.n_atoms}\n")
                f.write(f"{frame.comment}\n")
                for element, position in zip(frame.elements, frame.positions):
                    f.write(f"{element:<6s} {position[0]:18.10f} {position[1]:18.10f} {position[2]:18.10f}\n")

    def _slice_index_sanity_check(self, start, stop, step):
        start = self._index_sanity_check(start)
        if start >= len(self):
            raise ValueError(f"Start index is too large, start={start} >= {len(self)}")
        if stop is None:
            stop = len(self)
        stop = self._index_sanity_check(stop)
        
        if step <= 0:
            raise ValueError("Step must be positive integer!")
        elif step > stop - start:
            raise ValueError("Step is too large for the given slice!")
        return start, stop, step

    def _index_sanity_check(self, index):
        if index < -len(self):
            raise IndexError(f"Index out of range, {index} < -{len(self)}")
        elif index < 0:
            index = len(self) + index
        elif index > len(self):
            raise IndexError(f"Index out of range, {index} > {len(self)}")
        return index


class XyzWriter(object):
    """
    XyzWriter writes the trajectory file from the information of all frames in the trajectory file.
    
    Parameters
    ----------
    filename : str
        Name of the xyz output file
    verbose : bool, optional
        If True, print out more information. Mostly for debugging purpose.

    Attributes
    ----------
    n_frames : int
        The number of frames already written to the trajectory file.
    n_atoms : int
        The number of atoms in the simulation.

    """
    def __init__(self, filename, verbose = False) -> None:
        self.filename = filename
        self._verbose = verbose

        # check if the output file exists and whether to overwrite it, if it exists; if not or user chooses to overwrite, initialize the file object with write mode
        self._check_sanity_filename()

        # initialize the number of frames and number of atoms
        self._n_frames = 0
        self._n_atoms = 0


    def _check_sanity_filename(self):
        # check if the output file exists
        if os.path.isfile(self.filename):
            # warn the user that the file already exists
            print("Warning: file {} already exists and will be overwritten if you proceed!".format(self.filename))
            # let the user decide whether to overwrite the file
            user_input = input("Do you want to overwrite the file? (y/n)")
            if user_input.lower() not in ["y", "yes"]:
                print("User chose not to overwrite the file. Exiting...")
                raise FileExistsError("File {} already exists!".format(self.filename))
            else:
                print("User chose to overwrite the file. Proceeding...")
        else:
            pass
        self._file = open(self.filename, 'w')

    def write(self, n_atoms, elements, positions, comment=None):
        """
        Write the frame to the trajectory file

        Parameters
        ----------
        n_atoms : int
            The number of atoms in the frame
        elements : list of str
            The element names of the atoms in the frame
        positions : np.array
            The coordinates of the atoms in the frame with shape (n_atoms, 3)
        comment : str, optional
            The comment line of the frame, default is None
        """

        # check the sanity of the input
        if n_atoms != len(positions) or n_atoms != len(elements):
            raise ValueError("The number of atoms, elements, and positions do not match!")
        if self._n_frames == 0:
            self._n_atoms = n_atoms
        else:
            if n_atoms != self._n_atoms:
                raise ValueError("The number of atoms in this frame does not match the previous frames!")
        if comment is None:
            comment = f"Frame {self._n_frames}"

        # write the frame to the trajectory file
        self._write_frame(n_atoms, comment, elements, positions)

        # if the frame is written successfully, update the number of frames
        self._n_frames += 1

    def _write_frame(self, n_atoms, comment, elements, positions):
        self._file.write(f"{n_atoms}\n")
        self._file.write(f"{comment.strip()}\n")
        for element, position in zip(elements, positions):
            self._file.write(f"{element:<6s} {position[0]:18.6f} {position[1]:18.6f} {position[2]:18.6f}\n")

    def close(self):
        if hasattr(self, '_file'):
            self._file.close()


def test(filename):
    reader = XyzReader(filename,verbose=True)
    pass

if __name__ == "__main__":
    file=""
    test(file)