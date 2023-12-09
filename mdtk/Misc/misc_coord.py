
"""
Some miscellaneous functions for coordinate conversions.
"""

import numpy as np

def wrap_coordinates(coordinates, boundary):
    """
    Parameters
    ------------
    coordinates: ndarray, shape = (n1, n2, ..., # of dimensions)
    boundary: ndarray, shape = (# of dimensions, 2)
    
    Returns
    ------------
    wrapped_coordinates: ndarray, shape = (n1, n2, ..., # of dimensions)
    """
    dimension = coordinates.shape[-1]
    if dimension != boundary.shape[0]:
        raise ValueError("Dimension of coordinates and boundary do not match")
    moved_coordinates = coordinates - boundary[:, 0]
    wrapped_coordinates = moved_coordinates % np.abs(boundary[:, 1] - boundary[:, 0])
    wrapped_coordinates += boundary[:, 0]
    return wrapped_coordinates

def read_mols_file(filename, num_ion_pairs, max_layer):
    """
    Parameters
    ------------
    filename: str
    num_ion_pairs: int, the number of ion pairs
    max_layer: int
    
    Returns
    ------------
    cation_mols_set_list: a list of set. Contains set of cation mol ids in each frame, shape = (maxlayer*2, # of frames), upper layer is [0:maxlayer, : ], lower layer is [maxlayer:, :]
    anion_mols_set_list: a list of set. Contains set of anion mol ids in each frame, shape = (maxlayer*2, # of frames), upper layer is [0:maxlayer, : ], lower layer is [maxlayer:, :]
    
    Note
    ------------
    The mol id begins with 1, instead of 0;
    After the reading, both cation and anion mol id range from 1 (included) to numionpairs (included);
    this function is designed for reading mols file from ITIM analysis;
    this function is used for pure ionic liquid systems with only one type of cation and one type of anion;
    this function assumes that the mol id of cation is smaller than the mol id of anion;

    """
    # prepare a list to store mol ids
    print("Now reading Mols file {}".format(filename.split('/')[-1]))
    cation_mols_set_list=[]
    anion_mols_set_list=[]
    for ilayer in range(max_layer):
        cation_mols_set_list.append([])
        cation_mols_set_list.append([])
        anion_mols_set_list.append([])
        anion_mols_set_list.append([])
    # read mols file
    with open(filename,'r') as file:
        for line in file:
            if line.startswith("#"):
                continue
            sline=list(map(lambda x: int(x), line.split()))
            CationSet=set()
            AnionSet=set()
            for molid in sline[2:]:
                if molid <= num_ion_pairs:
                    CationSet.add(molid)
                else:
                    AnionSet.add(molid - num_ion_pairs)
            cation_mols_set_list[sline[0]*2+sline[1]].append(CationSet)
            anion_mols_set_list[sline[0]*2+sline[1]].append(AnionSet)
    return cation_mols_set_list,anion_mols_set_list

def read_com_file(filename, num_com = -1, num_dimensions = -1):
    """
    Parameters
    ------------
    filename: str
    num_com: int, the number of center of mass in one frame; default = -1, which means the number of com is determined by the first frame
    num_dimensions: int, the number of dimensions; default = -1, which means the number of dimensions is determined by the first frame

    Returns
    ------------
    frame_timestep_list: a list of int, shape = (# of frames)
    com_list: ndarray, shape = (# of frames, num_com, # of dimensions)
    """
    # prepare a list to store com
    print("Now reading COM file {}".format(filename.split('/')[-1]))
    frame_timestep_list=[]
    com_list=[]
    if num_com == -1:
        with open(filename,'r') as file:
            counter = 0
            trigger = False
            for line in file:
                if line.startswith("#"):
                    if trigger:
                        break
                    else:
                        trigger = True
                        continue
                counter += 1
        num_com = counter
    print(f"There are supposed to have {num_com} in one frame")

    if num_dimensions == -1:
        with open(filename,'r') as file:
            for line in file:
                if line.startswith("#"):
                    continue
                sline=line.strip().split()
                num_dimensions = len(sline) - 1
                break
    # read com file
    with open(filename,'r') as file:
        for line in file:
            if line.startswith("#"):
                frame_timestep_list.append(int(line.strip().split()[-1]))
                continue
            sline=list(map(lambda x: float(x), line.split()))
            com_list.append(sline[1:])
    com_list=np.array(com_list).reshape((-1,num_com,num_dimensions))
    if com_list.shape[0] != len(frame_timestep_list):
        raise ValueError("The number of frames in COM file does not match the number of timesteps in COM file")
    return frame_timestep_list, com_list
