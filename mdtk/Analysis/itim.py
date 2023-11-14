# Python 3.9.6
# Last Update: 07/12/2022
# itim analysis

from __future__ import print_function
from cmath import sqrt
import platform
from multiprocessing import Process, Queue
import numpy as np
from scipy.spatial import cKDTree

#from . import messages
#from . import utilities
#from .surface import SurfaceFlatInterface as Surface



class ItimAnalysis():

    @property
    def layers(self):
        """
        Access the layers as numpy arrays of atom indices, shape = (2, max_layer, # of atoms in the layer)
        
        """
        return self._layers_origin_indices

    @property
    def surfaces(self):
        """
        Access the surfaces positions as numpy arrays of z positions, shape = (2, max_layers, # of meshpoints)
        Note that the position of the surface is offseted by the user-defined value
        
        """
        return self._surfaces

    @property
    def surfmols(self):
        """
        Access the surfaces molecules' ids as numpy arrays, shape = (2, max_layers, # of molecules in the layer)
        
        """
        return self._layers_mol_indices
        
    def __init__(self, coords, types, dimensions, boundaries = None, mols=None, group = None, alpha = 1.5, molecular = False,max_layers=1,radii_dict=None, cluster_cut = None,side=None,mesh=1.0,autoassign=True, multiproc = False,offset=0.0, **kwargs):
        """
        The interface MUST be normal to z axis. The z coodinates of upper layers MUST be positive, the z coodinates of lower layers MUST be negative. Use the argument squ`offset` to offset the z coordinates to achieve the requirement

        Parameters
        -----
        - coords: np.ndarray, coordinates of all atoms, shape = (# of atoms,3) (required)
        - types: np.ndarray, atom types of all atoms, shape = (# of atoms) (required)
        - dimensions: list/np.ndarray of float, dimension in each direction, shape = (3, ) (required)
        - boundaries: list/np.ndarray of float, lower and upper bondaries in each direction, shape = (3, 2) (default = None, if default, the boundaries will be set to [[0,x],[0,y],[0,z]])
        - mols: np.ndarray/None, atom mols of all atoms, shape = (# of atoms) (required if molecular = True) (default = None)
        - group: list/np.ndarray/None, desired analysis atom types, shape = (# of analysis atom types); None means all (defalut = None)
        - alpha: float, radii of probe sphere (default = 1.5)
        - molecular: bool, Switches between search of interfacial molecules / atoms (default = False)
        - max_layers: int, the number of layers to be identified (default = 1)
        - radii_dict: dict, dictionary with the atomic radii of the elements in the group. (required) # TODO: finish guess from mass 
        - mesh: float, the distance between two testing lines in one dimension (default = 1.0)
        - autoassign: bool,
        - multiproc: bool, if multiprocessing, (default = False)
        - offset: float, the value to offset the z coordinates (default = 0.0). `offset = 1.0` means all atoms' z coordinates will add 1.0

        """
        self.autoassign=autoassign
        self.system=platform.system()
        self.multiproc = multiproc
        self.all_atoms_coords = np.array(coords)
        self.all_atoms_types = np.array(types)
        # create box and boundary
        self.dimensions = np.array(dimensions)
        if boundaries is None:
            self.boundaries = np.array([[0,self.dimensions[0]],[0,self.dimensions[1]],[0,self.dimensions[2]]])
        else:
            self.boundaries = np.array(boundaries)
        self.offset=offset
        # analysis groups
        self.analysis_types = group
        # assign analysis groups
        if self.analysis_types is None:
            self.analysis_types = np.unique(self.all_atoms_types)
            self.analysis_atoms_coords = self.all_atoms_coords
            self.analysis_atoms_number = len(self.analysis_atoms_coords)
            self.analysis_atoms_ids = np.arange(self.analysis_atoms_number)
            self.analysis_atoms_types = self.all_atoms_types
        elif isinstance(self.analysis_types,(list,np.ndarray)):
            analysis_ornot=np.sum(list(map(lambda x: self.all_atoms_types == x, self.analysis_types)),axis=0,dtype=bool) # find atoms that belong to desired atom types
            self.analysis_atoms_ids = np.flatnonzero(analysis_ornot)
            self.analysis_atoms_coords = self.all_atoms_coords[analysis_ornot]
            self.analysis_atoms_types = self.all_atoms_types[analysis_ornot]
            self.analysis_atoms_number = len(self.analysis_atoms_ids)
        else:
            raise RuntimeError("Wrong analysis group")
        if (len(self.analysis_atoms_coords) == 0):
            raise RuntimeError("Undefined analysis group")
        # switch molecular and check mols existency
        self.molecular=molecular
        self.all_atoms_mols = mols
        if self.molecular==True:
            try:
                self.analysis_atoms_mols = self.all_atoms_mols[self.analysis_atoms_ids]
            except:
                raise TypeError("Undefined molcular numbers")
        # assign radii
        if radii_dict == None:
            """
            TODO: guess the atom radii based on atom type
            """
            raise TypeError("Undefined atom radii")
        else:
            self.radii_dict = radii_dict
        self.analysis_atoms_radii = self._assign_radii(self.radii_dict,self.analysis_atoms_types,self.analysis_types)
        # alpha
        if (alpha > 0) and (alpha < np.amin(self.dimensions)/2) :
            self.alpha = alpha
        else:
            raise ValueError("alpha (={}) is too large or too small with box dimensions ({})".format(alpha,self.dimensions))
        # mesh
        self.target_mesh = mesh
        self.grid = None
        # initialize layers
        self.max_layers=max_layers
        self._layers_analysis_indices = [[],[]] # indices in the analysis group
        self._layers_origin_indices = [[],[]] # indices in all atom group
        self._layers_mol_indices = [[],[]] # indices of molecule ids
        #self._surfaces = np.zeros([2,self.max_layers])

        # calculate ITIM
        self._assign_layers()

        

    @staticmethod
    def _assign_radii(radii_dict,atoms_types,analysis_types):
        # check if all type radii are present
        for type in analysis_types:
            if not type in radii_dict.keys():
                raise ValueError("The radius of type {} doesn't exist!".format(type))

        # to efficiently replace type with radii, check `https://stackoverflow.com/questions/55949809/efficiently-replace-elements-in-array-based-on-dictionary-numpy-python`
        # In lammps, types are all integers
        # Otherwise, replace with intergers
        radii_types = np.array(list(radii_dict.keys()))
        radii_radii = np.array(list(radii_dict.values()))
        mapping_ar = np.zeros(radii_types.max()+1,dtype=radii_radii.dtype)
        mapping_ar[radii_types] = radii_radii
        return mapping_ar[atoms_types]

    def _create_mesh(self):
        """
        Mesh assignment, create the grid and intialize a cKDTree object to search gridpoints touched by molecules
        """
        dimensions = self.dimensions
        n = np.array([np.ceil(b/self.target_mesh) for b in dimensions],dtype=np.int_)
        d = dimensions/n
        self.mesh_nx = n[0]
        self.mesh_ny = n[1]
        self.mesh_dx = d[0]
        self.mesh_dy = d[1]

        _x = np.linspace(0,dimensions[0],num=self.mesh_nx,endpoint=False)
        _y = np.linspace(0,dimensions[1],num=self.mesh_ny,endpoint=False)
        _X, _Y = np.meshgrid(_x,_y)
        self.meshpoints=np.array([_X.ravel(),_Y.ravel()]).T
        self.meshtree = cKDTree(self.meshpoints,boxsize=dimensions[:2])

    def _wrap_coords(self):
        # move all atom x,y coordinates to starting with 0 and then wrap them into pbc box; offset the z coordinate
        all_atoms_coords_dists = self.all_atoms_coords - self.boundaries[:,0]
        wrapped_all_atoms_coords_dists = all_atoms_coords_dists % self.dimensions
        self.all_atoms_coords = wrapped_all_atoms_coords_dists + np.array([0.0, 0.0, self.boundaries[2,0] + self.offset])
        self.analysis_atoms_coords = self.all_atoms_coords[self.analysis_atoms_ids]

    def _prepare_box(self):
        self.original_all_atoms_coords = self.all_atoms_coords
        self.original_analysis_atoms_coords = self.analysis_atoms_coords
        self._wrap_coords()
    
    def _prepare_layers_assignment(self):
        self._create_mesh()
        size=(2,int(self.max_layers),int(self.mesh_nx)*int(self.mesh_ny))
        self.mask=np.zeros(size,dtype=int)
        self._prepare_box()
        self._surfaces = np.zeros(size)
    
    def _reset_labels(self):
        self.all_atoms_layers = np.zeros_like(self.all_atoms_mols) - 1 # -1 means in bulk
        self.all_atoms_layersside = np.zeros_like(self.all_atoms_mols) - 1 # -1 means not assigned

    def _touch_lines(self,atom,_x,_y,_z,_radii):
        return self.meshtree.query_ball_point([_x[atom],_y[atom]],_radii[atom]+self.alpha)

    def _append_layers(self,uplow,layer,layers):
        inlayer_indices = np.flatnonzero(self._seen[uplow] == layer+1)
        
        if self.molecular is True:
            inlayer_mol_ids=np.unique(self.analysis_atoms_mols[inlayer_indices])
            self._layers_mol_indices[uplow].append(inlayer_mol_ids) # append mol ids of molecules in the layer
            indices=np.sum(list(map(lambda x: self.analysis_atoms_mols == x, inlayer_mol_ids)),axis=0,dtype=bool)
            self._seen[uplow][indices] = layer + 1
        else:
            self._seen[uplow][inlayer_indices] = layer + 1
            indices=inlayer_indices

        if np.sum(indices) == 0:
            raise Exception("Empty layer {}".format(layer+1))

        layers.append(indices)



    def _assign_one_side(self,uplow,sorted_indices,_x,_y,_z,_radii,queue = None):
        layers = []
        uplowsign=[1.0,-1.0]
        for layer in range(0,self.max_layers):
            # test if the test lines have been touched or not
            mask=self.mask[uplow][layer]
            for atom in sorted_indices:
                if self._seen[uplow][atom] != 0:
                    continue

                touched_lines = self._touch_lines(atom,_x,_y,_z,_radii)
                _submask = mask[touched_lines]

                # find surface of the points
                for touched_line in touched_lines:
                    #for each line, find it topest stop point
                    # assume z of top layers > 0 and z of low layers < 0
                    abs_touched_point = (abs(_z[atom])+(sqrt((_radii[atom]+self.alpha)**2 - (_x[atom]-self.meshpoints[touched_line][0])**2 - (_y[atom]-self.meshpoints[touched_line][1])**2)-self.alpha)).real
                    """if touched_line% 2000 == 0:
                        print("touched_line of uplow {}: {}; self.meshpoints[touched_line]: {}; _z[atom]: {}; _radii[atom]: {}; _x[atom]: {}; _y[atom]: {}; abs_touched_point: {}; old surface: {}".format(uplow,touched_line,self.meshpoints[touched_line],_z[atom],_radii[atom],_x[atom],_y[atom],abs_touched_point,self._surfaces[uplow][layer][touched_line]))"""
                    if abs_touched_point > abs(self._surfaces[uplow][layer][touched_line]):
                        self._surfaces[uplow][layer][touched_line] = uplowsign[uplow] * abs_touched_point
                        """if touched_line% 2000 == 0:
                            print("new touched point of {}: {}".format(uplow,self._surfaces[uplow][layer][touched_line]))"""

                if (len(_submask[_submask==0]) == 0):
                    # no new contact, move to next atom
                    continue

                # mark touched lines
                mask[touched_lines] = 1

                # mark sorted atoms
                self._seen[uplow][atom] = 1+layer

                # if all lines have been touched, create a group that includes all atoms in this layer
                if np.sum(mask) == len(mask):
                    self._append_layers(uplow,layer,layers)
                    break

        if (queue is None):
            return layers
        else:
            queue.put(layers)


    def _assign_layers(self):
        up, low = 0,1
        self._reset_labels()
        self._prepare_layers_assignment()

        _radii = self.analysis_atoms_radii
        size = self.analysis_atoms_number
        self._seen = [np.zeros(size, dtype=np.int8), np.zeros(size, dtype=np.int8)]
        _x = self.analysis_atoms_coords[:,0]
        _y = self.analysis_atoms_coords[:,1]
        _z = self.analysis_atoms_coords[:,2]

        # sort z+radius to touch lines from top/bottom to center
        sort = np.argsort(_z + _radii * np.sign(_z))

        if self.multiproc :
            proc,queue = [None, None],[Queue(),Queue()]
            proc[up] = Process(target=self._assign_one_side,args=(up,sort[::-1],_x,_y,_z,_radii,queue[up]))
            proc[low] = Process(target=self._assign_one_side,args=(up,sort[::],_x,_y,_z,_radii,queue[low]))
            for p in proc:
                p.start()
            for uplow in [up,low]:
                for index,group in enumerate(queue[uplow].get()):
                    self._layers_analysis_indices[uplow].append(group)
                    self._layers_origin_indices[uplow].append(self.analysis_atoms_ids[group])
            for p in proc:
                p.join()
            for q in queue:
                q.close()
        else:
            for index, group in enumerate(self._assign_one_side(up, sort[::-1], _x, _y, _z, _radii)):
                self._layers_analysis_indices[up].append(group)
                self._layers_origin_indices[up].append(self.analysis_atoms_ids[group])
            for index, group in enumerate(self._assign_one_side(low, sort[::], _x, _y, _z, _radii)):
                self._layers_analysis_indices[low].append(group)
                self._layers_origin_indices[low].append(self.analysis_atoms_ids[group])
        # label all layers' atoms with corresponding layer number (1~max_layer) and side number (0 for up, 1 for low)
        for uplow in [0,1]:
            for nlayer,layer in enumerate(self._layers_origin_indices[uplow]):
                if layer is None:
                    pass
                else:
                    self.all_atoms_layers[layer] = nlayer+1
                    self.all_atoms_layersside[layer] = uplow


