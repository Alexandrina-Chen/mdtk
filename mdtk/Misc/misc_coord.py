
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