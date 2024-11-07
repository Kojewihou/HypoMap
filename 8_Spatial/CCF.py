"""[ChatGPT Generated Documentation]
This script provides functions for processing 3D coordinates within the context of the Common Coordinate Framework (CCF). 
It includes utilities to mirror coordinates across a specified midline, and to retrieve anatomical annotations from a 
3D annotation matrix based on the given points' coordinates.

Functions included:
1. `mirror_ccf_coords`: Mirrors 3D coordinates across a specified midline along a chosen axis.
2. `get_ccf_annotation`: Retrieves anatomical annotations from a 3D annotation matrix for a set of given points.

Example Usage:
--------------
# Mirror points across the brain's midline along the z-axis
mirrored_points = mirror_ccf_coords(points, midline=5.7, axis=2, side='below')

# Get anatomical annotations for points
annotations = get_ccf_annotation(points, resolution=25, annotation_matrix=annotation_volume)
"""

import numpy as np
import pandas as pd
from typing import Literal

def mirror_ccf_coords(
    points: np.ndarray,
    midline: float,
    axis: int,
    side: Literal['below', 'above'] = 'below'
) -> np.ndarray:
    """[ChatGPT Generated Documentation]
    Mirrors 3D coordinates across a specified midline along a chosen axis.

    Parameters:
    - points (np.ndarray): Array of shape (N, 3) containing the 3D coordinates to be mirrored.
    - midline (float): The coordinate of the midline across which the points will be mirrored.
    - axis (int): The axis (0 for x, 1 for y, 2 for z) along which the mirroring will occur.
    - side (Literal['below', 'above']): Specifies whether to mirror to the "below" (negative) or "above" (positive) 
      side of the midline. Defaults to "below".

    Returns:
    - np.ndarray: A new array of mirrored points with the same shape as `points`.

    Raises:
    - ValueError: If `side` is not "above" or "below".

    Notes:
    - For CCFv3, the midline on the z-axis is typically assumed to be at z=5.7.
    - This function copies the input points array to avoid modifying the original data.
    """
    # Create a copy of the points array to avoid modifying the original data
    mirrored_points = np.copy(points)
    
    # Calculate offsets from the midline along the specified axis
    offsets = mirrored_points[:, axis] - midline
    
    # Adjust points based on the specified side of the midline
    if side == "below":
        mirrored_points[:, axis] = midline - np.abs(offsets)
    elif side == "above":
        mirrored_points[:, axis] = midline + np.abs(offsets)
    else:
        raise ValueError("Parameter 'side' must be 'above' or 'below' the midline.")
    
    return mirrored_points

def get_ccf_annotation(
    points: np.ndarray,
    resolution: int,
    annotation_matrix: np.ndarray
) -> pd.Categorical:
    """[ChatGPT Generated Documentation]
    Retrieves anatomical annotations from a 3D annotation matrix for a set of given points.

    Parameters:
    - points (np.ndarray or pd.DataFrame): Array of shape (N, 3) containing the 3D coordinates in micrometers.
    - resolution (int): Resolution of the annotation matrix in micrometers per voxel.
    - annotation_matrix (np.ndarray): 3D array representing anatomical annotations.

    Returns:
    - pd.Categorical: Categorical array of annotations for each point.

    Notes:
    - The function scales coordinates to voxel indices based on the provided `resolution`.
    - Coordinates outside the bounds of `annotation_matrix` are assigned a default annotation of 0.
    """
    # Convert DataFrame to numpy array if necessary
    if isinstance(points, pd.DataFrame):
        points = points.values
    
    # Calculate voxel coordinates for each axis
    voxels = ((points * 1000) // resolution).astype(int)
    
    # Check if voxel coordinates are within the bounds of the annotation matrix
    inbounds_mask = (
        np.logical_and(voxels[:, 0] >= 0, voxels[:, 0] < annotation_matrix.shape[0]) &
        np.logical_and(voxels[:, 1] >= 0, voxels[:, 1] < annotation_matrix.shape[1]) &
        np.logical_and(voxels[:, 2] >= 0, voxels[:, 2] < annotation_matrix.shape[2]) 
    )
    
    # Initialize annotations with a default value of 0
    annotations = np.zeros(points.shape[0], dtype=int)
    
    # Assign annotation values only to points within matrix bounds
    annotations[inbounds_mask] = annotation_matrix[
        voxels[inbounds_mask, 0], 
        voxels[inbounds_mask, 1], 
        voxels[inbounds_mask, 2]
    ]

    return pd.Categorical(annotations)
