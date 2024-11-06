"""[ChatGPT Generated Documentation]
Voxel Aggregation and Smoothing Script

This script provides two main functions for processing 3D spatial data in a voxel grid format:
1. `aggregate_voxels`: Aggregates scores within a voxel grid by calculating the mean score for each voxel based on 3D point locations.
2. `smooth_across_voxels`: Applies Gaussian smoothing to a voxel matrix, creating a continuous field from discrete voxel data.

Example Usage:
--------------
points = np.array([[1.5, 2.0, 3.0], [4.0, 5.0, 6.0], ...])  # (N, 3) array of 3D points
scores = np.array([0.8, 1.2, ...])  # (N,) array of scores
voxel_size = 2

mean_scores_grid, counts_grid = aggregate_voxels(points, scores, voxel_size)
smoothed_scores_grid = smooth_across_voxels(mean_scores_grid, sigma=1.0)

Functions:
----------
- aggregate_voxels(points: NDArray, scores: NDArray, voxel_size: int) -> Tuple[NDArray, NDArray]
- smooth_across_voxels(voxel_matrix: NDArray, sigma: float) -> NDArray
"""

import numpy as np
from typing import Tuple
from numpy.typing import NDArray
from scipy.ndimage import gaussian_filter

def aggregate_voxels(
    points: NDArray[np.float64], 
    scores: NDArray[np.float64], 
    voxel_size: int
) -> Tuple[NDArray[np.float64], NDArray[np.int64]]:
    """[ChatGPT Generated Documentation]
    Aggregates scores within a 3D voxel grid based on given points, calculating mean scores for each voxel.

    Given a set of 3D points, this function groups the points into voxels of a specified size and aggregates the
    corresponding scores within each voxel. If multiple score metrics are provided, the mean score for each metric 
    is calculated independently per voxel.

    Parameters
    ----------
    points : NDArray[np.float64]
        An (N, 3) array of 3D coordinates, where N is the number of points.
    scores : NDArray[np.float64]
        An (N,) array of scores associated with each point, or an (N, M) array where each column represents a 
        different score metric.
    voxel_size : int
        The edge length of each voxel, which determines the binning resolution for the points.

    Returns
    -------
    mean_scores_grid : NDArray[np.float64]
        A 3D array of mean scores per voxel if there is a single metric, or a 4D array (X, Y, Z, M) if there are 
        multiple metrics. Each voxel contains the average score for all points within it.
    counts_grid : NDArray[np.int64]
        A 3D array (X, Y, Z) representing the number of points in each voxel.

    Notes
    -----
    - The `points` array should contain 3D coordinates; otherwise, results may be unpredictable.
    - The `scores` array can be either 1D or 2D. A 1D array is reshaped for compatibility, and the function
      handles it as a single score metric.
    - Voxel indices are determined by dividing the point coordinates by `voxel_size`.
    - If there are no points in a voxel, the mean score for that voxel is zero.
    - This function returns the mean score for each voxel, with the shape adjusted based on the number of metrics.
    """
    
    if scores.ndim == 1:
        scores = scores[:, None]  # Reshape (N,) to (N, 1)
        
    num_metrics = scores.shape[1]  # Number of score metrics
    
    min_x, min_y, min_z = (points.min(axis=0) // voxel_size).astype(int)
    
    # Calculate voxel indices for each point, relative to the minimum value
    x_indices = ((points[:, 0] - min_x) // voxel_size).astype(int)
    y_indices = ((points[:, 1] - min_y) // voxel_size).astype(int)
    z_indices = ((points[:, 2] - min_z) // voxel_size).astype(int)

    # Define the shape of the 3D grid based on bin sizes
    grid_shape = (
        x_indices.max() + 1,
        y_indices.max() + 1,
        z_indices.max() + 1
    )

    # Initialize grids for counts and scores
    counts_grid = np.zeros(grid_shape, dtype=int)
    scores_grid = np.zeros((*grid_shape, num_metrics), dtype=float)

    # Aggregate count and scores for each voxel
    np.add.at(counts_grid, (x_indices, y_indices, z_indices), 1)
    np.add.at(scores_grid, (x_indices, y_indices, z_indices, slice(None)), scores)

    # Calculate mean scores in each voxel where count > 0 for each metric
    mean_scores_grid = np.divide(scores_grid, counts_grid[..., None], where=(counts_grid[..., None] > 0))

    # If only one score metric, remove the last dimension for simplicity
    if num_metrics == 1:
        mean_scores_grid = mean_scores_grid[..., 0]
    
    return mean_scores_grid, counts_grid

def smooth_across_voxels(
    voxel_matrix: NDArray[np.float64], 
    sigma: float
) -> NDArray[np.float64]:
    """[ChatGPT Generated Documentation]
    Applies Gaussian smoothing across a voxel matrix, reducing discrete data noise and creating a continuous field.

    Parameters
    ----------
    voxel_matrix : NDArray[np.float64]
        A 3D or 4D voxel matrix to apply the Gaussian smoothing to. Dimensions represent (X, Y, Z) or (X, Y, Z, M).
    sigma : float
        Standard deviation for Gaussian kernel, controlling the degree of smoothing.

    Returns
    -------
    smoothed_matrix : NDArray[np.float64]
        A smoothed version of the voxel matrix, with the same shape as the input.
    """
    return gaussian_filter(voxel_matrix, sigma=sigma, truncate=3, mode='constant', axes=[0, 1, 2])
