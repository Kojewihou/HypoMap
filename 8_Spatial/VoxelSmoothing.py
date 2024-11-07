"""[ChatGPT Generated Documentation]
This script provides functions for aggregating and smoothing voxel data based on 3D point clouds and their associated scores. 
It includes utilities to divide spatial data into voxels, aggregate scores within each voxel, apply smoothing, and normalize 
scores along specific axes.

Functions included:
1. `aggregate_voxels`: Aggregates points and scores into a 3D voxel grid.
2. `smooth_across_voxels`: Applies Gaussian smoothing across the voxel grid, with options to ignore certain voxels.
3. `normalize_scores_axis`: Normalizes scores across a specified axis of the voxel grid.

Example Usage:
--------------
# Aggregate voxels from points and scores
mean_scores, counts = aggregate_voxels(points, scores, voxel_size=5)

# Smooth the aggregated voxel matrix
smoothed_voxels = smooth_across_voxels(mean_scores, sigma=1.0, score_threshold=0.1)

# Normalize scores across the last axis
normalized_voxels = normalize_scores_axis(smoothed_voxels)
"""

import numpy as np
from typing import Tuple, Optional
from numpy.typing import NDArray
from scipy.ndimage import gaussian_filter

def aggregate_voxels(
    points: NDArray[np.float64], 
    scores: NDArray[np.float64], 
    voxel_size: int
) -> Tuple[NDArray[np.float64], NDArray[np.int64]]:
    """[ChatGPT Generated Documentation]
    Aggregates points and scores into a 3D voxel grid, computing the average score in each voxel.

    Parameters:
    - points (NDArray[np.float64]): An array of shape (N, 3) containing 3D coordinates of the points.
    - scores (NDArray[np.float64]): An array of shape (N,) or (N, M) containing the scores associated with each point. 
    - voxel_size (int): Size of each voxel.

    Returns:
    - Tuple[NDArray[np.float64], NDArray[np.int64]]: 
        - `mean_scores_grid` (NDArray[np.float64]): A 3D grid of mean scores within each voxel, shape (X, Y, Z, M) or (X, Y, Z).
        - `counts_grid` (NDArray[np.int64]): A 3D grid of point counts within each voxel, shape (X, Y, Z).

    Notes:
    - If `scores` has only one metric, the output `mean_scores_grid` will have the last dimension removed for simplicity.
    - This function uses numpy's `add.at` for efficient aggregation over large voxel grids.
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
    sigma: float,
    ignore_zero_matrices: bool=True,
    score_threshold: Optional[float] = None,
) -> NDArray[np.float64]:
    """[ChatGPT Generated Documentation]
    Applies Gaussian smoothing across a voxel grid, with options to threshold and ignore certain matrices.

    Parameters:
    - voxel_matrix (NDArray[np.float64]): A 3D or 4D grid of scores representing the voxel data.
    - sigma (float): Standard deviation for the Gaussian kernel.
    - ignore_zero_matrices (bool): If True, applies smoothing only to non-zero matrices.
    - score_threshold (Optional[float]): Optional threshold; values below are set to zero before smoothing.

    Returns:
    - NDArray[np.float64]: The smoothed voxel grid with the same shape as `voxel_matrix`.

    Notes:
    - The Gaussian smoothing is applied along the first three dimensions of the voxel grid.
    - `ignore_zero_matrices` avoids smoothing over areas with no valid data if set to True.
    """
    smoothed_voxel_matrix = voxel_matrix.copy()
    
    # Optionally reduce noise by setting any score below a threshold to 0
    if score_threshold:
        smoothed_voxel_matrix[smoothed_voxel_matrix < score_threshold] = 0
    
    if ignore_zero_matrices: # Untested
        mask = np.count_nonzero(smoothed_voxel_matrix, axis=(0, 1, 2)) != 0
        smoothed_voxel_matrix[mask] = gaussian_filter(smoothed_voxel_matrix[mask], sigma=sigma, truncate=3, mode='constant', axes=[0, 1, 2])
    else:
        smoothed_voxel_matrix = gaussian_filter(smoothed_voxel_matrix, sigma=sigma, truncate=3, mode='constant', axes=[0, 1, 2])
    
    return smoothed_voxel_matrix

def normalize_scores_axis(voxel_matrix: NDArray[np.float64]) -> NDArray[np.float64]:
    """[ChatGPT Generated Documentation]
    Normalizes scores along the score axis (3) of the voxel grid.

    Parameters:
    - voxel_matrix (NDArray[np.float64]): A 4D grid of scores representing the voxel data.

    Returns:
    - NDArray[np.float64]: The voxel grid with normalized scores along the score axis (3).

    Notes:
    - The function prevents division by zero by normalizing only where score totals are positive.
    """
    
    assert voxel_matrix.ndim == 4, f'voxel_matrix is not 4D - (ndim={voxel_matrix.ndim})' # Untested
    
    score_totals = np.sum(voxel_matrix, axis=3, keepdims=True)
    return np.divide(voxel_matrix, score_totals, where=(score_totals > 0))
