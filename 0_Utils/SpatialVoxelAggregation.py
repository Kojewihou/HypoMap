from typing import Tuple
from numpy.typing import NDArray

def aggregate_voxels(
    points: NDArray,
    scores: NDArray,
    voxel_size: int
) -> Tuple[NDArray, NDArray]:
    """[ChatGPT Generated Docstring]
    Aggregates scores within a 3D voxel grid based on given points, calculating mean scores for each voxel.

    Given a set of 3D points, this function groups the points into voxels of a specified size and aggregates the
    corresponding scores within each voxel. If multiple score metrics are provided, the mean score for each metric 
    is calculated independently per voxel.

    Parameters
    ----------
    points : NDArray
        An (N, 3) array of 3D coordinates, where N is the number of points.
    scores : NDArray
        An (N,) array of scores associated with each point, or an (N, M) array where each column represents a 
        different score metric.
    voxel_size : int
        The edge length of each voxel, which determines the binning resolution for the points.

    Returns
    -------
    mean_scores_grid : NDArray
        A 3D array of mean scores per voxel if there is a single metric, or a 4D array (X, Y, Z, M) if there are 
        multiple metrics. Each voxel contains the average score for all points within it.
    counts_grid : NDArray
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
    
    # Calculate the range for bins in each dimension
    x_bins = np.arange(np.min(points[:, 0]) // voxel_size, np.max(points[:, 0]), step=voxel_size)
    y_bins = np.arange(np.min(points[:, 1]) // voxel_size, np.max(points[:, 1]), step=voxel_size)
    z_bins = np.arange(np.min(points[:, 2]) // voxel_size, np.max(points[:, 2]), step=voxel_size)

    # Determine voxel indices for each point
    x_indices = np.digitize(points[:, 0], x_bins) - 1
    y_indices = np.digitize(points[:, 1], y_bins) - 1
    z_indices = np.digitize(points[:, 2], z_bins) - 1

    # Define the shape of the 3D grid based on bin sizes
    grid_shape = (x_bins.size, y_bins.size, z_bins.size)
    num_metrics = scores.shape[1]  # Number of score metrics

    # Initialize grids for counts and scores
    counts_grid = np.zeros(grid_shape, dtype=int)
    scores_grid = np.zeros((*grid_shape, num_metrics), dtype=float)

    # Aggregate count and scores for each voxel and each score metric
    np.add.at(counts_grid, (x_indices, y_indices, z_indices), 1)
    for metric in range(num_metrics):
        np.add.at(scores_grid, (x_indices, y_indices, z_indices, metric), scores[:, metric])

    # Calculate mean scores in each voxel where count > 0 for each metric
    mean_scores_grid = np.zeros_like(scores_grid, dtype=float)
    occupied_voxels = counts_grid > 0
    mean_scores_grid[occupied_voxels] = scores_grid[occupied_voxels] / counts_grid[occupied_voxels, None]

    # If only one score metric, remove the last dimension for simplicity
    if num_metrics == 1:
        mean_scores_grid = mean_scores_grid[..., 0]
    
    return mean_scores_grid, counts_grid


# Gaussian Filter Output
# smoothed = gaussian_filter(mean_scores_grid, sigma=6, axes=[0,1,2])