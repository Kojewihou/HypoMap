"""[ChatGPT Generated Documentation]
This script provides an adapted Leiden clustering function optimized for multiprocessing and designed to be used without requiring an AnnData object from Scanpy. 

The function leverages the igraph-based Leiden community detection algorithm, with parameters to control clustering resolution and iteration count. 

Functions included:
1. `igraph_leiden`: Performs Leiden clustering on an adjacency matrix, returning cluster labels as a pandas Categorical type.

Example Usage:
--------------
# Run Leiden clustering on a sparse adjacency matrix with a specific resolution
cluster_labels = igraph_leiden(adjacency=adj_matrix, resolution=0.5, n_iterations=3)
"""

import pandas as pd
import numpy as np
from scanpy import _utils
from scipy import sparse
from typing import Optional

def igraph_leiden(
    adjacency: Optional[sparse.spmatrix] = None,
    resolution: float = 1.0,
    n_iterations: int = 2,
    random_state: _utils.AnyRandom = 0,
    **clustering_args,
) -> pd.Categorical:
    """[ChatGPT Generated Documentation]
    Performs Leiden clustering on a given adjacency matrix, with support for additional customization through keyword arguments.

    Parameters:
    - adjacency (Optional[sparse.spmatrix]): Sparse matrix of graph connectivities, such as `.obsp` from an AnnData object.
    - resolution (float): Resolution parameter for the Leiden algorithm, controlling the coarseness of the clustering. Defaults to 1.0.
    - n_iterations (int): Number of optimization iterations for the clustering algorithm. Defaults to 2.
    - random_state (_utils.AnyRandom): Seed or random state for reproducibility in the Leiden algorithm. Defaults to 0.
    - **clustering_args: Additional arguments passed directly to the `community_leiden` method.

    Returns:
    - pd.Categorical: A categorical array representing cluster labels for each node.

    Notes:
    - The function adapts `scanpy.tl.leiden` for multiprocessing by removing the requirement for an AnnData object.
    - Default clustering objective is set to "modularity" unless specified otherwise in `clustering_args`.
    """
    # Copy clustering arguments and set number of iterations
    clustering_args = dict(clustering_args)
    clustering_args["n_iterations"] = n_iterations

    # Convert adjacency to igraph object
    g = _utils.get_igraph_from_adjacency(adjacency, directed=False)
    
    # Set resolution if provided
    if resolution is not None:
        clustering_args["resolution"] = resolution

    # Set default objective function if not specified
    clustering_args.setdefault("objective_function", "modularity")
    
    # Run Leiden clustering with specified random state
    with _utils.set_igraph_random_state(random_state):
        part = g.community_leiden(**clustering_args)
        
    # Retrieve cluster labels as an array and convert to pandas Categorical
    groups = np.array(part.membership)
    membership = pd.Categorical(groups.astype("U"))

    return membership
