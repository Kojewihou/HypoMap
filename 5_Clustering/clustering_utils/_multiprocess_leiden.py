"""
[ChatGPT Generated Documentation]
Script for Leiden Clustering on an igraph Object

Overview:
    This script defines a single function, `igraph_leiden`, which performs Leiden clustering on an igraph object.
    It adapts the `scanpy.tl.leiden` function to allow for multiprocessing by removing the requirement for an AnnData object.

Key Function:
    - igraph_leiden: Executes the Leiden clustering algorithm with customizable resolution and iteration parameters.

Example Usage:
    ```python
    import igraph as ig
    from scanpy import _utils
    import pandas as pd

    # Example graph creation
    g = ig.Graph.Famous("Zachary")  # Use a famous graph dataset (e.g., Zachary's Karate Club)

    # Running the Leiden clustering
    membership = igraph_leiden(g, resolution=1.0, n_iterations=2)

    # Printing the membership
    print(membership)
    ```
    
Requirements:
    - pandas
    - scanpy

"""

import pandas as pd
from scanpy import _utils
from igraph import Graph
from typing import List, Dict, Any

def igraph_leiden(
    g: Graph,
    resolution: float = 1.0,
    n_iterations: int = 2,
    random_state: _utils.AnyRandom = 0,
    **clustering_args: Dict[str, Any]
) -> List[int]:
    """
    [ChatGPT Generated Documentation]
    Performs Leiden clustering on an igraph object with customizable parameters.

    Parameters:
    ----------
    g : igraph.Graph
        An igraph object representing the graph to be clustered.
        Typically derived using `sc._utils.get_igraph_from_adjacency(adata.obsp[connectivities_key])`
        after running `sc.pp.neighbors` on AnnData.

    resolution : float, optional (default=1.0)
        The resolution parameter controls the coarseness of the clustering. Higher values lead to more clusters.

    n_iterations : int, optional (default=2)
        The number of iterations for the Leiden algorithm. More iterations may result in better clustering at the cost of speed.

    random_state : _utils.AnyRandom, optional (default=0)
        The random seed or state for reproducibility. It can be an integer, a `numpy.random.RandomState` instance, or `None`.

    **clustering_args : Dict[str, Any]
        Additional keyword arguments passed directly to the `community_leiden` method of the igraph object.
        Common arguments include `objective_function`, which by default is set to "modularity".

    Returns:
    -------
    List[int]
        A list of integers representing the cluster membership of each node in the graph.

    Notes:
    ------
    - The function adapts `scanpy.tl.leiden` for use without an AnnData object, making it suitable for multiprocessing contexts.
    - Default clustering objective is set to "modularity" unless specified otherwise in `clustering_args`.
    - The function ensures reproducibility by setting the igraph random state using `scanpy._utils.set_igraph_random_state`.
    
    Example:
    --------
    ```python
    import igraph as ig
    from scanpy import _utils

    g = ig.Graph.Famous("Zachary")
    membership = igraph_leiden(g, resolution=0.5, n_iterations=3)
    print(membership)
    ```
    """
    
    # Copy clustering arguments and set number of iterations
    clustering_args = dict(clustering_args)
    clustering_args["n_iterations"] = n_iterations
    
    # Set resolution if provided
    if resolution is not None:
        clustering_args["resolution"] = resolution

    # Set default objective function if not specified
    clustering_args.setdefault("objective_function", "modularity")
    
    # Run Leiden clustering with specified random state
    with _utils.set_igraph_random_state(random_state):
        part = g.community_leiden(**clustering_args)

    return part.membership
