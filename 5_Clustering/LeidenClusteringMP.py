import pandas as pd
import numpy as np
from scanpy import _utils   
from scipy import sparse

# An adapted leiden function for multiprocessing
def igraph_leiden(
    adjacency: sparse.spmatrix | None = None,
    resolution: float = 1,
    n_iterations: int = 2,
    random_state: _utils.AnyRandom = 0,
    **clustering_args,
) -> pd.Categorical:
    """Leiden clustering functioned adapted from scanpy.tl.leiden - requiring no adata object. 
    Designed for use with multiprocessing.

    Args:
        adjacency (sparse.spmatrix | None, optional): Connectivies .obsp from adata obj. Defaults to None.
        resolution (float, optional): Resolution for leiden clustering. Defaults to 1.
        n_iterations (int, optional): Number of optimzation iterations. Defaults to 2.
        random_state (_utils.AnyRandom, optional): Random state for leiden algorithm. Defaults to 0.

    Returns:
        pd.Categorical: Leiden Cluster Labels
    """

    clustering_args = dict(clustering_args)
    clustering_args["n_iterations"] = n_iterations

    g = _utils.get_igraph_from_adjacency(adjacency, directed=False)
    
    if resolution is not None:
        clustering_args["resolution"] = resolution

    clustering_args.setdefault("objective_function", "modularity")
    
    with _utils.set_igraph_random_state(random_state):
        part = g.community_leiden(**clustering_args)
        
    groups = np.array(part.membership)
    
    membership = pd.Categorical(groups.astype("U"))

    return membership