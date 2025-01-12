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
    connectivities_key = adata.uns['neighbors']['connectivities_key']
    adata_igraph = sc._utils.get_igraph_from_adjacency(adata.obsp[connectivities_key])
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
