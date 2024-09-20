"""
This script defines two key functions for generating UMAP embeddings using RAPIDS and combining UMAP embeddings 
from reference and query datasets for visualization. The `get_umap` function computes UMAP on single-cell data 
and returns the UMAP object. The `add_query` function combines UMAP embeddings from reference and query datasets 
into a single object for visualization.

Overview:
- `get_umap`: Computes UMAP embeddings using RAPIDS for GPU-accelerated dimensionality reduction and returns the 
  UMAP object.
- `add_query`: Combines the UMAP embeddings of reference and query datasets into a single AnnData object for 
  combined visualization.

Dependencies:
- cupy: For GPU-accelerated computation.
- cuml: For RAPIDS UMAP computation.
- cupyx.scipy.sparse: For handling sparse matrices.
- scanpy: For handling single-cell data with AnnData objects.
- sklearn: For managing random state.
- joblib: For saving the UMAP object to disk.

Functions:
- get_umap(adata: AnnData, ..., copy: bool = False, neighbors_key: str | None = None) -> UMAP: 
  Computes UMAP embedding for the input AnnData object.
- add_query(reference: AnnData, query: AnnData, embed: str = 'X_embedded') -> AnnData: 
  Combines UMAP embeddings from reference and query datasets into a single object for joint plotting.

Example usage:
    umap = get_umap(adata, random_state=0, spread=1, min_dist=.05, init_pos='random')
    joblib.dump(umap, 'umap.pkl')
    
    combined_adata = add_query(reference_adata, query_adata)
"""

import cupy as cp
import numpy as np
from cuml import UMAP
from cuml.manifold.umap_utils import find_ab_params
from cupyx.scipy import sparse
from scanpy._utils import NeighborsView
from sklearn.utils import check_random_state
from rapids_singlecell.tools._utils import _choose_representation
from anndata import AnnData
import joblib
import pandas as pd

# Adapted version of RAPIDS UMAP - returns UMAP object
def get_umap(
    adata: AnnData,
    *,
    min_dist: float = 0.5,
    spread: float = 1.0,
    n_components: int = 2,
    maxiter: int | None = None,
    alpha: float = 1.0,
    negative_sample_rate: int = 5,
    init_pos="auto",
    random_state=0,
    a: float | None = None,
    b: float | None = None,
    copy: bool = False,
    neighbors_key: str | None = None,
):
    """
    Computes UMAP embeddings on the input AnnData object using RAPIDS for GPU acceleration. The function supports
    multiple UMAP parameters and leverages precomputed neighbors if available.

    Parameters
    ----------
    adata : AnnData
        The input AnnData object containing the data to be embedded.
    
    min_dist : float, optional
        The minimum distance between points in the low-dimensional UMAP space. Default is 0.5.
    
    spread : float, optional
        The effective scale of the embedded points. Default is 1.0.
    
    n_components : int, optional
        The number of dimensions in the UMAP embedding space. Default is 2.
    
    maxiter : int, optional
        The maximum number of iterations for optimizing the UMAP embedding. Default is None.
    
    alpha : float, optional
        The learning rate for UMAP optimization. Default is 1.0.
    
    negative_sample_rate : int, optional
        The rate of negative sampling used in UMAP optimization. Default is 5.
    
    init_pos : str, optional
        Initialization method for UMAP embedding ('random' or 'spectral'). Default is 'auto'.
    
    random_state : int, optional
        Seed for random number generation. Default is 0.
    
    a : float, optional
        Precomputed UMAP parameter a. If None, it's calculated based on spread and min_dist.
    
    b : float, optional
        Precomputed UMAP parameter b. If None, it's calculated based on spread and min_dist.
    
    copy : bool, optional
        If True, a copy of the AnnData object is made. Default is False.
    
    neighbors_key : str, optional
        The key in adata.uns where precomputed neighbors are stored. Default is None.

    Returns
    -------
    UMAP
        UMAP object after fitting the input data.
    """
    
    adata = adata.copy() if copy else adata

    if neighbors_key is None:
        neighbors_key = "neighbors"

    if neighbors_key not in adata.uns:
        raise ValueError(
            f'Did not find .uns["{neighbors_key}"]. Run `sc.pp.neighbors` first.'
        )

    neighbors = NeighborsView(adata, neighbors_key)

    if a is None or b is None:
        a, b = find_ab_params(spread, min_dist)
    adata.uns["umap"] = {"params": {"a": a, "b": b}}

    if random_state != 0:
        adata.uns["umap"]["params"]["random_state"] = random_state
    random_state = check_random_state(random_state)

    neigh_params = neighbors["params"]
    X = _choose_representation(
        adata,
        neigh_params.get("use_rep", None),
        neigh_params.get("n_pcs", None),
    )

    n_neighbors = neighbors["params"]["n_neighbors"]
    n_epochs = 500 if maxiter is None else maxiter  # Set default if maxiter is not provided
    metric = neigh_params.get("metric", "euclidean")

    if isinstance(X, cp.ndarray):
        X_contiguous = cp.ascontiguousarray(X, dtype=np.float32)
    elif isinstance(X, sparse.spmatrix):
        X_contiguous = X
    else:
        X_contiguous = np.ascontiguousarray(X, dtype=np.float32)

    n_obs = adata.shape[0]
    if neigh_params.get("method") == "rapids":
        knn_dist = neighbors["distances"].data.reshape(n_obs, n_neighbors)
        knn_indices = neighbors["distances"].indices.reshape(n_obs, n_neighbors)
        pre_knn = (knn_indices, knn_dist)
    else:
        pre_knn = None

    if init_pos == "auto":
        init_pos = "spectral" if n_obs < 1000000 else "random"

    umap = UMAP(
        n_neighbors=n_neighbors,
        n_components=n_components,
        metric=metric,
        n_epochs=n_epochs,
        learning_rate=alpha,
        init=init_pos,
        min_dist=min_dist,
        spread=spread,
        negative_sample_rate=negative_sample_rate,
        a=a,
        b=b,
        random_state=random_state,
        output_type="numpy",
        precomputed_knn=None,
    )

    adata.obsm["X_umap"] = umap.fit_transform(X_contiguous)
    
    return umap

# Example use
umap = get_umap(adata, random_state=0, spread=1, min_dist=.05, init_pos='random')

# Save UMAP object to disk
umap_savefile = 'out/path/umap.pkl'
joblib.dump(umap, umap_savefile)

# Load UMAP from disk
umap = joblib.load('out/path/umap.pkl')

# Transforming query into reference UMAP embedding
X_embed = query.obsm['X_embed']

chunk_size = 500_000 # UMAP has to be chunked due to memory errors
transformed_chunks = []
num_chunks = int(np.ceil(X_embed.shape[0] / chunk_size))

for i in range(num_chunks):
    start_idx = i * chunk_size
    end_idx = min((i + 1) * chunk_size, X_embed.shape[0])
    chunk = X_embed[start_idx:end_idx]

    transformed_chunk = UMAP.transform(chunk)
    transformed_chunks.append(transformed_chunk)
    
    torch.cuda.empty_cache() # Memory is sometimes not released correctly from GPU

query.obsm['X_umap'] = np.vstack(transformed_chunks)



def add_query(reference: AnnData, query: AnnData, embed: str = 'X_embedded') -> AnnData:
    """
    Combines the UMAP embeddings of a reference and query AnnData object for joint visualization.

    Parameters
    ----------
    reference : AnnData
        The reference AnnData object containing the original data and UMAP embeddings.
    
    query : AnnData
        The query AnnData object containing new data to be visualized alongside the reference data.
    
    embed : str, optional
        The key in `obsm` where the embeddings are stored. Default is 'X_embedded'.

    Returns
    -------
    AnnData
        A combined AnnData object containing both reference and query data with concatenated embeddings.
    """
    
    # Create a new AnnData object to hold combined reference and query data
    plotting_adata = AnnData(shape=(reference.shape[0] + query.shape[0], 1))
    plotting_adata.obs = pd.concat([reference.obs, query.obs])
    plotting_adata.obsm[embed] = np.vstack((reference.obsm[embed], query.obsm[embed]))
    
    # Combine UMAP embeddings if available
    if 'X_umap' in reference.obsm_keys() and 'X_umap' in query.obsm_keys():
        plotting_adata.obsm['X_umap'] = np.vstack((reference.obsm['X_umap'], query.obsm['X_umap']))
    
    return plotting_adata
