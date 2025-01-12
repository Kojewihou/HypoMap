from tqdm import tqdm

import numpy as np
import scanpy as sc
from ccHBGF import ccHBGF

from ._IGraphLeiden import igraph_leiden

def decompose(adata, cluster_obs, resolution, n_random_states=10):
    
    decomposed_clusters = adata.obs[cluster_obs].astype('object').copy()
    
    for cluster in tqdm(adata.obs[cluster_obs].unique(), desc='Decomposing Clusters'):
        mask = adata.obs[cluster_obs] == cluster
        subset = adata[mask]
        
        sc.pp.neighbors(subset, n_neighbors=15, use_rep='X_embed', random_state=0)
        
        igraph_subset = sc._utils.get_igraph_from_adjacency(subset.obsp['connectivities'])
        
        labelmatrix = np.zeros((subset.shape[0], n_random_states))
        for seed in range(n_random_states):
            labelmatrix[:, seed] = igraph_leiden(igraph_subset, resolution=resolution, random_state=seed)
        
        labels = ccHBGF(labelmatrix, random_state=0, verbose=False)
        prefixed_labels = [f"{cluster}_{label}" for label in labels]
        
        decomposed_clusters.loc[mask] = prefixed_labels
    
    return decomposed_clusters.astype('category')