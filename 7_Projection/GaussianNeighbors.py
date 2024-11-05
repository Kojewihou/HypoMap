"""
This function computes a Gaussian-weighted k-nearest neighbors (k-NN) matrix between a reference and query AnnData object. 
It uses a Gaussian kernel to weight the distances between the neighbors and returns the resulting sparse matrix.

Overview:
- `GaussianNeighbors`: Fits a k-NN model to a reference dataset's embeddings and applies a Gaussian kernel 
  to the distances between the query dataset's neighbors to produce a normalized, sparse matrix.

Dependencies:
- numpy: For numerical operations.
- scipy.sparse: For creating and manipulating sparse matrices.
- sklearn.neighbors: For computing k-nearest neighbors using brute-force search.

Functions:
- GaussianNeighbors(reference: sc.AnnData, query: sc.AnnData, embed: str = 'X_embed', n_neighbors: int = 15) -> sp.csr_matrix: 
  Computes the Gaussian-weighted k-NN matrix between the reference and query datasets.

Example usage:
    gauss_matrix = GaussianNeighbors(reference_adata, query_adata, embed='X_embed', n_neighbors=15)
"""

import numpy as np
import scipy.sparse as sp
from sklearn.neighbors import NearestNeighbors # replace with cuml.neighbors for GPU accelerated KNN for large datasets
import scanpy as sc

def GaussianNeighbors(reference: sc.AnnData, query: sc.AnnData, embed: str = 'X_embed', n_neighbors: int = 15) -> sp.csr_matrix:
    """
    Computes a Gaussian-weighted k-nearest neighbors matrix between a reference and query AnnData object based on their embeddings.

    Parameters
    ----------
    reference : sc.AnnData
        The reference AnnData object, containing the embedding of the reference dataset in its `obsm` attribute.
    
    query : sc.AnnData
        The query AnnData object, containing the embedding of the query dataset in its `obsm` attribute.
    
    embed : str, optional
        The key in the `obsm` attribute where the embeddings are stored. Default is 'X_embed'.
    
    n_neighbors : int, optional
        The number of nearest neighbors to compute for each point in the query dataset. Default is 15.

    Returns
    -------
    sp.csr_matrix
        A normalized sparse matrix of shape (n_query_samples, n_reference_samples) containing the Gaussian-weighted 
        k-nearest neighbors for each point in the query dataset.
    """

    # Initialize and fit the k-NN model on the reference embedding
    query_nn = NearestNeighbors(n_neighbors=n_neighbors, algorithm='brute', metric='euclidean')
    query_nn.fit(reference.obsm[embed])
    
    # Compute the k-nearest neighbors for the query data
    knn_dists, knn_indices = query_nn.kneighbors(query.obsm[embed])
    
    # Compute the standard deviation for each query point's neighbors
    knn_sd = np.sqrt((knn_dists**2).sum(axis=1) / n_neighbors)
    
    # Apply the Gaussian kernel to the neighbor distances
    knn_gauss = np.exp(-knn_dists / (2 / knn_sd[:, np.newaxis])**2)
    
    # Create a sparse matrix from the Gaussian-weighted distances and indices
    indptr = np.arange(0, (n_neighbors * query.shape[0]) + 1, n_neighbors)
    gauss_matrix = sp.csr_matrix((knn_gauss.ravel(), knn_indices.ravel(), indptr), shape=(query.shape[0], reference.shape[0]))
    
    # Normalize the matrix such that the neighbors sum to 1 across each row
    norm_gauss_matrix = (gauss_matrix / gauss_matrix.sum(axis=1)).tocsr()
    
    return norm_gauss_matrix
