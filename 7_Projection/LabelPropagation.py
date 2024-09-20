"""
This function propagates labels from a reference dataset to a query dataset using a Gaussian-weighted k-nearest neighbors 
approach. It computes the label probabilities for each query point based on the neighbors in the reference dataset and 
assigns the most probable label along with the confidence score.

Overview:
- `PropagateLabels`: Propagates labels from a reference dataset to a query dataset by leveraging Gaussian-weighted k-NN.
  The function calculates label probabilities for the query dataset and returns the predicted labels and their confidence.

Dependencies:
- GaussianNeighbors: For computing the Gaussian-weighted k-nearest neighbors.
- pandas: For handling label assignments and confidence scores.
- scipy.sparse: For efficient storage and manipulation of sparse matrices.

Functions:
- PropagateLabels(reference: sc.AnnData, query: sc.AnnData, label: str, embed: str = 'X_embed', n_neighbors: int = 15) 
  -> Tuple[pd.Series, pd.Series]: Propagates labels from the reference to the query dataset and returns predicted labels 
  and their confidence scores.

Example usage:
    labels, confidence = PropagateLabels(reference, query, 'major_features_propagated')
    labels[confidence < 0.65] = None  # Remove labels with low confidence
    query.obs['major_features_projected'] = labels
    query.obs['major_features_projected_confidence'] = confidence
"""

import numpy as np
import pandas as pd
import scipy.sparse as sp
from GaussianNeighbors import GaussianNeighbors
import scanpy as sc
from typing import Tuple

def PropagateLabels(reference: sc.AnnData, query: sc.AnnData, label: str, embed: str = 'X_embed', n_neighbors: int = 15) -> Tuple[pd.Series, pd.Series]:
    """
    Propagates labels from a reference dataset to a query dataset using Gaussian-weighted k-nearest neighbors. 
    The function calculates label probabilities based on the neighbors in the reference dataset and assigns 
    the most probable label to each query point.

    Parameters
    ----------
    reference : sc.AnnData
        The reference AnnData object containing the labeled data.
    
    query : sc.AnnData
        The query AnnData object to which labels will be propagated.
    
    label : str
        The column name in `reference.obs` that contains the labels to propagate.
    
    embed : str, optional
        The key in `obsm` where the embeddings are stored. Default is 'X_embed'.
    
    n_neighbors : int, optional
        The number of neighbors to use for k-nearest neighbors. Default is 15.

    Returns
    -------
    Tuple[pd.Series, pd.Series]
        A tuple containing:
        - Predicted labels as a pandas Series, indexed by the query's observation indices.
        - Confidence scores for the predicted labels as a pandas Series, indexed by the query's observation indices.
    """

    # Compute Gaussian-weighted nearest neighbors matrix
    gauss_mtx = GaussianNeighbors(reference, query, embed, n_neighbors)
    
    # Initialize a sparse matrix to store label probabilities
    unique_labels = reference.obs[label].unique()
    label_probs = sp.lil_matrix((query.shape[0], len(unique_labels)))

    # Calculate label probabilities by summing Gaussian weights for each label
    for idx, u_label in enumerate(unique_labels):
        cluster_mask = reference.obs[label] == u_label
        label_probs[:, idx] = gauss_mtx[:, cluster_mask].sum(axis=1).A1

    # Convert label probabilities to CSR format
    query_label_probs = label_probs.tocsr()
    
    # Find the label with the maximum probability for each query point
    top_label_idx = query_label_probs.argmax(axis=1).A1
    label_confidence = query_label_probs.max(axis=1).A.ravel()
    query_label = unique_labels[top_label_idx]
    
    # Convert labels and confidence scores to pandas Series
    query_label = pd.Series(query_label, index=query.obs.index, dtype="category")
    label_confidence = pd.Series(label_confidence, index=query.obs.index, dtype='float32')

    return query_label, label_confidence


# Example usage 1 - Propagate reference labels to query
labels, confidence = PropagateLabels(reference, query, 'label')

# Filter out labels with low confidence
labels[confidence < 0.65] = None

# Store the labels and confidence in the query AnnData object
query.obs['labels_projected'] = labels
query.obs['labels_projected_confidence'] = confidence


# Example usage 2 - Propagate partial reference labels across the all reference
mask = adata.obs['partial_label'].isna()

# Set reference to partial label subset and query to unlabelled subset
unlabelled_view = adata[mask]
labelled_view = adata[~mask]

# Propogate the same way
propagated_labels, propagated_confidence = PropagateLabels(labelled_view, unlabelled_view, label='partial_label')

# Comment out if all cells are to labelled regardless of confidence
propagated_labels[propagated_confidence < .5] = None

# Return labels to full reference
adata.obs['label_projected'] = adata.obs['partial_label']
adata.obs['label_projected_confidence'] = pd.Series(None, dtype=np.float32)

adata.obs.loc[propagated_labels.index, 'label_projected'] = propagated_labels
adata.obs.loc[propagated_labels.index, 'label_projected_confidence'] = propagated_confidence