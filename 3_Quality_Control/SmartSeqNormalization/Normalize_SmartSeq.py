"""
This script provides a function for normalizing single-cell RNA-seq data, specifically tailored 
for Smart-seq sequencing data, by correcting for gene length biases.

Overview:
- `normalize_smartseq`: Normalizes gene expression data in an AnnData object by adjusting for 
  gene length, accounting for biases inherent in Smart-seq data.

Dependencies:
- numpy: For numerical operations.
- pandas: For handling gene length data from a CSV file.
- scipy.sparse: For efficient storage and manipulation of sparse matrices.
- scanpy: For processing and analyzing single-cell RNA-seq data in AnnData format.

Functions:
- normalize_smartseq(adata: sc.AnnData, gene_length_url: str) -> sc.AnnData: Normalizes the 
  AnnData object's expression values by gene length using data from a provided URL.

Example usage:
    adata = sc.read_h5ad('example_data.h5ad')
    gene_length_url = 'https://example.com/gene_lengths.csv'
    adata = normalize_smartseq(adata, gene_length_url)
"""

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.sparse import csr_matrix
import scanpy as sc

def normalize_smartseq(adata: sc.AnnData, gene_length_url: str) -> sc.AnnData:
    """
    Normalizes gene expression data in an AnnData object for Smart-seq sequencing data 
    by correcting for gene length biases.

    This function reads a CSV file from a specified URL containing gene lengths, joins this 
    information with the `adata.var` DataFrame based on gene IDs, and normalizes the expression 
    matrix (`adata.X`) by dividing each gene's expression by its length. The expression values 
    are then scaled by the median gene length to account for sequencing depth.

    Parameters
    ----------
    adata : sc.AnnData
        Annotated data matrix from the `scanpy` library, containing gene expression data 
        to be normalized. The `adata.var` DataFrame should have a column named 'gene_ids' 
        to allow for the merging with the gene lengths.
    
    gene_length_url : str
        URL to the CSV file containing gene lengths. The file should have no header and 
        must contain gene IDs in the first column and gene lengths in the second column.

    Returns
    -------
    sc.AnnData
        The input AnnData object with normalized expression values in `adata.X`.

    Notes
    -----
    - The function assumes that the AnnData object uses raw expression values.
    - The normalization is done by adjusting the expression levels to account for gene length differences.

    Example
    -------
    >>> adata = sc.read_h5ad('example_data.h5ad')
    >>> gene_length_url = 'https://example.com/gene_lengths.csv'
    >>> adata = normalize_smartseq(adata, gene_length_url)
    """

    # Load the gene lengths from the provided URL
    gene_len = pd.read_csv(
        gene_length_url,
        delimiter=",",
        header=None,
        index_col=0
    )

    # Rename the index and columns to improve readability
    gene_len.index.name = 'gene_ids'
    gene_len = gene_len.rename(columns={1: 'canonical_transcript_length'})

    # Merge the gene length information with the 'var' DataFrame of the AnnData object
    var_temp = adata.var.join(gene_len, on='gene_ids')

    # Normalize the expression matrix by gene lengths and adjust to median gene length
    adata.X = adata.X / var_temp['canonical_transcript_length'].values * np.median(var_temp['canonical_transcript_length'].values)
    
    # Convert the normalized expression matrix to a sparse matrix format
    adata.X = csr_matrix(np.rint(adata.X))

    return adata
