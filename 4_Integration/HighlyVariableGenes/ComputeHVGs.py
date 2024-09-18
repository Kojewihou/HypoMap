"""
This script provides a function for calculating the ranking of highly variable genes (HVGs) in single-cell RNA-seq data.

Overview:
- `calculate_hvg_ranking`: Identifies and ranks highly variable genes in an AnnData object, 
  optionally excluding specific genes from the ranking.

Dependencies:
- pandas: For handling and manipulating data frames.
- scanpy: For processing and analyzing single-cell RNA-seq data, specifically for identifying HVGs.

Functions:
- calculate_hvg_ranking(adata: sc.AnnData, excluded_hvgs: str = None) -> pd.DataFrame: 
  Calculates the ranking of HVGs in the provided AnnData object, with an option to exclude certain genes.

Example usage:
    adata = sc.read_h5ad('example_data.h5ad')
    hvg_ranking = calculate_hvg_ranking(adata)
    hvg_ranking_with_exclusions = calculate_hvg_ranking(adata, excluded_hvgs='excluded_genes.csv')
"""

import pandas as pd
import scanpy as sc

def calculate_hvg_ranking(adata: sc.AnnData, excluded_hvgs: str = None) -> pd.DataFrame:
    """
    Identifies and ranks highly variable genes (HVGs) in an AnnData object. This function uses 
    the `scanpy` library to compute the HVGs and ranks them according to their variability. 
    It also provides an option to exclude specific genes from the ranking based on an external file.

    Parameters
    ----------
    adata : sc.AnnData
        Annotated data matrix from the `scanpy` library, containing gene expression data. 
        The function assumes that the data contains a 'Sample' column in `adata.obs` to use as 
        the `batch_key` for identifying HVGs.

    excluded_hvgs : str, optional
        Path to a CSV file containing gene names to be excluded from the HVG ranking. The file 
        should contain a single column without a header, listing the genes to be excluded. 
        If None, no genes are excluded. Default is None.

    Returns
    -------
    pd.DataFrame
        A DataFrame containing the HVG information, including the `highly_variable_rank` 
        column that ranks genes based on their variability. If `excluded_hvgs` is provided, 
        excluded genes are assigned a `None` rank.

    Notes
    -----
    - The function uses the 'seurat_v3_paper' flavor in `scanpy` for identifying HVGs.
    - The ranking process considers all genes in the `adata` object by setting `n_top_genes` to `adata.shape[1]`.
    - The function modifies the rank of excluded genes to `None` if `excluded_hvgs` is provided.

    Example
    -------
    >>> adata = sc.read_h5ad('example_data.h5ad')
    >>> hvg_ranking = calculate_hvg_ranking(adata)
    >>> hvg_ranking_with_exclusions = calculate_hvg_ranking(adata, excluded_hvgs='excluded_genes.csv')
    """

    # Identify and rank highly variable genes using scanpy
    hvg_df = sc.pp.highly_variable_genes(adata, 
                                         inplace=False, 
                                         n_top_genes=adata.shape[1],
                                         batch_key='Sample',
                                         flavor='seurat_v3_paper',
                                         span=1)
    
    # Sort the DataFrame by the 'highly_variable_rank' column
    hvg_df = hvg_df.sort_values('highly_variable_rank')

    # If excluded HVGs are provided, read and process the exclusion list
    if excluded_hvgs:
        excluded_features = pd.read_csv(excluded_hvgs, header=None)[0].to_list()
        
        # Identify genes to exclude and set their rank to None
        excluded_mask = hvg_df.index.isin(excluded_features)
        hvg_df.loc[excluded_mask, 'highly_variable_rank'] = None
        
        # Re-sort the DataFrame after modification
        hvg_df = hvg_df.sort_values('highly_variable_rank')

    return hvg_df
