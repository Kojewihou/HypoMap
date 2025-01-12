"""
This script provides a utility function to sanitize an AnnData object's raw count matrix.

Overview:
- The `sanitize_raw_anndata` function modifies the `X` attribute of an AnnData object by converting it to an integer type 
  and eliminating any zero entries in sparse matrices.

Dependencies:
- anndata: A package for handling annotated data matrices, commonly used in single-cell genomics.

Functions:
- sanitize_raw_anndata(adata): Converts the count matrix to integers and removes zeros from sparse matrices.

Example usage:
    import anndata as ad
    sanitize_raw_anndata(adata)
"""

import anndata as ad

def sanitize_anndata(adata: ad.AnnData) -> None:
    """
    Sanitizes the raw count matrix in an AnnData object by converting its data type to integers and 
    eliminating zero entries if the matrix is sparse.

    Args:
        adata (ad.AnnData): An AnnData object containing the data to be sanitized.

    Returns:
        None: The function modifies the AnnData object in place.
    """

    # Convert the count matrix to integer type
    adata.X = adata.X.astype(int)
    
    # Eliminate zero entries if the matrix is sparse
    if hasattr(adata.X, 'eliminate_zeros'):
        adata.X.eliminate_zeros()
