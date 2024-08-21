"""
This script provides functionality to reformat an AnnData object that was mapped to an incorrect genome assembly. 
The `convert_to_GRCm39` function remaps the dataset to align with the reference genome assembly used in a template file.

Overview:
- The script loads a template AnnData object containing the correct genome assembly.
- The `convert_to_GRCm39` function creates a new AnnData object, reformatting the input AnnData to match the template's genome assembly.

Dependencies:
- numpy: For numerical operations and array handling.
- pandas: For handling indexed data.
- scanpy: For working with single-cell data in AnnData format.
- scipy: For sparse matrix operations.

Functions:
- convert_to_GRCm39(adata): Remaps the input AnnData object to the correct genome assembly based on a template.

Example usage:
    adata_corrected = convert_to_GRCm39(adata)
"""

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import lil_matrix, csr_matrix

# Load the template file with the correct genome assembly
template_file = '../results/filtered/v1-03.h5ad'
template = sc.read_h5ad(template_file, backed='r')

def convert_to_GRCm39(adata: sc.AnnData) -> sc.AnnData:
    """
    Converts an AnnData object mapped to an incorrect genome assembly to match the genome assembly of a template AnnData object.

    Args:
        adata (sc.AnnData): An AnnData object that needs to be remapped to the correct genome assembly.

    Returns:
        sc.AnnData: A new AnnData object with the data reformatted to match the genome assembly of the template.
    """

    print("Incorrect genome build identified: remapping to match template file")

    # Create a new AnnData object with the correct shape and variable names
    X_shape = (adata.shape[0], template.shape[1])
    reformatted_adata = sc.AnnData(
        X=lil_matrix(X_shape, dtype=np.float32),  # Initialize with sparse matrix
        obs=adata.obs, 
        var=template.var
    )
    
    # Get the gene IDs from both the input AnnData and the template
    adata_gene_ids = pd.Index(adata.var['gene_ids'].values)
    reformatted_adata_gene_ids = pd.Index(reformatted_adata.var['gene_ids'].values)

    # Find common genes between the input AnnData and the template
    common_genes = adata_gene_ids.intersection(reformatted_adata_gene_ids)

    # Get indices of common genes in both the input and reformatted AnnData objects
    common_gene_indices_data = [np.where(adata_gene_ids == gene)[0][0] for gene in common_genes]
    common_gene_indices_reformatted_data = [np.where(reformatted_adata_gene_ids == gene)[0][0] for gene in common_genes]

    # Map the data from the input AnnData to the reformatted AnnData based on common genes
    reformatted_adata[:, common_gene_indices_reformatted_data].X = adata[:, common_gene_indices_data].X

    # Convert the X matrix to CSR format for efficiency
    adata = reformatted_adata
    adata.X = csr_matrix(adata.X)

    return adata
