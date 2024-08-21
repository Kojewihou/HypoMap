"""
This script defines a function to merge multiple AnnData objects using the Scanpy library.

Overview:
- The `merge` function concatenates a list of AnnData objects, reattaches the `.var` attribute if the variables are consistent,
  and returns the merged AnnData object.

Dependencies:
- scanpy: A library for single-cell data analysis, used for handling AnnData objects.

Functions:
- merge(adatas): Concatenates a list of AnnData objects and reattaches the `.var` attribute if appropriate.

Example usage:
    import scanpy as sc
    merged_data = merge(adatas)
"""

import scanpy as sc
from typing import List
from anndata import AnnData

def merge(adatas: List[AnnData]) -> AnnData:
    """
    Merges a list of AnnData objects by concatenating them. If all the input objects share the same `.var` index,
    the `.var` attribute is reattached from the first AnnData object after concatenation.

    Args:
        adatas (List[AnnData]): A list of AnnData objects to be concatenated.

    Returns:
        AnnData: The merged AnnData object.
    """

    # Concatenate the AnnData objects with an outer join on observations and variables
    merged = sc.concat(adatas, join='outer')

    # Re-attach the .var attribute if all input AnnData objects have the same variable index
    if all(merged.var.index == adatas[0].var.index):
        merged.var = adatas[0].var
        print('Re-attached .var')

    # Non obs fields must be manually reattached similar to .var e.g. uns obsm

    return merged
