"""
This script performs consensus clustering using the ccHBGF method on an AnnData object, which contains 
Leiden clustering solutions stored in its `obsm` attribute. The consensus clustering results are saved 
in the `obs` attribute of the AnnData object.

Overview:
- `main`: Main function that loads an AnnData object, identifies Leiden clustering matrices, performs 
  consensus clustering using ccHBGF, and saves the results in the AnnData object.

Dependencies:
- argparse: For parsing command-line arguments.
- os: For file and directory checking.
- tqdm: For displaying progress bars.
- ccHBGF: For performing consensus clustering.
- numpy: For numerical operations.
- scanpy: For processing and analyzing single-cell RNA-seq data.

Functions:
- main(filepath: str, outfile: str) -> None: Loads an AnnData object, performs consensus clustering 
  on the Leiden clustering matrices, and saves the results in the specified output file.

Example usage:
    python script.py dataset.h5ad output.h5ad
"""

import argparse
import os
from tqdm import tqdm
from ccHBGF import ccHBGF
import numpy as np
import scanpy as sc

def main(filepath: str, outfile: str) -> None:
    """
    Perform consensus clustering on Leiden clustering solutions stored in the AnnData object.

    Parameters
    ----------
    filepath : str
        Path to the input .h5ad file containing the AnnData object with Leiden clustering matrices.
    
    outfile : str
        Path to the output .h5ad file where the consensus clustering results will be saved.

    Returns
    -------
    None
    """
   
    # Load the AnnData object from the specified file
    adata = sc.read_h5ad(filepath)
    
    # Identify keys in the `obsm` attribute containing Leiden clustering matrices
    valid_keys = [key for key in adata.obsm_keys() if '_leiden_matrix' in key]
    
    # Perform consensus clustering on each matrix
    for key in tqdm(valid_keys, desc="Resolutions", unit="solution"): 
        rtag = key.split('_')[0]
        
        leiden_matrix = adata.obsm[key]
        
        # Apply consensus clustering using ccHBGF
        adata.obs[f'{rtag}_consensus'] = ccHBGF(leiden_matrix, random_state=0, verbose=False)
        
    # Save the modified AnnData object to the specified output file
    adata.write_h5ad(outfile)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Script to perform consensus clustering on an AnnData object with Leiden solutions.")
    
    # Required arguments
    parser.add_argument("dataset", type=str, help=".h5ad file containing the AnnData object.")
    parser.add_argument("outfile", type=str, help="Output .h5ad file to save the results.")
    
    args = parser.parse_args()
        
    filepath = args.dataset
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"The dataset file {filepath} does not exist.")
    
    outfile = args.outfile
    if not os.path.exists(os.path.dirname(outfile)):
        raise FileNotFoundError(f"The directory for the output file {os.path.dirname(outfile)} does not exist.")
    
    main(filepath, outfile)
