"""
This script performs Leiden clustering on an AnnData object for multiple resolutions and seeds.

Overview:
- `main`: Main function that loads an AnnData object, performs Leiden clustering at multiple resolutions 
  and random seeds, and stores the results in the AnnData object's `obsm` attribute.
- `format_number`: Helper function to format a numeric string into an integer or float.

Dependencies:
- argparse: For parsing command-line arguments.
- os: For file and directory checking.
- tqdm: For displaying progress bars.
- numpy: For numerical operations.
- scanpy: For processing and analyzing single-cell RNA-seq data.
- rapids_singlecell: For GPU-accelerated computation of neighbors and Leiden clustering.
- typing: For type hinting.

Functions:
- main(filepath: str, resolutions: List[Union[int, float]], outfile: str, n_seed: int) -> None: 
  Loads an AnnData object from a file, performs Leiden clustering for each resolution and seed, 
  and saves the results to an output file.
- format_number(num: str) -> Union[int, float]: Formats a string as an integer or float.

Example usage:
    python script.py dataset.h5ad output.h5ad 0.1 0.5 1.0 -n 10
"""

import argparse
import os
from tqdm import tqdm
import numpy as np
import scanpy as sc
import rapids_singlecell as rsc
from typing import List, Union

def main(filepath: str, resolutions: List[Union[int, float]], outfile: str, n_seed: int) -> None:
    """
    Perform Leiden clustering on an AnnData object for multiple resolutions and seeds.

    Parameters
    ----------
    filepath : str
        Path to the input .h5ad file containing the AnnData object.
    
    resolutions : List[Union[int, float]]
        List of Leiden clustering resolutions to apply.
    
    outfile : str
        Path to the output .h5ad file where the results will be saved.
    
    n_seed : int
        Number of random seeds to use for the Leiden clustering.

    Returns
    -------
    None
    """
    
    # Load the AnnData object from the specified file
    adata = sc.read_h5ad(filepath)
    
    # Preprocess: Compute the neighborhood graph using the specified embedding
    rsc.pp.neighbors(adata, n_neighbors=15, use_rep='X_embed', metric='euclidean')
    
    # Iterate over each resolution
    for r in tqdm(resolutions, desc="Resolutions", unit="res"):
        matrix = np.zeros(shape=(adata.shape[0], n_seed))
        
        # Inner loop for seeds
        for seed in tqdm(range(n_seed), desc=f"Seeds for r={r}", unit="seed", leave=False):
            rsc.tl.leiden(adata, resolution=r, random_state=seed)
            matrix[:, seed] = adata.obs['leiden'].values
            del adata.obs['leiden']
        
        # Store the resulting matrix in the AnnData object
        adata.obsm[f'r{r}_leiden_matrix'] = matrix
    
    # Save the modified AnnData object to the specified output file
    adata.write_h5ad(outfile)

def format_number(num: str) -> Union[int, float]:
    """
    Converts a string representing a number to an integer or float.

    Parameters
    ----------
    num : str
        String to be converted to a number.
    
    Returns
    -------
    Union[int, float]
        The input string converted to an integer if it's a whole number, otherwise a float.
    """
    num = float(num)
    if num.is_integer():
        return int(num)
    else:
        return num

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Script to perform Leiden clustering on an AnnData object with specified resolutions.")
    
    # Required arguments
    parser.add_argument("dataset", type=str, help=".h5ad file containing the AnnData object.")
    parser.add_argument("outfile", type=str, help="Output .h5ad file to save the results.")
    parser.add_argument('resolutions', nargs='+', type=format_number, help="List of Leiden resolutions to apply.")
    
    # Optional argument
    parser.add_argument('-n', '--n_seed', type=int, default=10, help="Number of random seeds for the Leiden clustering (default: 5).")
    
    args = parser.parse_args()
    resolutions = args.resolutions
    
    filepath = args.dataset
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"The dataset file {filepath} does not exist.")
    
    outfile = args.outfile
    if not os.path.exists(os.path.dirname(outfile)):
        raise FileNotFoundError(f"The directory for the output file {os.path.dirname(outfile)} does not exist.")
    
    n_seed = args.n_seed
    
    main(filepath, resolutions, outfile, n_seed)
