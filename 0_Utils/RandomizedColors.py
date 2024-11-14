"""
[ChatGPT Generated Documentation]
Script for Setting Randomized Colors in AnnData Object

Overview:
    This script defines a function `set_randomized_colors`, which is designed to correct color assignments
    for high cluster numbers in `sc.pl.embedding` plots by adding randomized colors to the
    `adata.uns['{color}_colors']` attribute of an AnnData object.

Key Function:
    - set_randomized_colors: Randomizes and sets a list of colors for the specified categorical attribute.

Example Usage:
    ```python
    import scanpy as sc
    import numpy as np

    # Example data creation
    adata = sc.datasets.pbmc3k()  # Load a sample dataset
    sc.pp.neighbors(adata)
    sc.tl.leiden(adata, resolution=0.5)

    # Randomizing colors for 'leiden' clusters
    set_randomized_colors(adata, 'leiden', palette='tab20')

    # Plot with the new color settings
    sc.pl.umap(adata, color='leiden')
    ```

Requirements:
    - matplotlib
    - numpy
    - scanpy

"""

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
from typing import Any

def set_randomized_colors(
    adata: Any,
    color: str,
    palette: str = 'tab20'
) -> None:
    """
    [ChatGPT Generated Documentation]
    Randomizes and sets colors for a specified categorical attribute in an AnnData object.

    Parameters:
    ----------
    adata : AnnData
        The AnnData object containing the data. This is typically used in single-cell analysis with Scanpy.

    color : str
        The key in `adata.obs` representing the categorical attribute (e.g., 'leiden', 'louvain')
        for which colors need to be randomized.

    palette : str, optional (default='tab20')
        The color palette to be used for generating random colors. This should be a valid Matplotlib colormap name.

    Returns:
    -------
    None
        This function modifies the `adata.uns` attribute in-place by assigning a list of randomized colors
        to `adata.uns['{color}_colors']`.

    Notes:
    ------
    - This function is particularly useful when plotting embeddings with a high number of clusters,
      ensuring that each cluster gets a distinct, randomized color.
    - The function uses the specified Matplotlib color palette and converts the colors to hexadecimal format.

    Example:
    --------
    ```python
    import scanpy as sc

    adata = sc.datasets.pbmc3k()
    set_randomized_colors(adata, 'louvain', palette='tab20')
    sc.pl.umap(adata, color='louvain')
    ```
    """

    # Ensure the color column is treated as a categorical type
    adata.obs[color] = adata.obs[color].astype('category')
    n_cats = adata.obs[color].nunique()

    # Get colors from the specified Matplotlib palette
    colors = plt.get_cmap(palette).colors
    hex_colors = np.array([mcolors.to_hex(color) for color in colors], dtype='object')

    # Randomly assign colors for each category
    adata.uns[f"{color}_colors"] = np.random.choice(hex_colors, n_cats, replace=False)
