import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

import numpy as np
import scanpy as sc
import seaborn as sns

from anndata import AnnData
from typing import Any, Optional, Literal
from scipy.optimize import linear_sum_assignment

def plot_clusters_overlap(original_labels, new_labels):
    original_labels = original_labels.astype(str)
    new_labels = new_labels.astype(str)

    ctmx = sc.metrics.confusion_matrix(original_labels, new_labels)

    row_idx, col_idx = linear_sum_assignment(-ctmx) #-cm to maximize, o/w it minimizes 

    n_merges = [len(x.split('_')) for x in ctmx.columns]

    # Create the figure and axes objects
    fig, ax = plt.subplots(figsize=(10, 8))

    # Plot the heatmap
    sns.heatmap(ctmx.iloc[row_idx, col_idx], annot=False, cmap='viridis', ax=ax)

    # Create a second axis for the bar plot (to place it below the heatmap)
    ax_bar = ax.inset_axes([0, 1.01, 1, 0.05])  # [left, bottom, width, height]

    # Plot the bar for n_merges below the heatmap
    ax_bar.bar(np.arange(len(n_merges)), n_merges, color='red')

    # Hide the x and y axes of the bar plot for a cleaner look
    ax_bar.axis('off')

    ax.get_xaxis().set_ticklabels([])

    # Show the plot
    plt.show()

def plot_adata(
    adata: AnnData, 
    plot_type: Literal['tsne', 'umap'] = 'tsne', 
    **kwargs
) -> None:
    if plot_type not in {'tsne', 'umap'}:
        raise ValueError("Invalid plot_type. Choose 'tsne' or 'umap'.")

    # Get palette from kwargs or use default
    palette: str = kwargs.get('palette', 'tab20')

    # Normalize 'color' argument to a list
    if 'color' in kwargs:
        if isinstance(kwargs['color'], str):
            kwargs['color'] = [kwargs['color']]
        
        # Randomize colors for specified annotations
        for color in kwargs['color']:
            if f"{color}_colors" in adata.uns:
                del adata.uns[f"{color}_colors"]
            _set_randomized_colors(adata, color, palette)

    # Remove 'palette' from kwargs to avoid conflicts
    kwargs.pop('palette', None)

    # Call the appropriate plotting function
    if plot_type == 'tsne':
        sc.pl.tsne(adata, **kwargs)
    elif plot_type == 'umap':
        sc.pl.umap(adata, **kwargs)
        
def _set_randomized_colors(
    adata: Any,
    color: str,
    palette: str = 'tab20'
) -> None:

    # Ensure the color column is treated as a categorical type
    n_cats = adata.obs[color].astype('category').cat.categories.__len__()

    # Get colors from the specified Matplotlib palette
    colors = plt.get_cmap(palette).colors
    hex_colors = np.array([mcolors.to_hex(color) for color in colors], dtype='object')

    if n_cats <= len(hex_colors):
        adata.uns[f"{color}_colors"] = np.array(hex_colors[:n_cats])
    else:
        adata.uns[f"{color}_colors"] = np.random.choice(hex_colors, n_cats, replace=True)

