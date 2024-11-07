"""[ChatGPT Generated Documentation]
This script provides a function to visualize spatial data from an AnnData object, allowing for slice selection, annotation overlay, 
and optional zooming into specific annotated areas. It supports flexible plotting with options for color selection, 
plot customization, and figure control.

Functions included:
1. `plot_spatial`: Plots spatial annotations from an AnnData object with flexible slice, selection, and zoom options.

Example Usage:
--------------
# Plot spatial data with selected annotation IDs and zoom
plot_spatial(adata, label="cluster_label", annotation_res=25, annotation_matrix=annotation_volume, selected_ids={1, 2, 3}, zoom_selection=True)
"""

import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from typing import Optional, Dict, Union

def plot_spatial(
    adata: pd.DataFrame,
    label: str,
    annotation_res: int,
    annotation_matrix: np.ndarray,
    slice_idx: Optional[int] = None,
    selected_ids: Union[set, Dict[int, str]] = {},
    return_fig: bool = False,
    show: bool = True,
    zoom_selection: bool = False,
    **kwargs
) -> Optional[plt.Figure]:
    """[ChatGPT Generated Documentation]
    Plots spatial data with annotations from an AnnData object, allowing for slice selection, region highlighting, 
    and optional zooming on selected annotations.

    Parameters:
    - adata (pd.DataFrame): The AnnData object containing spatial information.
    - label (str): Column name in `adata.obs` used for coloring the points.
    - annotation_res (int): Resolution of the annotation matrix in micrometers per voxel.
    - annotation_matrix (np.ndarray): 3D annotation matrix for spatial labeling.
    - slice_idx (Optional[int]): Specific index for the slice along the x-axis. If None, it selects the slice with 
      the highest sum of the label values. Defaults to None.
    - selected_ids (Union[set, Dict[int, str]]): IDs to be highlighted in the annotation matrix. Defaults to an empty set.
    - return_fig (bool): If True, returns the matplotlib Figure object. Defaults to False.
    - show (bool): If True, displays the plot. If False, it suppresses the display. Defaults to True.
    - zoom_selection (bool): If True, zooms into the area around the selected IDs. Defaults to False.
    - **kwargs: Additional keyword arguments passed to `ax.scatter` for customizing the point plot.

    Returns:
    - Optional[plt.Figure]: The matplotlib Figure if `return_fig` is True; otherwise, None.

    Notes:
    - This function modifies spatial coordinates based on `annotation_res` to match the annotation matrix.
    - If `zoom_selection` is enabled, the plot limits adjust to show the selected area with a margin.
    """
    # Prepare voxel data based on spatial coordinates and resolution
    voxel_df = adata.obs[[label, 'x_ccf', 'y_ccf', 'z_ccf']]
    voxel_df.loc[:, 'x_ccf'] = (voxel_df['x_ccf'] * 1000 // annotation_res).astype(int)
    voxel_df.loc[:, 'y_ccf'] = (voxel_df['y_ccf'] * 1000 // annotation_res).astype(int)
    voxel_df.loc[:, 'z_ccf'] = (voxel_df['z_ccf'] * 1000 // annotation_res).astype(int)

    # Automatically select the slice with the highest sum of the label values if `slice_idx` is not provided
    if slice_idx is None:
        grouped_df = voxel_df.groupby('x_ccf')
        slice_idx = int(grouped_df[label].sum().idxmax())
    
    # Extract the selected slice from the annotation matrix
    slice = annotation_matrix[slice_idx]
    
    # Create a mask for areas with no annotations in the selected slice
    empty_mask = slice == 0
    
    # Create selection and unselected masks based on `selected_ids`
    if selected_ids:
        selected_ids = {int(id) for id in selected_ids}
        selection_mask = np.vectorize(selected_ids.__contains__)(slice)
        unselected_mask = ~(selection_mask | empty_mask)
    else:
        unselected_mask = ~empty_mask

    # Initialize the image with a white background for empty areas
    image = np.zeros((*slice.shape, 3), dtype=np.uint8)
    image[empty_mask] = [255, 255, 255]
    image[unselected_mask] = [211, 211, 211]

    # Apply colors to selected IDs if any are provided
    if selected_ids:
        unique_ids = np.unique(slice[selection_mask])
        color_map = plt.cm.Dark2.colors  # Example cyclic colormap
        color_map = [[int(255 * v) for v in colors.to_rgb(color)] for color in color_map]

        for i, id in enumerate(unique_ids):
            color_idx = i % len(color_map)
            image[slice == id] = color_map[color_idx]

    # Create the plot
    fig, ax = plt.subplots()

    # Scatter plot of the voxel data with specified color mapping
    ax.scatter(
        voxel_df['z_ccf'], voxel_df['y_ccf'], s=1, c=adata.obs[label], cmap='viridis', alpha=1, zorder=0, **kwargs
    )
    ax.imshow(image, alpha=0.4, zorder=1)
    
    # Zoom into the selection if `zoom_selection` is enabled and there are selected points
    if selected_ids and zoom_selection:
        bb = np.argwhere(selection_mask)
        if bb.size > 0:
            (y_min, z_min), (y_max, z_max) = bb.min(0), bb.max(0)
            height, width = y_max - y_min, z_max - z_min
            margin_y, margin_z = int(height * 0.2), int(width * 0.2)
            ax.set_ylim(y_max + margin_y, y_min - margin_y)
            ax.set_xlim(z_max + margin_z, z_min - margin_z)
    
    # Turn off axis for a clean plot
    ax.axis('off')
    
    # Display, return, or close the plot based on user options
    if return_fig:
        return fig
    elif show:
        plt.show()
    else:
        plt.close(fig)
