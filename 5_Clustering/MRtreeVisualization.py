"""
This script converts a hierarchical label matrix into a Newick tree format and visualizes it using a circular tree layout. 
It applies gradient colors to tree branches based on an external clustering representation, and uses a custom color palette 
for first-level nodes. The tree is visualized using the `ete3` library and styled to highlight node connections.

Overview:
- `labelmat_to_newick`: Converts a label matrix into Newick format, creating a tree structure from the hierarchy.
- `get_gradient_color`: Generates a color gradient between light grey and red based on input values.
- Tree creation and visualization using the `ete3` library with a circular layout and customizable node styles.

Dependencies:
- numpy: For matrix operations and hierarchical label extraction.
- ete3: For tree structure creation and visualization.
- matplotlib: For color interpolation and conversion to hex format.

Functions:
- labelmat_to_newick(labelmat: np.ndarray, add_root: bool = False, wrap_quotes: bool = True) -> str: 
  Converts a label matrix into Newick tree format.
- get_gradient_color(value: float) -> str: Generates a color gradient between grey and red based on the input value.

Example usage:
    newick_tree = labelmat_to_newick(labelmat_matrix, add_root=True)
    t = Tree(newick_tree, format=8)
    t.render("%%inline", tree_style=ts, h=900, w=900)
"""

import numpy as np
from ete3 import Tree, TreeStyle, NodeStyle
import matplotlib

def labelmat_to_newick(labelmat: np.ndarray, add_root: bool = False, wrap_quotes: bool = True) -> str:
    """
    Converts a label matrix to Newick format. The label matrix contains hierarchical labels, and 
    this function constructs a Newick tree by iterating through the matrix levels.

    Parameters
    ----------
    labelmat : np.ndarray
        The hierarchical label matrix with each column representing a hierarchical level.
    
    add_root : bool, optional
        Whether to add a root to the Newick string. Default is False.
    
    wrap_quotes : bool, optional
        Whether to wrap labels in quotes. Default is True.

    Returns
    -------
    str
        Newick formatted string representing the hierarchy as a tree.
    """
    first_layer = labelmat[:, 0]
    first_unique = np.unique(first_layer)
    newick_format = f'{",".join(first_unique)}'
    
    for lvl in range(labelmat.shape[1] - 1):
        layer = labelmat[:, lvl]
        next_layer = labelmat[:, lvl + 1]
        
        unique = np.unique(layer)
        
        for val in unique:
            unique_mask = layer == val
            next_unique = np.unique(next_layer[unique_mask])
            
            updated_val = f'({",".join(next_unique)}){val}'        
            newick_format = newick_format.replace(val, updated_val)
    
    if add_root:
        return f'({newick_format});'
    else:
        return f'{newick_format};'

# Define a color palette for first-level nodes
color_palette = ['LightSalmon', 'LightPink', 'Plum', 'MediumAquamarine', 'Cornsilk', 'Cyan']

# Convert label matrix to Newick format and load the tree
newick_tree = labelmat_to_newick(labelmat, add_root=True)
t = Tree(newick_tree, format=8)

# Define custom TreeStyle for circular layout
ts = TreeStyle()
ts.show_leaf_name = False
ts.show_branch_length = False
ts.show_scale = False
ts.mode = "c"  # Circular layout
ts.root_opening_factor = 0.05

# Initialize a counter for cycling through the color palette
color_counter = 0

# Traverse the tree to customize nodes and apply styles
for node in t.traverse():

    # Define the default node style
    style = NodeStyle()
    style["shape"] = "circle"    # Node shape
    style["size"] = 0            # Hide the node circle by setting size to 0
    style["fgcolor"] = "white"   # Node foreground color

    # Apply the style to the node
    node.set_style(style)

# Render the tree inline with customized dimensions
t.render("%%inline", tree_style=ts, h=900, w=900)
