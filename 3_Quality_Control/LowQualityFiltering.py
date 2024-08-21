"""
This script provides a set of functions for preprocessing single-cell RNA-seq data using Scanpy.
It includes functions to calculate quality control (QC) metrics, normalize the data, score marker genes, and filter cells based on specific thresholds.

Overview:
- `compute_QC_metrics`: Computes and stores quality control metrics in the AnnData object.
- `normalize`: Normalizes the AnnData object and applies logarithmic transformation.
- `compute_marker_gene_scores`: Scores marker genes for specific cell types and optionally displays UMAP plots.
- `filter_cells`: Filters the AnnData object based on minimum gene counts, cell counts, and mitochondrial content.

Dependencies:
- json: For loading marker gene lists from a file.
- scanpy: For processing and analyzing single-cell RNA-seq data in AnnData format.

Functions:
- compute_QC_metrics(adata: sc.AnnData) -> None: Adds QC metrics to the AnnData object.
- normalize(adata: sc.AnnData) -> None: Normalizes and log-transforms the AnnData object.
- compute_marker_gene_scores(adata: sc.AnnData, fp: str, keys: Optional[List[str]] = None, show: bool = False) -> None: Scores marker genes and optionally plots UMAP.
- filter_cells(adata: sc.AnnData, min_genes: int, min_counts: int, pct_mito: float) -> sc.AnnData: Filters cells based on provided thresholds.

Example usage:
    adata = sc.read_h5ad('sample_data.h5ad')
    compute_QC_metrics(adata)
    normalize(adata)
    compute_marker_gene_scores(adata, 'marker_genes.json')
    adata_filtered = filter_cells(adata, min_genes=450, min_counts=2000, pct_mito=5)
"""

import json
import scanpy as sc
from typing import List, Optional

def compute_QC_metrics(adata: sc.AnnData) -> None:
    """
    Computes and stores quality control metrics in the AnnData object.

    Args:
        adata (sc.AnnData): The AnnData object containing single-cell RNA-seq data.

    Returns:
        None: The function modifies the AnnData object in place.
    """
    qc_metrics = sc.pp.calculate_qc_metrics(adata, percent_top=False, log1p=False, return_dict=True)
    adata.obs['total_counts'] = qc_metrics['total_counts']
    adata.var['mt'] = adata.var_names.str.startswith('mt-')
    adata.obs['pct_mito'] = qc_metrics['pct_counts_mt']

def normalize(adata: sc.AnnData) -> None:
    """
    Normalizes the AnnData object to a total count per cell and applies logarithmic transformation.

    Args:
        adata (sc.AnnData): The AnnData object to be normalized.

    Returns:
        None: The function modifies the AnnData object in place.
    """
    sc.pp.normalize_total(adata, target_sum=1e6)
    sc.pp.log1p(adata)

def compute_marker_gene_scores(adata: sc.AnnData, fp: str, keys: Optional[List[str]] = None, show: bool = False) -> None:
    """
    Scores marker genes for specific cell types and optionally displays UMAP plots.

    Args:
        adata (sc.AnnData): The AnnData object containing single-cell RNA-seq data.
        fp (str): File path to a JSON file containing marker genes.
        keys (Optional[List[str]]): Specific cell types to score. If None, scores all types in the JSON.
        show (bool): Whether to display UMAP plots of the scores.

    Returns:
        None: The function modifies the AnnData object in place and optionally displays plots.
    """
    with open(fp, 'r') as file:
        cell_types = json.load(file)
        
    valid_keys = cell_types.keys()

    if keys:
        valid_keys = [k for k in valid_keys if k in keys]
    
    for key in valid_keys:
        score_name = key + "_score"
        sc.tl.score_genes(adata, gene_list=cell_types[key], score_name=score_name)
    
    if show:
        vmin = [-0.1 for _ in valid_keys]
        figure = sc.pl.umap(adata, color=valid_keys, ncols=3, legend_loc=None, return_fig=True, show=False, size=20, vmin=vmin, frameon=False)
        figure.tight_layout()
        figure.show()

def filter_cells(adata: sc.AnnData, min_genes: int, min_counts: int, pct_mito: float) -> sc.AnnData:
    """
    Filters cells based on minimum gene counts, cell counts, and mitochondrial content.

    Args:
        adata (sc.AnnData): The AnnData object to be filtered.
        min_genes (int): Minimum number of genes expressed per cell.
        min_counts (int): Minimum number of counts per cell.
        pct_mito (float): Maximum allowed percentage of mitochondrial counts.

    Returns:
        sc.AnnData: The filtered AnnData object.
    """
    sc.pp.filter_cells(adata, min_counts=min_counts)
    adata.uns['min_count_cutoff'] = min_counts

    sc.pp.filter_cells(adata, min_genes=min_genes)
    adata.uns['min_gene_cutoff'] = min_genes

    adata = adata[adata.obs.pct_mito < pct_mito]
    adata.uns['min_pct_mito_cutoff'] = pct_mito

    return adata
