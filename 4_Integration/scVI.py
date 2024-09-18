import argparse
import json
import logging
import os
from shutil import copy2, SameFileError
import warnings

logger = logging.getLogger(__name__)
warnings.filterwarnings('ignore', category=DeprecationWarning)
warnings.filterwarnings("ignore", category=UserWarning)

import numpy as np
import pandas as pd
import scanpy as sc
import scvi
from torch import set_float32_matmul_precision

scvi.settings.seed=0
set_float32_matmul_precision('medium')

global PARAMETERS
global OUTDIR
EMBEDDING_KEY = 'X_embed'

def main(args: argparse.Namespace) -> None:
    
    global PARAMETERS
    global OUTDIR 
    
    with open(args.parameter_file) as parameter_file:
        PARAMETERS = json.load(parameter_file)
    
    OUTDIR = os.path.join(PARAMETERS['outdir'], PARAMETERS['model_name'])

    logging.info(f'Initializing {OUTDIR}')
    
    initialize_model_directory(args.parameter_file)
    
    adata = sc.read(PARAMETERS['dataset'])
    logger.info(f'Successfully loaded raw dataset.')
    
    hvgs = pd.read_csv(PARAMETERS['hvgs'], index_col=0, header=None)
    logger.info(f'Subsetting to highly variable genes ({len(hvgs)}).')
    adata = adata[:, hvgs.index].copy()

    logger.info(adata)

    logger.info('Getting model embedding.')
    adata.obsm[EMBEDDING_KEY] = get_model_embedding(adata, PARAMETERS, get_scvi_embedding)

    del adata.X

    logger.info('Calculating neighbors.')
    neighbors_parameters = PARAMETERS.get('neighbors', {})
    sc_neighbors_kwargs = {"use_rep": EMBEDDING_KEY} | neighbors_parameters.get('kwargs', {})
    rsc_neighbors_kwargs = {"use_rep": EMBEDDING_KEY} | neighbors_parameters.get('kwargs', {})
    if neighbors_parameters.get('gpu', False):
        try:
            import rapids_singlecell as rsc
            rsc.pp.neighbors(adata, **rsc_neighbors_kwargs)
        except ModuleNotFoundError as e:
            logger.warning(f'RAPIDS gpu accelration unavailable, falling back to cpu. \n\t!> {e}')
            sc.pp.neighbors(adata, **sc_neighbors_kwargs)   
    else:
        sc.pp.neighbors(adata, **sc_neighbors_kwargs)
        
    logger.info('Calculating UMAP.')
    umap_parameters = PARAMETERS.get('umap', {})
    sc_umap_kwargs = {} | umap_parameters.get('kwargs', {})
    rsc_umap_kwargs = {} | umap_parameters.get('kwargs', {})
    if umap_parameters.get('gpu', False):
        try:
            import rapids_singlecell as rsc
            rsc.tl.umap(adata, **rsc_umap_kwargs)
        except ModuleNotFoundError as e:
            logger.warning(f'RAPIDS gpu accelration unavailable, falling back to cpu. \n\t!> {e}')
            sc.tl.umap(adata, **sc_umap_kwargs)
    else:
        sc.tl.umap(adata, **sc_umap_kwargs)
        
    logger.info('Calculating t-SNE.')
    tsne_parameters = PARAMETERS.get('tsne', {})
    sc_tsne_kwargs = {"perplexity": 10} | tsne_parameters.get('kwargs', {})
    rsc_tsne_kwargs = {"perplexity": 10} | tsne_parameters.get('kwargs', {})
    if tsne_parameters.get('gpu', False):
        try:
            import rapids_singlecell as rsc
            rsc.tl.tsne(adata, use_rep=EMBEDDING_KEY, **rsc_tsne_kwargs)
        except ModuleNotFoundError as e:
            logger.warning(f'RAPIDS gpu accelration unavailable, falling back to cpu. \n\t! {e} !')
            sc.tl.tsne(adata, use_rep=EMBEDDING_KEY, **sc_tsne_kwargs)
    else:
        sc.tl.tsne(adata, use_rep=EMBEDDING_KEY, **sc_tsne_kwargs)
        
    logger.info('Calculating leiden.')
    leiden_parameters = PARAMETERS.get('leiden', {})
    sc_leiden_kwargs = {"flavor": "igraph", "n_iterations": 2} | leiden_parameters.get('kwargs', {})
    rsc_leiden_kwargs = {} | leiden_parameters.get('kwargs', {})
    if leiden_parameters.get('gpu', False):
        try:
            import rapids_singlecell as rsc
            rsc.tl.leiden(adata, **rsc_leiden_kwargs)
        except ModuleNotFoundError as e:
            print(f'RAPIDS gpu accelration unavailable, falling back to cpu. \n\t! {e} !')
            sc.tl.leiden(adata, **sc_leiden_kwargs)
    else:
        sc.tl.leiden(adata, **sc_leiden_kwargs)
    
    embedding_only_file = os.path.join(OUTDIR, 'datasets/embed_only.h5ad')
    adata.write_h5ad(embedding_only_file)
    logger.info(f'Saving embedding to {embedding_only_file}')

def initialize_model_directory(parameter_file: str) -> None:
    """Create directory for new model (WORKDIR), hard linking in datasets and copying the parameter file.

    Args:
        parameter_file (str): Filepath to parameter file.
    """
    global PARAMETERS
    global OUTDIR
    
    os.makedirs(os.path.join(OUTDIR, 'datasets'), exist_ok=True)
    os.makedirs(os.path.join(OUTDIR, 'metrics'), exist_ok=True)
    os.makedirs(os.path.join(OUTDIR, 'models'), exist_ok=True)
    os.makedirs(os.path.join(OUTDIR, 'plots'), exist_ok=True)
    logger.info(f'Subdirectories created in {OUTDIR}')

    try:
        copy2(PARAMETERS['hvgs'],
                os.path.join(OUTDIR, 'metrics/hvgs.csv'))
        copy2(PARAMETERS['dataset'],
                os.path.join(OUTDIR, 'datasets/raw.h5ad'))
        logger.info(f'Copied datasets into {OUTDIR}')
    except SameFileError:
        logger.warn(f'Datasets already present in {OUTDIR}')
        pass
    
    try:
        copy2(parameter_file, os.path.join(OUTDIR,'parameters.json'))
    except SameFileError:
        logger.warn(f'Parameter file is already present in {OUTDIR}')
        pass

def get_model_embedding(adata: sc.AnnData, parameters: dict, model_func) -> np.array:
    return model_func(adata, parameters)

def get_scvi_embedding(adata: sc.AnnData, parameters: dict) -> np.array:  
    global OUTDIR
    
    model_dir = os.path.join(OUTDIR, 'models/scvi')
    model_parameters = parameters['model']
    logger.info(model_parameters)
    
    if os.path.exists(model_dir):
        model = scvi.model.SCVI.load(model_dir, adata=adata)
    else:
        scvi.model.SCVI.setup_anndata(adata, **model_parameters.get('setup_anndata_kwargs', {}))

        model = scvi.model.SCVI(adata, **model_parameters.get('model_kwargs', {}), )

        model.train(**model_parameters.get('train_kwargs', {}))

        model.save(model_dir)
    
    fig = plot_model_history(model.history)
    fig.savefig(os.path.join(OUTDIR, 'plots/model_history.png'))

    return model.get_latent_representation()

def plot_model_history(history):
    import matplotlib.pyplot as plt
    
    fig = plt.figure(figsize=(10,20))
    gs = fig.add_gridspec(3, hspace=0)
    axs = gs.subplots(sharex=True)

    axs[0].plot(history['elbo_train'], label='Train', color='black')
    axs[0].plot(history['elbo_validation'], label='Validation', color='red')
    axs[0].set_ylabel('ELBO')
    axs[0].legend()

    axs[1].plot(history['reconstruction_loss_train'], label='Train', color='black')
    axs[1].plot(history['reconstruction_loss_validation'], label='Validation', color='red')
    axs[1].set_ylabel('Reconstruction Loss')
    axs[1].legend()

    axs[2].plot(history['kl_weight'], label='Training KL Weight', color='purple')
    axs[2].set_ylabel('KL Weight')
    axs[2].set_xlabel('Epoch')
    axs[2].xaxis.set_major_locator(plt.MaxNLocator(integer=True))
    axs[2].legend()

    for ax in axs:
        ax.label_outer()
        
    plt.close()

    return fig

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Setup raw data for integration model setup')

    parser.add_argument('parameter_file', metavar='[ Integration Parameters (JSON) ]', help='Integration Parameters')
    parser.add_argument('--verbose', '-v', action='count', default=0, help='Increase output verbosity (e.g., -vv for more)')
    
    args = parser.parse_args()

    log_level={0: logging.WARNING, 1: logging.INFO, 2: logging.DEBUG}
    logging.basicConfig(level=log_level.get(args.verbose, logging.DEBUG), format='%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

    main(args)
