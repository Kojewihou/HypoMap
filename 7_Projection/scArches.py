"""
This function embeds a query AnnData object using a pre-trained scVI model. It prepares the query 
data, trains a query model, and stores the latent embedding in the AnnData object's `obsm` attribute.

Overview:
- `embed_query`: Prepares query data, trains a scVI model, and obtains the latent representation (embedding) 
  of the query dataset.

Dependencies:
- scvi: For loading and training the scVI model on query data.
- scanpy: For handling AnnData objects and pre-processing data.

Functions:
- embed_query(query: sc.AnnData, model_path: str, embed: str = 'X_embed') -> None: Prepares the query dataset, 
  loads and trains a query scVI model, and stores the resulting latent embeddings in `obsm`.

Example usage:
    query = sc.read_h5ad('query_data.h5ad')
    embed_query(query, 'path/to/model', embed='X_embed')
"""

import scanpy as sc
import scvi

def embed_query(query: sc.AnnData, model_path: str, embed: str = 'X_embed') -> None:
    """
    Embeds a query AnnData object using a pre-trained scVI model. The function prepares the query data, 
    loads the model, trains a query-specific scVI model, and stores the latent embedding in the `obsm` 
    attribute of the query AnnData object.

    Parameters
    ----------
    query : sc.AnnData
        The query AnnData object containing single-cell RNA-seq data to be embedded. The object is assumed 
        to have the necessary preprocessing done, including a `brain_section_label` in `query.obs`.
    
    model_path : str
        Path to the pre-trained scVI model to be used for embedding the query data.
    
    embed : str, optional
        The name of the key in `query.obsm` where the latent embedding will be stored. Default is 'X_embed'.

    Returns
    -------
    None
    """

    # Prepare the query metadata for scVI query model setup
    query.obs['Sample'] = query.obs['brain_section_label']  # Copy brain section label to 'Sample'
    query.obs['Suspension_bin'] = 0  # Add 'Suspension_bin' to metadata
    query.obs['total_counts'] = query.X.sum(axis=1).A1  # Calculate total counts and store in 'total_counts'

    # Prepare the query data for use with the pre-trained scVI model
    scvi.model.SCVI.prepare_query_anndata(query, model_path)

    # Load the pre-trained scVI model with the query data
    vae_q = scvi.model.SCVI.load_query_data(
        query,
        model_path,
    )

    # Train the query-specific scVI model
    vae_q.train(max_epochs=25, batch_size=1024, plan_kwargs={'weight_decay': 0.0}, accelerator='gpu')

    # Get the latent representation of the query data and store it in the obsm attribute
    query.obsm[embed] = vae_q.get_latent_representation()

