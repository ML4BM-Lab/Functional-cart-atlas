###############################################################################
###############################################################################

# Program: 5.1.1_Change_metadata_for_shiny.py
# Author: Sergio Cámara Peña
# Date: 12/05/2025
# Version: FINAL

###############################################################################
###############################################################################

# %% Import libraries
import os
import scanpy as sc

# %% Read files
os.chdir("/home/scamara/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Raw_Atlas")
adata_V4 = sc.read_h5ad("Python_scVI_adata_big_V4_state4_Normalized.h5ad")

# %% Change metadata columns
adata_V4.obs.rename(columns={'ScFv': 'CAR_Type'}, inplace=True)

# %% Remove raw and genes that have few counts
adata_V4.raw = None
sc.pp.filter_genes(adata_V4, min_counts=20)
sc.pp.filter_genes(adata_V4, min_cells=20)

# %% Remove extra columns
cols_to_remove = ["Marked_cells", "barcode", "timepoint", "barcode_timepoint", "manual_celltype_annotation_high_V3", "manual_celltype_annotation_low", "pct_counts_in_top_20_genes", "total_counts_ribo", "log1p_total_counts_ribo", "pct_counts_ribo", "total_counts_hb", "log1p_total_counts_hb", "pct_counts_hb", "score_ery"]
adata_V4.obs.drop(columns=cols_to_remove, inplace=True)

keep_uns = ["log1p", "umap"]

for key in list(adata_V4.uns.keys()):
    if key not in keep_uns:
        del adata_V4.uns[key]

for key in ["_scvi_extra_categorical_covs"]:
    if key in adata_V4.obsm:
        del adata_V4.obsm[key]

for key in ["X_pca"]:
    if key in adata_V4.obsm:
        del adata_V4.obsm[key]

adata_V4.var.drop(columns=["n_counts", "n_cells"], inplace=True)

for key in ["distances"]:
    if key in adata_V4.obsp:
        del adata_V4.obsp[key]

# %% Save files
os.chdir("/home/scamara/data/scamara/Atlas_Mieloma_Multiple/Shiny_app")
adata_V4.write("Python_scVI_adata_big_V4_state4_Normalized_for_shiny.h5ad")

# %% End of script