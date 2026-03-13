###############################################################################
###############################################################################

# Program: Github_Reduced.py
# Author: Sergio Cámara Peña
# Date: 12/03/2026
# Version: FINAL

###############################################################################
###############################################################################

# %% Import libraries
import os
import scanpy as sc
import random

# %% Read files
os.chdir("/home/scamara/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Raw_Atlas")
adata_V4 = sc.read_h5ad("Python_scVI_adata_big_V4_state4_Normalized.h5ad")

# %% Remove raw and genes that have few counts
adata_V4.raw = None
sc.pp.filter_genes(adata_V4, min_counts=20)
sc.pp.filter_genes(adata_V4, min_cells=20)

# %% Remove extra columns
cols_to_remove = ["Marked_cells", "barcode", "timepoint", "barcode_timepoint", "celltypist_predicted_labels", "celltypist_over_clustering", "celltypist_majority_voting", "celltypist_conf_score", "Product", "_scvi_batch", "_scvi_labels"]
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

for key in ["distances", "connectivities"]:
    if key in adata_V4.obsp:
        del adata_V4.obsp[key]

# %% Random subset of 1000 cells
random.seed(2504) # Set seed for reproducibility

n_cells = 1000
idx = random.sample(range(adata_V4.n_obs), n_cells)

adata_mini = adata_V4[idx, :].copy()

# %% Save mini object
adata_mini.write("Atlas_DEMO.h5ad")

# %% End of script
