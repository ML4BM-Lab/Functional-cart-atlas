###############################################################################
###############################################################################

# Program: Change_metadata_for_UCSC.py
# Author: Sergio Cámara Peña
# Date: 27/01/2026
# Version: FINAL

###############################################################################
###############################################################################

# %% Import libraries
import os
import scanpy as sc

# %% Read files
os.chdir("/home/scamara/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Raw_Atlas")
adata_V5 = sc.read_h5ad("Atlas_integ_scArches_FINAL_V5_Normalized.h5ad")
adata_V5_orig = adata_V5.copy()

# %% Change metadata columns
adata_V5.obs.rename(columns={'ScFv': 'CAR_Type'}, inplace=True)

# %% Remove raw and genes that have few counts
adata_V5.raw = None
sc.pp.filter_genes(adata_V5, min_counts=20)
sc.pp.filter_genes(adata_V5, min_cells=20)

# %% Remove extra columns
cols_to_remove = ["Marked_cells", "barcode", "timepoint", "barcode_timepoint", "celltypist_predicted_labels", "celltypist_over_clustering", "celltypist_majority_voting", "celltypist_conf_score", "Product", "_scvi_batch", "_scvi_labels", "leiden"]
adata_V5.obs.drop(columns=cols_to_remove, inplace=True)

keep_uns = ["log1p", "umap"]

for key in list(adata_V5.uns.keys()):
    if key not in keep_uns:
        del adata_V5.uns[key]

for key in ["_scvi_extra_categorical_covs"]:
    if key in adata_V5.obsm:
        del adata_V5.obsm[key]

for key in ["X_pca"]:
    if key in adata_V5.obsm:
        del adata_V5.obsm[key]

adata_V5.var.drop(columns=["n_counts", "n_cells"], inplace=True)

for key in ["distances", "connectivities"]:
    if key in adata_V5.obsp:
        del adata_V5.obsp[key]

# %% Save files
os.chdir("/home/scamara/data/scamara/Atlas_Mieloma_Multiple/UCSC_Browser")
adata_V5.write("Atlas_integ_scArches_FINAL_V5_Normalized_for_UCSC.h5ad")

# %% End of script