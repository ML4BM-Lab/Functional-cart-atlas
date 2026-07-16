###############################################################################
###############################################################################

# Program: 2.2_Celltypist_annotation.py
# Author: Sergio Cámara Peña
# Date: 04/09/2023
# Version: V4

###############################################################################
###############################################################################

# %% Import needed packages + Download and select celltypist model
import celltypist
from celltypist import models
import scanpy as sc
import os

Online=True
if Online:
    models.download_models(force_update=True)
model_to_use = models.Model.load(model="Immune_All_Low.pkl")


####################################################################################################################
####################################################################################################################
####################################################################################################################

# %% Select (.csv) file on where to do the predictions
dir_merged = "/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/Python-Celltypist/V4/Seurat_merged.h5ad"
adata_Seurat_merged = sc.read(dir_merged)
sc.pp.normalize_per_cell(adata_Seurat_merged, counts_per_cell_after=10000)
sc.pp.log1p(adata_Seurat_merged)

# %% Predictions generation
predictions_merged = celltypist.annotate(
    filename=adata_Seurat_merged,
    model=model_to_use,
    majority_voting=True,
    transpose_input=True,
)
print(predictions_merged.predicted_labels)

# %% Choose destination path
adata = predictions_merged.to_adata()
os.chdir(
    "/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/Python-Celltypist/V4"
)

# Save metadata table to a file
adata.obs.to_csv("celltypist_metadata_table.csv", sep=",", index=True, header=True)

# %% End of program
