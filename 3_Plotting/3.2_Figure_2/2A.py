###############################################################################
###############################################################################

# Program: 2A.py
# Author: Sergio Cámara Peña
# Date: 04/08/2025
# Version: FINAL

###############################################################################
###############################################################################

# %% Import libraries
import os
import scanpy as sc
import anndata
import numpy as np
import random
import pandas as pd

# %% Set random seed
random.seed(2504)

# %% Load data
path = "/home/scamara/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Raw_Atlas"
file = f"{path}/Python_scVI_adata_big_V4_state4.h5ad"
adata = anndata.read_h5ad(file)
print(adata.shape[0])

# %% Normalize
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Optionally, set raw to keep a copy of the normalized data
adata.raw = adata

# %% Filter object
# Antigen == "Blood"
adata_filtered = adata[adata.obs['Antigen'] == "Blood"].copy()
print(adata_filtered.shape[0])

# manual_celltype_annotation_high == "CD8 cytotoxic"
adata_filtered = adata_filtered[adata_filtered.obs['manual_celltype_annotation_high'] == "CD8 cytotoxic"].copy()
print(adata_filtered.shape[0])

# STATUS is not "Healthy"
adata_filtered = adata_filtered[adata_filtered.obs['STATUS'] != "HEALTHY"].copy()
print(adata_filtered.shape[0])

# Remove rows where Max_Response is NA
adata_filtered = adata_filtered[~adata_filtered.obs['Max_Response'].isna()].copy()

# Remove cells where Time_Point_Ranges == "Infusion_Product" and Stimulated == "YES"
mask = ~((adata_filtered.obs['Time_Point_Ranges'] == "Infusion_Product") & (adata_filtered.obs['Stimulated'] == "YES"))
adata_filtered = adata_filtered[mask].copy()
print(adata_filtered.shape[0])

adata_filtered_bis = adata_filtered.copy()

# %% Plot figure 2A
color_maps = {
    'Max_Response': {
        'CR': '#006400',
        'PR': '#66c2a5',
        'NR': '#fc8d62'
    }
}

# Convert to an ordered categorical factor
adata_filtered_bis.obs['Max_Response'] = pd.Categorical(
    adata_filtered_bis.obs['Max_Response'],
    categories=['CR', 'PR', 'NR'],
    ordered=True
)

# Visualize using the predefined color mapping
os.chdir("/home/scamara/data/scamara/Atlas_Mieloma_Multiple/Resultados_Figuras/Figura_2")
sc.pl.umap(
    adata_filtered_bis,
    color="Max_Response",
    palette=[color_maps['Max_Response'][cat] for cat in adata_filtered_bis.obs['Max_Response'].cat.categories],
    frameon=False,
    save="Figura_2A.pdf"
)

# %% End of script
