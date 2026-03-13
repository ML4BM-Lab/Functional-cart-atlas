###############################################################################
###############################################################################

# Program: 4A.py
# Author: Sergio Cámara Peña
# Date: 27/01/2025
# Version: FINAL

###############################################################################
###############################################################################
###############################################################################
###############################################################################

# %% Load all the needed libraries
import scanpy as sc
import os
import random
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import matplotlib.patches as mpatches

# %% Set a random seed
random.seed(2504)

# %% Read needed files
os.chdir("/home/scamara/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Raw_Atlas")
adata = sc.read_h5ad("Python_scVI_adata_big_V4_state4.h5ad")
print(adata.shape[0])  # Print the number of cells

# %% Filter object
filtered_adata = adata[adata.obs["Antigen"] == "Blood"].copy()
print(filtered_adata.shape[0])  # Print the number of cells

# Display the distribution of 'manual_celltype_annotation_high'
print(filtered_adata.obs["manual_celltype_annotation_high"].value_counts())

# %% Add a new column 'IACs_dummy' based on 'manual_celltype_annotation_high'
filtered_adata.obs["IACs_dummy"] = np.where(
    filtered_adata.obs["manual_celltype_annotation_high"] == "Monocyte-like T cells",
    "IACs",
    "Rest_of_Cells"
)

# Display the distribution of 'IACs_dummy'
print(filtered_adata.obs["IACs_dummy"].value_counts())

# %% Plot UMAP - Figure 4Aa
# Define colors for the categories
os.chdir("/home/scamara/data/scamara/Atlas_Mieloma_Multiple/Resultados_Figuras/Figura_4")
custom_colors = {
    "Rest_of_Cells": "lightgrey",
    "IACs": "#9D00FF"
}

# Map colors to the `IACs_dummy` column
filtered_adata.obs["IACs_color"] = filtered_adata.obs["IACs_dummy"].map(custom_colors)

# Figure
fig, ax = plt.subplots(figsize=(8, 8))
ax.scatter(
    filtered_adata.obsm["X_umap"][:, 0],  # UMAP 1
    filtered_adata.obsm["X_umap"][:, 1],  # UMAP 2
    c=filtered_adata.obs["IACs_color"],  # Color by IACs_dummy mapping
    s=10,  # Size of the points
    alpha=0.8  # Transparency
)
ax.axis("off")  # Remove axes
plt.tight_layout()
plt.savefig("Figura_4Aa.pdf")
plt.show()

# %% Plot UMAP - Figure 4Ab
# Filter and only keep IACs
adata2 = adata.copy()
adata2 = adata2[adata2.obs["manual_celltype_annotation_high"] == "Monocyte-like T cells"]
adata2 = adata2[adata2.obs["STATUS"] == "DISEASE"]
filtered_adata_2 = adata2.copy()

# %% Add a new column 'IACs_dummy' based on 'Time_Point_Ranges'
filtered_adata_2.obs["IACs_dummy"] = np.where(
    filtered_adata_2.obs["Time_Point_Ranges"] == "Infusion_Product",
    "IP",
    "Post"
)

# Set working directory
os.chdir("/home/scamara/data/scamara/Atlas_Mieloma_Multiple/Resultados_Figuras/Figura_4")

# Define colors for the categories
custom_colors = {
    "IP": "#fdae61",
    "Post": "#2c7bb6"
}

# Map colors to the `IACs_dummy` column
filtered_adata_2.obs["IACs_color"] = filtered_adata_2.obs["IACs_dummy"].map(custom_colors)

# Figure
fig, ax = plt.subplots(figsize=(8, 8))
ax.scatter(
    filtered_adata_2.obsm["X_umap"][:, 0],  # UMAP 1
    filtered_adata_2.obsm["X_umap"][:, 1],  # UMAP 2
    c=filtered_adata_2.obs["IACs_color"],   # Color by IACs_dummy mapping
    s=10,
    alpha=0.8
)
ax.axis("off")  # Remove axes
plt.tight_layout()

# Add legend
legend_handles = [
    mpatches.Patch(color=color, label=label)
    for label, color in custom_colors.items()
]
ax.legend(handles=legend_handles, title="IACs", loc="lower left", frameon=False)

# Save and show
plt.savefig("Figura_4Ab.pdf")
plt.show()

# %% End of script