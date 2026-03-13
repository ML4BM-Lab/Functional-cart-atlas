###############################################################################
###############################################################################
###############################################################################
###############################################################################

# Program: 4B.py
# Author: Sergio Cámara Peña
# Date: 29/05/2025
# Version: FINAL

###############################################################################
###############################################################################
###############################################################################
###############################################################################

# %% Load all the needed libraries
import scanpy as sc
import os
import random
import seaborn as sns
import pandas as pd

# %% Set a random seed
random.seed(2504)

# %% Read data
os.chdir("/home/scamara/data_a/scamara/Atlas/Input")
adata = sc.read("Python_scVI_adata_big_V4_state4.h5ad")

# %% PATH to save figs
os.chdir("/home/scamara/data_a/scamara/Atlas/Figuras/Figura_4")

# %% Normalize and log data
# 1. Normalize total counts per cell
sc.pp.normalize_total(adata, target_sum=1e4)

# 2. Log-transform the data
sc.pp.log1p(adata)

# 3. Optionally, set raw to keep a copy of the normalized data
adata.raw = adata
adata_orig = adata.copy()

# %% Set genes to plot
Genes_IACs_vs_Rest = ["CD3D", "CD8A", "CD4", "CD14", "LYZ", "CXCL16", "CFD", "GRN", "TYROBP"]

# %% Figure 4Ba - Create the dummy column to separate IACs from the rest
adata.obs["IACs_dummy"] = adata.obs["manual_celltype_annotation_high"].apply(
    lambda x: "IACs" if x == "Monocyte-like T cells" else "Rest_of_Cells"
)

print(adata.obs["IACs_dummy"].value_counts())
print(adata.obs["manual_celltype_annotation_high"].value_counts())

# %% Figure 4Ba - Create and save plots
sc.pl.stacked_violin(adata, var_names=Genes_IACs_vs_Rest, groupby='IACs_dummy', title="", use_raw=True, save="4Ba.pdf", vmax=3)

# %% Figure 4Bb - Filter and only keep IACs
adata2 = adata_orig.copy()
adata2 = adata2[adata2.obs["manual_celltype_annotation_high"] == "Monocyte-like T cells"]
adata2

# %% Figure 4Bb - Create the dummy column to separate pre- and post-Infusion
adata2.obs["IACs_dummy"] = adata2.obs["Time_Point_Ranges"].apply(
    lambda x: "Infusion Product" if x == "Infusion_Product" else "Post Infusion"
)

print(adata2.obs["IACs_dummy"].value_counts())
print(adata2.obs["Time_Point_Ranges"].value_counts())

# %% Figure 4Bb - Create and save plots
sc.pl.stacked_violin(adata2, var_names=Genes_IACs_vs_Rest, groupby='IACs_dummy', title="", use_raw=True, save="4Bb.pdf", vmax=3)

# %% End of script
