###############################################################################
###############################################################################
###############################################################################
###############################################################################

# Program: 2H.py
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

# %% Normalize and log data
# 1. Normalize total counts per cell
sc.pp.normalize_total(adata, target_sum=1e4)

# 2. Log-transform the data
sc.pp.log1p(adata)

# 3. Optionally, set raw to keep a copy of the normalized data
adata.raw = adata

# %% Filter and only keep IACs
adata2 = adata.copy()
adata2 = adata2[adata2.obs["manual_celltype_annotation_high"] == "Monocyte-like T cells"]
adata2 = adata2[adata2.obs["STATUS"] == "DISEASE"]
adata2
del adata

# %% Create the dummy column to separate IACs from the rest
adata2.obs["IACs_dummy"] = adata2.obs["Time_Point_Ranges"].apply(
    lambda x: "Infusion Product" if x == "Infusion_Product" else "Post Infusion"
)

print(adata2.obs["IACs_dummy"].value_counts())
print(adata2.obs["Time_Point_Ranges"].value_counts())

# %% Create and save plots
os.chdir("/home/scamara/data_a/scamara/Atlas/Figuras/Figura_2")

Genes_IACs_clusters = ["AIF1", "LTA4H", "LST1"]

sc.pl.stacked_violin(adata2, var_names=Genes_IACs_clusters, groupby='IACs_dummy', title="", use_raw=True, save="2H.pdf")

# %% End of script