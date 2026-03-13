###############################################################################
###############################################################################
###############################################################################
###############################################################################

# Program: 5E.py
# Author: Sergio Cámara Peña
# Date: 04/08/2025
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
adata = sc.read("Atlas_integ_scArches_FINAL_V5.h5ad")

# %% Normalize and log data
# 1. Normalize total counts per cell
sc.pp.normalize_total(adata, target_sum=1e4)

# 2. Log-transform the data
sc.pp.log1p(adata)

# 3. Optionally, set raw to keep a copy of the normalized data
adata.raw = adata

# %% Filters:
# 1. Filter only cells where 'Antigen' is "Blood"
adata_filtered = adata[adata.obs['Antigen'] == "Blood"].copy()
print(adata_filtered.n_obs)

# 2. Filter to remove cells where 'Time_Point_Ranges' == "Infusion_Product" and 'Stimulated' == "YES"
mask = ~((adata_filtered.obs['Time_Point_Ranges'] == "Infusion_Product") & (adata_filtered.obs['Stimulated'] == "YES"))
adata_filtered = adata_filtered[mask].copy()
print(adata_filtered.n_obs)

# 3. Filter to remove cells with NA value in 'Max_Response'
adata_filtered = adata_filtered[~adata_filtered.obs['Max_Response'].isna()].copy()

# 4. Filter to keep only cells with 'Max_Response' == "CR"
adata_filtered = adata_filtered[adata_filtered.obs['Max_Response'] == "CR"].copy()
print(adata_filtered.n_obs)

# 5. Filter to keep only cells with 'STATUS' == "DISEASE" - NOT NEEDED
adata_filtered = adata_filtered[adata_filtered.obs['STATUS'] == "DISEASE"].copy()
print(adata_filtered.n_obs)

# 6. Filter to remove cells with 'Stimulation_Location' == "In_vitro" - NOT NEEDED
adata_filtered = adata_filtered[adata_filtered.obs['Stimulation_Location'] != "In_vitro"].copy()
print(adata_filtered.n_obs)

# 7. Filter to keep only cells with 'Time_Point_Ranges' == "2_weeks-3_months"
adata_filtered = adata_filtered[adata_filtered.obs['Time_Point_Ranges'] == "2_weeks-3_months"].copy()
print(adata_filtered.n_obs)

# %% Create and save plots
os.chdir("/home/scamara/data_a/scamara/Atlas/Figuras/Figura_5")

Genes_cyto = ["CD7", "IFNG", "GZMK", "GZMH", "CD8A", "GZMA", "CCL4", "NKG7", "CCL5", "GNLY", "PRF1", "GZMB"]

# Filter scFv == "BCMA"
adata_scFv1 = adata_filtered[adata_filtered.obs['ScFv'] == "BCMA"].copy()
sc.pl.stacked_violin(adata_scFv1, var_names=Genes_cyto, groupby='Product_norm', use_raw=True, vmax=5, save="5E_BCMA.pdf")

# Filter scFv == "CD19"
adata_scFv2 = adata_filtered[adata_filtered.obs['ScFv'] == "CD19"].copy()
sc.pl.stacked_violin(adata_scFv2, var_names=Genes_cyto, groupby='Product_norm', use_raw=True, vmax=5, save="5E_CD19.pdf")

# %% End of script
