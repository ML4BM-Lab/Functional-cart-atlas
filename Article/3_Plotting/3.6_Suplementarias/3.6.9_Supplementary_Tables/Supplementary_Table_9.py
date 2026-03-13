###############################################################################
###############################################################################
###############################################################################
###############################################################################

# Program: Supplementary_Table_9.py
# Author: Sergio Cámara Peña
# Date: 15/09/2025
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

Genes_IACs_clusters = ["AIF1", "LTA4H", "LST1"]

# %% Apply statistics
from scipy.stats import mannwhitneyu

def get_significance(p):
    if p < 0.001:
        return '***'
    elif p < 0.01:
        return '**'
    elif p < 0.05:
        return '*'
    else:
        return 'ns'

results = []

for gene in Genes_IACs_clusters:
    expr = adata2.raw[:, gene].X.toarray().flatten()
    groups = adata2.obs["IACs_dummy"]
    
    # Separate the groups
    group1 = expr[groups == "Infusion Product"]
    group2 = expr[groups == "Post Infusion"]
    
    # Mann-Whitney test
    stat, p = mannwhitneyu(group1, group2, alternative='two-sided')
    
    results.append({
        "Gene": gene,
        "U statistic": stat,
        "p-value": p
    })

# Convert into DataFrame for better visualization.
results_df = pd.DataFrame(results)
results_df["Significance"] = results_df["p-value"].apply(get_significance)
print(results_df)

# %% End of script