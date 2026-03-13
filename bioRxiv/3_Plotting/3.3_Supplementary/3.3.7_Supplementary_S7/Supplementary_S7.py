###############################################################################
###############################################################################

# Program: Supplementary_S7.py
# Author: Sergio Cámara Peña
# Date: 28/01/2025
# Version: FINAL

###############################################################################
###############################################################################

# %% Load all the needed libraries
import scanpy as sc
import os
import random
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from plotnine import *
import milopy.core as milo
import milopy.plot as milopl
import milopy.utils
import scipy.sparse as sp
import anndata as ad

# %% Set a random seed
random.seed(2504)

# %% Set palette
orig_ident_palette = {
    'Bai_et_al': "#F8766D",
    'Boroughs_et_al': "#00BA38",
    'Lynn_et_al': "#FF9900",
    'Rodriguez_Marquez_et_al': "#619CFF",
    'Wang_et_al': "#F564E3",
    'Xhangolli_et_al': "#B79F00",
    'Deng_et_al': "#A81818",
    'Good_et_al': "#9D00FF",
    'Haradvala_et_al': "#006400",
    'Li_X_Cancer_Cell_letter_et_al': "#00BFC4",
    'Li_X_et_al': "#00008A",
    'Melenhorst_et_al': "#753900",
    'Sheih_et_al': "#EBFF0F"
}

# %% Read ATLAS data
os.chdir(
    "/home/scamara/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Raw_Atlas"
)
adata_orig = sc.read_h5ad("Python_scVI_adata_big_V4_state4.h5ad")

# %% Set PATH to save figs
os.chdir("/home/scamara/data/scamara/Atlas_Mieloma_Multiple/Resultados_Figuras/Suplementarias")

# %% Normalize and log data
# 1. Normalize total counts per cell
sc.pp.normalize_total(adata_orig, target_sum=1e4)

# 2. Log-transform the data
sc.pp.log1p(adata_orig)

# 3. Optionally, set raw to keep a copy of the normalized data
adata_orig.raw = adata_orig

########################################################################################################################
########################################################################################################################

# %% S7A - Complete short time point comparison heatmap

##### In "Supplementary_S7_AE.R" script #####

########################################################################################################################
########################################################################################################################

# %% S7B - Artificial doublets vs IACs
adataB = adata_orig.copy()

# Select raw count matrix
if "counts" in adataB.layers:
    X = adataB.layers["counts"]
else:
    X = adataB.raw.X

print("Matrix used for doublets: (cells × genes):", X.shape)

# Subset
adata_IACs   = adataB[adataB.obs["manual_celltype_annotation_high"] == "Monocyte-like T cells"].copy()
adata_noIACs = adataB[adataB.obs["manual_celltype_annotation_high"] != "Monocyte-like T cells"].copy()

# Extract matrix
X_IACs   = adata_IACs.layers["counts"] if "counts" in adata_IACs.layers else adata_IACs.raw.X
X_noIACs = adata_noIACs.layers["counts"] if "counts" in adata_noIACs.layers else adata_noIACs.raw.X

print("IACs matrix shape:", X_IACs.shape)      # (n_IAC_cells, n_genes)
print("No-IACs matrix shape:", X_noIACs.shape) # (n_noIAC_cells, n_genes)

# Number of doublets to generate
n_doublets = 6000
rng = np.random.default_rng(2504)

doublets = []
for i in range(n_doublets):
    # Choose 2 random cells (rows)
    r1, r2 = rng.integers(0, X_noIACs.shape[0], size=2)
    v1 = X_noIACs[r1, :]  # full row = cell
    v2 = X_noIACs[r2, :]
    dbl = (v1 + v2) * 0.5
    dbl = sp.csr_matrix(dbl)
    dbl.data = np.rint(dbl.data)
    dbl.eliminate_zeros()
    doublets.append(dbl)

# Combine rows in only one matrix (artifitial doublets as new cells)
doublets_matrix = sp.vstack(doublets, format="csr")
print("Doublets matrix shape:", doublets_matrix.shape)

# Create AnnData for doublets
var = adataB.var.copy()
obs_doublets = pd.DataFrame(index=[f"Doublet_{i+1}" for i in range(n_doublets)])
obs_doublets["Group"] = "Artificial doublet"

adata_doublets = ad.AnnData(
    X=doublets_matrix,
    var=var,
    obs=obs_doublets
)

# Add tag group to true IACs
adata_IACs.obs["Group"] = "IACs"

# Combine IACs + Doublets
adata_comb = ad.concat([adata_IACs, adata_doublets], axis=0, join="outer")
print("Combined shape:", adata_comb.shape)

# Proprocess + UMAP
sc.pp.normalize_total(adata_comb, target_sum=1e4)
sc.pp.log1p(adata_comb)
sc.pp.highly_variable_genes(adata_comb, n_top_genes=3000, flavor="seurat_v3")
adata_comb = adata_comb[:, adata_comb.var["highly_variable"]].copy()

sc.pp.scale(adata_comb, max_value=10)
sc.tl.pca(adata_comb, n_comps=30)
sc.pp.neighbors(adata_comb, n_neighbors=15, n_pcs=15)
sc.tl.umap(adata_comb, min_dist=0.3)

adata_comb.obs["Group"] = adata_comb.obs["Group"].astype("category")

# Plot UMAP
sc.pl.umap(
    adata_comb,
    color="Group",
    size=5,
    frameon=False,
    legend_loc="right margin",
    title=None,
    save="_S7B_IACs_vs_Doublets.pdf"
)

########################################################################################################################
########################################################################################################################

# %% S7C - IACs Cluster Proportion
## Filter object
adata = adata_orig.copy()

# Filter by 'Antigen' == 'Blood'
adata = adata[adata.obs['Antigen'] == 'Blood']
print(adata.shape[0])

# Further filter by 'manual_celltype_annotation_high' == 'Monocyte-like T cells'
adata = adata[adata.obs['manual_celltype_annotation_high'] == 'Monocyte-like T cells']
print(adata.shape[0])

## Add and modify 'IACs_differenciation' column
# Initialize the column with 'Post_Infusion'
adata.obs['IACs_differenciation'] = 'Post_Infusion'

# Update values to 'Infusion_Product' where 'Time_Point_Ranges' == 'Infusion_Product'
adata.obs.loc[
    adata.obs['Time_Point_Ranges'] == 'Infusion_Product', 'IACs_differenciation'
] = 'Infusion_Product'

## Cluster proportion calculation and plotting
metadata = adata.obs.copy()

proportions = (
    metadata.groupby(
        ["IACs_differenciation", "orig_ident"]
    )
    .size()
    .reset_index(name="count")
)
proportions["proportion"] = proportions.groupby(["IACs_differenciation"])[
    "count"
].transform(lambda x: x / x.sum())

plot = (
    ggplot(
        proportions,
        aes(x="IACs_differenciation", y="proportion", fill="orig_ident"),
    )
    + geom_bar(stat="identity")
    + theme_classic()
    + ggtitle("Cluster Proportion IACs")
    + xlab("")
    + ylab("Proportion")
    + scale_fill_manual(values=orig_ident_palette)
)

## Save it
ggsave(plot, "S7C_IACs_Cluster_Proportion.pdf")

########################################################################################################################
########################################################################################################################

# %% S7D - Dotplot of CXCL10 AND CCL19
adata2 = adata_orig.copy()

## Create the dummy column to separate IACs from the rest
adata2.obs["IACs_dummy"] = adata2.obs["manual_celltype_annotation_high"].apply(
    lambda x: "IACs" if x == "Monocyte-like T cells" else "Rest_of_Cells"
)

print(adata2.obs["IACs_dummy"].value_counts())
print(adata2.obs["manual_celltype_annotation_high"].value_counts())

markers_V2= ["CXCL10", "CCL19"]
sc.pl.dotplot(adata2, markers_V2, use_raw=False, groupby='IACs_dummy', dendrogram=True, figsize=(17, 4.76), save="S7D_Dotplot.pdf")

########################################################################################################################
########################################################################################################################

# %% S7E - IACs vs Rest summarized comparison heatmap

##### In "Supplementary_S7_AE.R" script #####

########################################################################################################################
########################################################################################################################

# %% S7F - GSEA GO dotplot of IACs IP vs Post

##### In "Supplementary_S7F.R" script #####

# %% End of script