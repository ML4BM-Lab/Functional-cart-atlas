###############################################################################
###############################################################################

# Program: Supplementary_S3.py
# Author: Sergio Cámara Peña
# Date: 05/06/2025
# Version: FINAL

###############################################################################
###############################################################################

# %% Load all the needed libraries
import scanpy as sc
import os
import random
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import milopy.core as milo
import milopy.plot as milopl
import milopy.utils
import numpy as np
from plotnine import *

from matplotlib.colors import to_rgb
import matplotlib.patches as mpatches
from collections import OrderedDict

os.chdir(
    "/home/scamara/data/scamara/Atlas_Mieloma_Multiple/Codigo/Codigo_datasets_atlas/Datasets_Integration/Initial_Version_Atlas"
)
from scanpy_cluster_proportions import get_cluster_proportions, plot_cluster_proportions


# %% Set a random seed
random.seed(2504)

# %% Load data
os.chdir("/home/scamara/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Raw_Atlas")
adata_orig = sc.read_h5ad("Python_scVI_adata_big_V4_state4.h5ad")
adata_orig_normalized = sc.read_h5ad("Python_scVI_adata_big_V4_state4_Normalized.h5ad")

# %% Define a 12 color palette
palette = [
"#F8766D",
"#00BA38",
"#B79F00",
"#FF9900",
"#619CFF",
"#F564E3",
"#A81818",
"#9D00FF",
"#006400",
"#00BFC4",
"#00008A"
]

sns.palplot(palette)

# %% Set path to save figs
os.chdir("/home/scamara/data/scamara/Atlas_Mieloma_Multiple/Resultados_Figuras/Suplementarias")

########################################################################################################################
########################################################################################################################

# %% S3A - UMAPs
# Adjust global font size
plt.rcParams.update({
    'font.size': 14,
    'axes.titlesize': 16,
    'axes.labelsize': 14,
    'legend.fontsize': 12,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12
})

adata2 = adata_orig_normalized.copy()

# CD3D
fig = sc.pl.umap(adata2, use_raw=False, color=["CD3D"], frameon=False, ncols=3, vmin=0, vmax="p99", show=False)
fig.axes.title.set_fontsize(16)
fig.axes.xaxis.label.set_fontsize(14)
fig.axes.yaxis.label.set_fontsize(14)
fig.axes.tick_params(labelsize=12)
plt.savefig("S3A_umap_CD3.pdf")
plt.close()

# CD4
fig = sc.pl.umap(adata2, use_raw=False, color=["CD4"], frameon=False, ncols=3, vmin=0, vmax="p99", show=False)
fig.axes.title.set_fontsize(16)
fig.axes.xaxis.label.set_fontsize(14)
fig.axes.yaxis.label.set_fontsize(14)
fig.axes.tick_params(labelsize=12)
plt.savefig("S3A_umap_CD4.pdf")
plt.close()

# CD8A
fig = sc.pl.umap(adata2, use_raw=False, color=["CD8A"], frameon=False, ncols=3, vmin=0, vmax="p99", show=False)
fig.axes.title.set_fontsize(16)
fig.axes.xaxis.label.set_fontsize(14)
fig.axes.yaxis.label.set_fontsize(14)
fig.axes.tick_params(labelsize=12)
plt.savefig("S3A_umap_CD8.pdf")
plt.close()

########################################################################################################################
########################################################################################################################

# %% S3B - Dotplot of curated markers
markers_V2 = ["CD3D", "CD4", "CD8A", "NKG7", "GNLY", "PRF1", "GZMK", "GZMB", "CD69", "CD27", "CD28", "HLA-DRA", "TNFRSF9", "TCF7", "CCR7", "SELL", "FOXP3", "IL2RA", "CTLA4", "IKZF2", "ZWINT", "MKI67", "TUBA1B", "TUBB", "CD68", "CD14", "LYZ"]
sc.pl.dotplot(adata2, markers_V2, use_raw=False, groupby='manual_celltype_annotation_high', dendrogram=True, figsize=(16, 4), save="S3B_Annotation.pdf")

########################################################################################################################
########################################################################################################################

# %% S3C - Stacked barplot with cell cycle phases
# Cluster proportion calculation and plotting
metadata = adata2.obs.copy()

proportions = (
    metadata.groupby(
        ["manual_celltype_annotation_high", "Phase"]
    )
    .size()
    .reset_index(name="count")
)
proportions["proportion"] = proportions.groupby(["manual_celltype_annotation_high"])[
    "count"
].transform(lambda x: x / x.sum())

plot = (
    ggplot(
        proportions,
        aes(x="manual_celltype_annotation_high", y="proportion", fill="Phase"),
    )
    + geom_bar(stat="identity")
    + theme_classic()
    + ggtitle("")
    + xlab("")
    + ylab("Proportion")
    + theme(axis_text_x=element_text(rotation=90, ha='right'))
)

# Save it
ggsave(plot, "S3C_Stacked_barplot_cellcycle_phase_Proportion.pdf")

########################################################################################################################
########################################################################################################################

# %% S3D - Stacked barplot with contribuition of each cluster
## Cluster proportion calculation and plotting
sns.palplot(palette)

categorias = [
    "Apoptotic T cells",
    "CD4 central memory",
    "CD4 cytotoxic",
    "CD4 effector memory",
    "CD8 cytotoxic",
    "CD8 effector memory",
    "CD8 memory",
    "Monocyte-like T cells",
    "Proliferative T cells",
    "Regulatory T cells",
    "Ribosomal enriched"
]

# Fixed colour dictionary
palette_fija = dict(zip(categorias, palette))

color_dict = palette_fija

metadata = adata2.obs.copy()

proportions = (
    metadata.groupby(
        ["manual_celltype_annotation_high", "orig_ident"]
    )
    .size()
    .reset_index(name="count")
)
proportions["proportion"] = proportions.groupby(["orig_ident"])[
    "count"
].transform(lambda x: x / x.sum())

plot2 = (
    ggplot(
        proportions,
        aes(x="orig_ident", y="proportion", fill="manual_celltype_annotation_high"),
    )
    + geom_bar(stat="identity")
    + theme_classic()
    + ggtitle("")
    + xlab("")
    + ylab("Proportion")
    + theme(axis_text_x=element_text(rotation=90, ha='right'))
    + scale_fill_manual(values=color_dict)
)

# Save it
ggsave(plot2, "S3D_Stacked_barplot_cluster_contribution.pdf")

# %% End of script
