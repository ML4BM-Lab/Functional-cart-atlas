###############################################################################
###############################################################################

# Program: 1F.py
# Author: Sergio Cámara Peña
# Date: 07/01/2026
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

# Fixed color dict
palette_fija = dict(zip(categorias, palette))

# %% Set path to save figs
os.chdir("/home/scamara/data/scamara/Atlas_Mieloma_Multiple/Resultados_Figuras/Figura_1")

########################################################################################################################
########################################################################################################################

# %% 1F - Stacked barplot with cell cycle phases
# Cluster proportion calculation and plotting
adata2 = adata_orig_normalized.copy()
metadata = adata2.obs.copy()

proportions = (
    metadata.groupby(
        ["ScFv", "manual_celltype_annotation_high"]
    )
    .size()
    .reset_index(name="count")
)
proportions["proportion"] = proportions.groupby(["ScFv"])[
    "count"
].transform(lambda x: x / x.sum())

proportions["manual_celltype_annotation_high"] = pd.Categorical(
    proportions["manual_celltype_annotation_high"],
    categories=list(palette_fija.keys()),
    ordered=True
)

scfv_order = ["BCMA", "CD19", "APRIL", "GD2", "HER2"]
proportions["ScFv"] = pd.Categorical(
    proportions["ScFv"],
    categories=scfv_order,
    ordered=True
)

plot = (
    ggplot(
        proportions,
        aes(x="ScFv", y="proportion", fill="manual_celltype_annotation_high"),
    )
    + geom_bar(stat="identity")
    + scale_fill_manual(values=palette_fija)
    + theme_classic()
    + ggtitle("")
    + xlab("")
    + ylab("Proportion")
    + theme(axis_text_x=element_text(rotation=90, ha='right'))
)

# Save it
ggsave(plot, "1F_HEALTHY_Stacked_barplot_ScFv.pdf")

########################################################################################################################
########################################################################################################################

# %% End of script
