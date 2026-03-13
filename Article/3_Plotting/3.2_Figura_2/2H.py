###############################################################################
###############################################################################

# Program: 2H.py
# Author: Sergio Cámara Peña
# Date: 10/06/2025
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

# %% Set a random seed
random.seed(2504)

# %% Load data
os.chdir("/home/scamara/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Raw_Atlas")
adata_orig = sc.read_h5ad("Python_scVI_adata_big_V4_state4.h5ad")
adata_orig_normalized = sc.read_h5ad("Python_scVI_adata_big_V4_state4_Normalized.h5ad")

# %% Set path to save figs
os.chdir("/home/scamara/data/scamara/Atlas_Mieloma_Multiple/Resultados_Figuras/Figura_2")

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

# Fixed colour dictionary
palette_fija = dict(zip(categorias, palette))

# %% Figure 2H - Celltype composition: Female vs Male
adata4 = adata_orig.copy()
adata4 = adata4[adata4.obs["Antigen"] == "Blood"]

metadata = adata4.obs.copy()
metadata = metadata[metadata["Time_Point_Ranges"] == "Infusion_Product"]
metadata = metadata[metadata["STATUS"] == "DISEASE"]
metadata = metadata[metadata["Stimulated"] == "NO"]
metadata = metadata[metadata["Max_Response"] != "PR"]
metadata["Max_Response"] = pd.Categorical(
    metadata["Max_Response"], categories=["NR", "CR"], ordered=True
)

## General graph
metadata0 = metadata.copy()

mask = (
    (metadata0["Age_Range"] == "20-40")
    | (metadata0["Age_Range"] == "<20")
)

metadata0 = metadata0[~mask]

proportions_General = (
    metadata0.groupby(["manual_celltype_annotation_high", "Max_Response", "Sex"])
    .size()
    .reset_index(name="count")
)

proportions_General["proportion"] = proportions_General.groupby(["Sex", "Max_Response"])[
    "count"
].transform(lambda x: x / x.sum())

p_General = (
    ggplot(
        proportions_General,
        aes(x="Max_Response", y="proportion", fill="manual_celltype_annotation_high"),
    )
    + geom_bar(stat="identity")
    + facet_wrap("Sex")
    + theme_classic()
    + ggtitle('')
    + xlab("Max Response")
    + ylab("Proportion")
    + guides(fill=guide_legend(title="Celltype Annotation"))
    + scale_fill_manual(values=palette_fija)
)

os.chdir(
    "/home/scamara/data/scamara/Atlas_Mieloma_Multiple/Resultados_Figuras/Figura_2"
)

ggsave(p_General,"2H.pdf")

print(p_General)

## Get the values in proportions just to check
p_General # Transcriptomics of general group
NR_females_General = proportions_General[(proportions_General['Sex'] == 'F') & (proportions_General['Max_Response'] == 'NR')]
CR_females_General = proportions_General[(proportions_General['Sex'] == 'F') & (proportions_General['Max_Response'] == 'CR')]
NR_males_General = proportions_General[(proportions_General['Sex'] == 'M') & (proportions_General['Max_Response'] == 'NR')]
CR_males_General = proportions_General[(proportions_General['Sex'] == 'M') & (proportions_General['Max_Response'] == 'CR')]

# %% End of script
