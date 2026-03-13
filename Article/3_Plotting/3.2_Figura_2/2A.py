###############################################################################
###############################################################################

# Program: 2A.py
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

# Anndata filtering
adata2 = adata_orig_normalized.copy()
adata2 = adata2[adata2.obs["Antigen"] == "Blood"]

adata2 = adata2[
    (adata2.obs["Age_Range"] == "40-60") | (adata2.obs["Age_Range"] == ">60")
]
adata2 = adata2[(adata2.obs["STATUS"] == "DISEASE")]
adata2 = adata2[(adata2.obs["Time_Point_Ranges"] == "Infusion_Product")]
adata2 = adata2[
    (adata2.obs["Max_Response"] == "CR") | (adata2.obs["Max_Response"] == "NR")
]
adata2 = adata2[(adata2.obs["Stimulated"] == "NO")]

rownames = adata2.obs.index.tolist()

########################################################################################################################
########################################################################################################################

# %% Figure 2A - UMAP of cell type composition of the red zone
sc.pl.umap(
    adata2, 
    color=["manual_celltype_annotation_high"],
    frameon=False,
    palette=palette_fija,
    save="_2Aa.pdf"
)

sc.pl.umap(
    adata2,
    color="Max_Response",
    frameon=False,
    palette=["#006400", "#fc8d62"],
    save="_2Ab.pdf"
)

# %% End of script
