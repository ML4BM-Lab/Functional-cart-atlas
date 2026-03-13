###############################################################################
###############################################################################

# Program: 1C.py
# Author: Sergio Cámara Peña
# Date: 06/03/2025
# Version: V FINAL

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
import scProportionTest as pt

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
os.chdir("/home/scamara/data/scamara/Atlas_Mieloma_Multiple/Resultados_Figuras/Figura_1")

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

########################################################################################################################
########################################################################################################################

# %% Figure 1C - UMAP cell types
adata1 = adata_orig.copy()

sc.pl.umap(
    adata1, 
    color=["manual_celltype_annotation_high"],
    frameon=False,
    palette=palette_fija,
    save="Figura1C.pdf"
)

# %% End of script
