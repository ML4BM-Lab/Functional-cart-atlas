###############################################################################
###############################################################################

# Program: 1C.py
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
os.chdir("/home/scamara/data/scamara/Atlas_Mieloma_Multiple/Resultados_Figuras/Figura_1")

########################################################################################################################
########################################################################################################################

# %% Figure 1C - Dotplot of curated markers
adata2 = adata_orig_normalized.copy()
markers_V2 = ["CD3D", "CD4", "CD8A", "NKG7", "GNLY", "PRF1", "GZMK", "GZMB", "CD69", "CD27", "CD28", "HLA-DRA", "TNFRSF9", "TCF7", "CCR7", "SELL", "FOXP3", "IL2RA", "CTLA4", "IKZF2", "ZWINT", "MKI67", "TUBA1B", "TUBB", "CD68", "CD14", "LYZ"]
sc.pl.dotplot(adata2, markers_V2, use_raw=False, groupby='manual_celltype_annotation_high', dendrogram=True, figsize=(16, 4), save="1C_Annotation.pdf")

# %% End of script
