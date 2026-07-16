###############################################################################
###############################################################################

# Program: 2B.py
# Author: Sergio Cámara Peña
# Date: 02/07/2025
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

from matplotlib.colors import to_rgb
import matplotlib.patches as mpatches
from collections import OrderedDict

# %% Set a random seed
random.seed(2504)

# %% Load data
os.chdir("/home/scamara/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Raw_Atlas")
adata_orig_normalized = sc.read_h5ad("Python_scVI_adata_big_V4_state4_Normalized.h5ad")

# %% Set path to save figs
os.chdir("/home/scamara/data/scamara/Atlas_Mieloma_Multiple/Resultados_Figuras/Figura_2")

# %% Figure 2B
# Make a copy of the original object
adata2 = adata_orig_normalized.copy()

# Filter
adata2 = adata2[adata2.obs["manual_celltype_annotation_high"] == "CD8 cytotoxic"]
adata2 = adata2[adata2.obs["STATUS"] != "HEALTHY"]

# Subset the color dictionary to the groups of interest
color_maps = {
    'Time_Point_Ranges': {
        'Infusion_Product': '#fdae61',
        '<2_weeks': '#fee08b',
        '2_weeks-3_months': 'lightgrey',
        '>3_months': 'lightgrey'
    }
}

# Desired order
desired_order = list(color_maps['Time_Point_Ranges'].keys())

# Reassign as an ordered categorical variable
adata2.obs['Time_Point_Ranges'] = pd.Categorical(
    adata2.obs['Time_Point_Ranges'],
    categories=desired_order,
    ordered=True
)

# Color palette in order
palette = [color_maps['Time_Point_Ranges'][cat] for cat in desired_order]

# Plot
sc.pl.umap(
    adata2,
    color="Time_Point_Ranges",
    palette=palette,
    title="",
    save="_2B.pdf"
)

# %% End of script
