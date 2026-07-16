###############################################################################
###############################################################################

# Program: 1I.py
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

# %% Figure 1I - Proportion Test CR Female vs CR Male
adata3 = adata_orig.copy()

## Filter the object
adata3 = adata3[adata3.obs["Antigen"] == "Blood"]
adata_Second = adata3.copy()

# Apply filters directly on adata_Second.obs
adata_Second = adata_Second[
    (adata_Second.obs["Time_Point_Ranges"] == "Infusion_Product") &
    (adata_Second.obs["Stimulated"] == "NO") &
    (adata_Second.obs["STATUS"] == "DISEASE")
].copy()

adata_Second = adata_Second[(adata_Second.obs["Max_Response"] == "CR")].copy()

# Create a mask for the age ranges to exclude
mask = (
    (adata_Second.obs["Age_Range"] == "20-40") |
    (adata_Second.obs["Age_Range"] == "<20")
)

# Apply the mask to filter cells
adata_Second = adata_Second[~mask].copy()

## Make and Save/Load the proportions test
os.chdir("/home/scamara/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Post_Integration/Data")

if False:
    results_1 = pt.permutation_test(adata_Second,"F","M",group_col='Sex',cell_type_col='manual_celltype_annotation_high',nperm=10000,alpha=0.05,n_bootstrap=10000,verbose=True)
    results_1.to_pickle("F_CR_vs_M_CR.pkl")
else:
    results_1 = pd.read_pickle("F_CR_vs_M_CR.pkl")

## Create plots
os.chdir("/home/scamara/data/scamara/Atlas_Mieloma_Multiple/Resultados_Figuras/Figura_1")
pt.point_range_plot(results_1, figsize=(10,8))
plot_1 = plt.gcf()
plot_1.savefig(
    "1I_ProportionTest_F_CR_vs_M_CR.pdf",
    bbox_inches="tight",
    dpi=300
)

# %% End of script
