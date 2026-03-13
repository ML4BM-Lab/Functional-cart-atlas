###############################################################################
###############################################################################

# Program: 3A.py
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
os.chdir("/home/scamara/data/scamara/Atlas_Mieloma_Multiple/Resultados_Figuras/Figura_3")

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

# %% Figure 3A:
##### UMAP general grey - CD8 cytotoxic #####
adata1 = adata_orig_normalized.copy()

adata1.obs["dummy_variable"] = "No"
adata1.obs.loc[
    adata1.obs["manual_celltype_annotation_high"] == "CD8 cytotoxic",
    "dummy_variable"
] = "Yes"

# Make sure the categories are in the correct order
adata1.obs["dummy_variable"] = adata1.obs["dummy_variable"].astype("category")
adata1.obs["dummy_variable"].cat.set_categories(["No", "Yes"], inplace=True)

## Plot UMAP
# Get UMAP coordinates
umap = adata1.obsm["X_umap"]
df = adata1.obs.copy()
df["UMAP1"] = umap[:, 0]
df["UMAP2"] = umap[:, 1]

# Plot grey dots ("No") first
plt.figure(figsize=(6, 6))
df_no = df[df["dummy_variable"] == "No"]
plt.scatter(df_no["UMAP1"], df_no["UMAP2"], c="lightgrey", s=0.3, label="No", alpha=0.6)

# Plot red dots ("Yes") on top
df_yes = df[df["dummy_variable"] == "Yes"]
plt.scatter(df_yes["UMAP1"], df_yes["UMAP2"], c=palette[4], s=0.3, label="Yes", alpha=0.9)

plt.tight_layout()
plt.savefig("Figura_3Aa.pdf", dpi=300)
plt.show()

##### UMAP CD8 cytotoxic - Time Point #####
adata2 = adata_orig_normalized.copy()
adata2 = adata2[adata2.obs["manual_celltype_annotation_high"] == "CD8 cytotoxic"]
adata2 = adata2[adata2.obs["STATUS"] == "DISEASE"]

color_maps = {
    'Time_Point_Ranges': {
        'Infusion_Product': '#fdae61',
        '<2_weeks': '#fee08b',
        '2_weeks-3_months': '#f46d43',
        '>3_months': "#c80d038a"
    }
}

# Define the correct order based on the dictionary
desired_order = list(color_maps['Time_Point_Ranges'].keys())

# Make sure the variable is categorical and with the desired order
adata2.obs['Time_Point_Ranges'] = pd.Categorical(
    adata2.obs['Time_Point_Ranges'],
    categories=desired_order,
    ordered=True
)

# Extract the color palette in the same order
palette = [color_maps['Time_Point_Ranges'][cat] for cat in desired_order]

# UMAP with the palette in the correct order
sc.pl.umap(
    adata2,
    color="Time_Point_Ranges",
    palette=palette,
    title="",
    save="_Figura_3Ab.pdf"
)

# %% End of script
