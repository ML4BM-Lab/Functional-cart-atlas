###############################################################################
###############################################################################

# Program: 2B&C.py
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

# Fixed color dict
palette_fija = dict(zip(categorias, palette))

########################################################################################################################
########################################################################################################################

# %% Figure 2B and 2C - Milopy Complete Response and Memory Phenotype
# https://nbviewer.org/github/emdann/milopy/blob/master/notebooks/milopy_example.ipynb

adata2 = adata_orig_normalized.copy()
adata2 = adata2[adata2.obs["Antigen"] == "Blood"]

adata2 = adata2[
    (adata2.obs["Age_Range"] == "40-60") | (adata2.obs["Age_Range"] == ">60")
]

adata2 = adata2[(adata2.obs["Time_Point_Ranges"] == "Infusion_Product")]

adata2 = adata2[
    (adata2.obs["Max_Response"] == "CR") | (adata2.obs["Max_Response"] == "NR")
]

adata2 = adata2[(adata2.obs["Stimulated"] == "NO")]
print(adata2)

d = 30
k = 60

print(adata2.obs["Age_Range"].value_counts())
print(adata2.obs["Sex"].value_counts())
print(adata2.obs["Max_Response"].value_counts())

sc.pp.neighbors(adata2, use_rep="X_scVI", n_neighbors=k, n_pcs=d)

milo.make_nhoods(adata2)

nhood_size = np.array(adata2.obsm["nhoods"].sum(0)).ravel()
plt.hist(
    nhood_size, bins=100
)  # The peak should be around X samples we are using × 3, to ensure at least 3 cells per sample per neighbourhood
plt.axvline(
    len(adata2.obs["Product_norm"].unique()) * 3,
    color="black",
    linestyle="dashed",
    linewidth=1,
)

milo.count_nhoods(adata2, sample_col="Product_norm")

##### Run milopy ~Age_Range + Sex + Max_Response #####
milo.DA_nhoods(
    adata2, design="~Age_Range + Sex + Max_Response"
)  # Design: A formula or model.matrix object describing the experimental design for differential abundance testing. The last component of the formula or last column of the model matrix are by default the test variable. This behaviour can be overridden by setting the model.contrasts argument
adata2.uns["nhood_adata"].obs

# Diagnostic plots
old_figsize = plt.rcParams["figure.figsize"]
plt.rcParams["figure.figsize"] = [10, 5]
plt.subplot(1, 2, 1)
plt.hist(adata2.uns["nhood_adata"].obs.PValue, bins=50)
plt.xlabel("P-Vals")
plt.subplot(1, 2, 2)
plt.plot(
    adata2.uns["nhood_adata"].obs.logFC,
    -np.log10(adata2.uns["nhood_adata"].obs.SpatialFDR),
    ".",
)
plt.xlabel("log-Fold Change")
plt.ylabel("- log10(Spatial FDR)")
plt.tight_layout()
plt.rcParams["figure.figsize"] = old_figsize

milopy.utils.build_nhood_graph(adata2)

# Visualize results on embedding
plt.rcParams["figure.figsize"] = [10, 10]

milopl.plot_nhood_graph(
    adata2, alpha=0.01, min_size=2, save="Figura2B.pdf"  ## SpatialFDR level (1%) --- Size of smallest dot
)

# Visualize result by celltype
milopy.utils.annotate_nhoods(adata2, anno_col="manual_celltype_annotation_high")

# We can see that for the majority of neighbourhoods, almost all cells have the same neighbourhood. We can rename neighbourhoods where less than 60% of the cells have the top label as "Mixed"
plt.hist(adata2.uns["nhood_adata"].obs["nhood_annotation_frac"])
plt.xlabel("celltype fraction")
adata2.uns["nhood_adata"].obs.loc[
    adata2.uns["nhood_adata"].obs["nhood_annotation_frac"] < 0.6, "nhood_annotation"
] = "Mixed"

# Copy the working DataFrame and filter
df = adata2.uns["nhood_adata"].obs.copy()
excluir = ["Mixed", "Monocyte-like T cells", "Ribosomal enriched"]
df = df[~df["nhood_annotation"].isin(excluir)]

# Delete non-used levels
df["nhood_annotation"] = df["nhood_annotation"].astype("category")
df["nhood_annotation"] = df["nhood_annotation"].cat.remove_unused_categories()

# Sort (optional)
order = df.groupby("nhood_annotation")["logFC"].median().sort_values().index

# Plot
plt.figure(figsize=(8, 6))
sns.violinplot(
    data=df,
    y="nhood_annotation",
    x="logFC",
    order=order,
    cut=0,
    scale="width",
    linewidth=1,
    inner=None,
    palette=palette_fija
)

sns.stripplot(x="logFC", y="nhood_annotation", data=df, jitter=True, size = 1.5, color="black", order=order)

plt.axvline(x=0, color="black", linestyle="--")
plt.xlabel("logFC")
plt.ylabel("")
plt.tight_layout()
plt.savefig("Figura2C.pdf", dpi=300)

# %% End of script
