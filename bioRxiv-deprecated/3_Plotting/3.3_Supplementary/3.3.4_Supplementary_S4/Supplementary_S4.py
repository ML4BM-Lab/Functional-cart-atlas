###############################################################################
###############################################################################

# Program: Supplementary_S4.py
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
os.chdir("/home/scamara/data/scamara/Atlas_Mieloma_Multiple/Resultados_Figuras/Suplementarias")

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

########################################################################################################################
########################################################################################################################

# %% S4A - UMAP Red/gray IP analysis
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

adata1 = adata_orig_normalized.copy()
adata1.obs["dummy_variable"] = "No"
adata1.obs.loc[adata1.obs.index.isin(rownames), "dummy_variable"] = "Yes"

# Define colors: strong red for "Yes", grey for "No"
colors = ["grey", "red"]

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
plt.scatter(df_yes["UMAP1"], df_yes["UMAP2"], c="darkred", s=0.3, label="Yes", alpha=0.9)

plt.tight_layout()
plt.savefig("Suplementaria_S4A.pdf", dpi=300)
plt.show()

########################################################################################################################
########################################################################################################################

# %% S4B - UMAP of cell type composition of the red zone
sc.pl.umap(
    adata2, 
    color=["manual_celltype_annotation_high"],
    frameon=False,
    palette=palette_fija,
    save="Suplementaria_S4B.pdf"
)

########################################################################################################################
########################################################################################################################

# %% S4C - Splitted violins of each cell type impairment produced by IL10

##### IN Supplementary_S4C.R script #####

########################################################################################################################
########################################################################################################################

# %% S4D - Celltype composition: Female vs Male
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
    "/home/scamara/data/scamara/Atlas_Mieloma_Multiple/Resultados_Figuras/Suplementarias"
)

ggsave(p_General,"Suplementaria_S4D.pdf")

print(p_General)

## Get the values in proportions just to check
p_General # Transcriptomics of general group
NR_females_General = proportions_General[(proportions_General['Sex'] == 'F') & (proportions_General['Max_Response'] == 'NR')]
CR_females_General = proportions_General[(proportions_General['Sex'] == 'F') & (proportions_General['Max_Response'] == 'CR')]
NR_males_General = proportions_General[(proportions_General['Sex'] == 'M') & (proportions_General['Max_Response'] == 'NR')]
CR_males_General = proportions_General[(proportions_General['Sex'] == 'M') & (proportions_General['Max_Response'] == 'CR')]

########################################################################################################################
########################################################################################################################

# %% S4E - scProportion test summary
##### IN Supplementary_S4E.R script #####

# %% End of script