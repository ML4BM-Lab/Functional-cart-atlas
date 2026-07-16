###############################################################################
###############################################################################

# Program: Supplementary_S1.py
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

from matplotlib.colors import to_rgb
import matplotlib.patches as mpatches
from collections import OrderedDict

# %% Set a random seed
random.seed(2504)

# %% Load data
os.chdir("/home/scamara/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Raw_Atlas")
adata_orig = sc.read_h5ad("Python_scVI_adata_big_V4_state4.h5ad")

# %% Set path to save figs
os.chdir("/home/scamara/data/scamara/Atlas_Mieloma_Multiple/Resultados_Figuras/Suplementarias")

# %% Set palette
color_maps = {
    'STATUS': {'HEALTHY': '#a6d854', 'DISEASE': '#fc8d62'},
    'Time_Point_Ranges': {
        'Infusion_Product': '#fdae61',
        '<2_weeks': '#fee08b',
        '2_weeks-3_months': '#f46d43',
        '>3_months': '#d73027'
    },
    'Age_Range': {
    '<20':   '#a6d854',
    '20-40': '#66c2a5',
    '40-60': '#fc8d62',
    '>60':   '#e31a1c'
    },
    'Sex': {'F': '#e5c1ff', 'M': '#9ecae1'},
    'ScFv': {
        'CD19': '#66c2a5',
        'BCMA': '#fc8d62',
        'APRIL': '#8da0cb',
        'GD2': '#e78ac3',
        'HER2': '#a6d854'
    },
    'Max_Response': {
        'CR': '#006400',
        'PR': '#66c2a5',
        'NR': '#fc8d62'
    },
    'ICANS_Grade_Range': {
        '0': '#a6d854',
        '1-2': '#fdae61',
        '3-4': '#d73027'
    },
    'orig_ident': {
        'Bai_et_al': "#F8766D",
        'Boroughs_et_al': "#00BA38",
        'Lynn_et_al': "#FF9900",
        'Rodriguez_Marquez_et_al': "#619CFF",
        'Wang_et_al': "#F564E3",
        'Xhangolli_et_al': "#B79F00",
        'Deng_et_al': "#A81818",
        'Good_et_al': "#9D00FF",
        'Haradvala_et_al': "#006400",
        'Li_X_Cancer_Cell_letter_et_al': "#00BFC4",
        'Li_X_et_al': "#00008A",
        'Melenhorst_et_al': "#753900",
        'Sheih_et_al': "#EBFF0F"

}
}

# Desired order
status_order = ["DISEASE", "HEALTHY"]
time_order = [
    "Infusion_Product",
    "<2_weeks",
    "2_weeks-3_months",
    ">3_months"
]
age_order = [
    "<20", "20-40", "40-60", ">60"
]

response_order = [
    "NR", "PR", "CR"
]

# Transform into ordered categories
adata_orig.obs["STATUS"] = pd.Categorical(adata_orig.obs["STATUS"], categories=status_order, ordered=True)
adata_orig.obs["Time_Point_Ranges"] = pd.Categorical(adata_orig.obs["Time_Point_Ranges"], categories=time_order, ordered=True)
adata_orig.obs["Age_Range"] = pd.Categorical(adata_orig.obs["Age_Range"], categories=age_order, ordered=True)
adata_orig.obs["Max_Response"] = pd.Categorical(adata_orig.obs["Max_Response"], categories=response_order, ordered=True)

########################################################################################################################
########################################################################################################################

# %% Plot UMAPs
sc.pl.umap(
    adata_orig, 
    color=["orig_ident"],
    frameon=False,
    palette=color_maps["orig_ident"],
    save="_Suplementaria_S1_Orig_Ident.pdf")

sc.pl.umap(
    adata_orig, 
    color=["STATUS"],
    frameon=False,
    palette=color_maps["STATUS"],
    save="_Suplementaria_S1_STATUS.pdf")

sc.pl.umap(
    adata_orig, 
    color=["Time_Point_Ranges"],
    frameon=False,
    palette=color_maps["Time_Point_Ranges"],
    save="_Suplementaria_S1_Time_Point_Ranges.pdf")

sc.pl.umap(
    adata_orig, 
    color=["ScFv"],
    frameon=False,
    palette=color_maps["ScFv"],
    title="CAR Type",
    save="_Suplementaria_S1_CAR_Type.pdf")

sc.pl.umap(
    adata_orig, 
    color=["Max_Response"],
    frameon=False,
    palette=color_maps["Max_Response"],
    save="_Suplementaria_S1_Max_Response.pdf")

sc.pl.umap(
    adata_orig, 
    color=["ICANS_Grade_Range"],
    frameon=False,
    palette=color_maps["ICANS_Grade_Range"],
    save="_Suplementaria_S1_ICANS_Grade_Range.pdf")

sc.pl.umap(
    adata_orig, 
    color=["Age_Range"],
    frameon=False,
    palette=color_maps["Age_Range"],
    save="_Suplementaria_S1_Age_Range.pdf")

sc.pl.umap(
    adata_orig, 
    color=["Sex"],
    frameon=False,
    palette=color_maps["Sex"],
    save="_Suplementaria_S1_Sex.pdf")

# %% End of script