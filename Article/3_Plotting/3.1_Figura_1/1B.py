###############################################################################
###############################################################################

# Program: 1B.py
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

# %% Figure 1Ba - Categorical heatmap
adata1_bis = adata_orig.copy()
adata1_bis.obs.rename(columns={"ScFv": "CAR_Type"}, inplace=True)

adata1_bis.obs = adata1_bis.obs[[
    "Product_norm", "STATUS", "Time_Point_Ranges", "Age_Range",
    "Sex", "CAR_Type", "Max_Response", "ICANS_Grade_Range"
]]

Df_heatmap = adata1_bis.obs.copy()

# 1. Reset the index so cell codes are no longer used as the index
Df_heatmap.reset_index(inplace=True)

# 2. Set 'Product_norm' as the new index
Df_heatmap.set_index('Product_norm', inplace=True)

# 3. Remove duplicates in the index (keeping only the first occurrence)
Df_heatmap = Df_heatmap[~Df_heatmap.index.duplicated(keep='first')]

cols = ["STATUS", "Time_Point_Ranges", "CAR_Type", "Max_Response",
        "ICANS_Grade_Range", "Age_Range", "Sex"]

color_maps = {
    'STATUS': {'HEALTHY': '#a6d854', 'DISEASE': '#fc8d62'},
    'Time_Point_Ranges': {
        'Infusion_Product': '#fdae61',
        '<2_weeks': '#fee08b',
        '2_weeks-3_months': '#f46d43',
        '>3_months': '#d73027'
    },
    'Age_Range': {
        '<20': '#c6dbef',
        '20-40': '#9ecae1',
        '40-60': '#4292c6',
        '>60': '#08519c'
    },
    'Sex': {'F': '#e5c1ff', 'M': '#9ecae1'},
    'CAR_Type': {
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
    }
}

# Desired order
status_order = ["DISEASE", "HEALTHY"]
time_order = ["Infusion_Product", "<2_weeks", "2_weeks-3_months", ">3_months"]
age_order = ["<20", "20-40", "40-60", ">60"]
response_order = ["NR", "PR", "CR"]

# Convert to ordered categorical
Df_heatmap["STATUS"] = pd.Categorical(Df_heatmap["STATUS"], categories=status_order, ordered=True)
Df_heatmap["Time_Point_Ranges"] = pd.Categorical(Df_heatmap["Time_Point_Ranges"], categories=time_order, ordered=True)
Df_heatmap["Age_Range"] = pd.Categorical(Df_heatmap["Age_Range"], categories=age_order, ordered=True)
Df_heatmap["Max_Response"] = pd.Categorical(Df_heatmap["Max_Response"], categories=response_order, ordered=True)

Df_heatmap = Df_heatmap.rename(columns={'index': 'original_index'})
Df_heatmap = Df_heatmap.reset_index()
Df_heatmap = Df_heatmap.rename(columns={'index': 'Product_norm'})

# Order the DataFrame
Df_heatmap = Df_heatmap.sort_values(
    ["STATUS", "Time_Point_Ranges", "CAR_Type", "Max_Response",
     "ICANS_Grade_Range", "Age_Range", "Sex"]
).reset_index(drop=True)

# Convert categorical columns to string to avoid errors
for col in cols:
    Df_heatmap[col] = Df_heatmap[col].astype(str)

# Map colors, values not present in the dictionary will return NaN
color_df = pd.DataFrame(index=Df_heatmap.index)
for col in cols:
    color_df[col] = Df_heatmap[col].map(color_maps[col])

# Replace NaN values in color_df with light gray
COLOR_NA = '#d9d9d9'
color_df = color_df.fillna(COLOR_NA)

# Convert hex colors to RGB arrays
rgb_matrix = color_df[cols].T.applymap(to_rgb)
rgb_array = np.array([[list(x) for x in row] for row in rgb_matrix.to_numpy()], dtype=float)

# Plot heatmap
fig, ax = plt.subplots(figsize=(len(Df_heatmap) * 1.2, len(cols) * 0.3 + 2))
ax.imshow(rgb_array, aspect=3, interpolation='nearest', rasterized=True)  # 👈 rasterized evita PDF corrupto

ax.set_yticks(range(len(cols)))
ax.set_yticklabels(cols, fontsize=10)
ax.set_xticks(range(len(Df_heatmap)))

# Clean borders
ax.set_xticks([], minor=True)
ax.set_yticks([], minor=True)
ax.set_xticklabels([])
ax.spines[:].set_visible(False)

from matplotlib.offsetbox import AnchoredOffsetbox, HPacker, VPacker, DrawingArea, TextArea

# Create sub-legends by column
col_blocks = []
for col in cols:
    items = []
    label_title = TextArea(f"{col}", textprops=dict(size=9, weight='bold'))

    unique_vals = pd.Series(Df_heatmap[col].unique())
    valid_vals = [v for v in unique_vals if pd.notna(v) and v in color_maps[col]]
    has_missing = unique_vals.isna().any() or any(v not in color_maps[col] for v in unique_vals if pd.notna(v))

    for val in valid_vals:
        clr = color_maps[col][val]
        da = DrawingArea(15, 10, 0, 0)
        rect = plt.Rectangle((0, 0), 15, 10, facecolor=clr, edgecolor='black', lw=0.3)
        da.add_artist(rect)
        label = TextArea(f" {val}", textprops={'size': 8})
        packed = HPacker(children=[da, label], align="center", pad=0, sep=2)
        items.append(packed)

    if has_missing:
        da = DrawingArea(15, 10, 0, 0)
        rect = plt.Rectangle((0, 0), 15, 10, facecolor=COLOR_NA, edgecolor='black', lw=0.3)
        da.add_artist(rect)
        label = TextArea(" NA", textprops={'size': 8})
        packed = HPacker(children=[da, label], align="center", pad=0, sep=2)
        items.append(packed)

    block = VPacker(children=[label_title] + items, align="left", pad=0, sep=2)
    col_blocks.append(block)

legend_box = HPacker(children=col_blocks, align="top", pad=5, sep=20)
anchored_legend = AnchoredOffsetbox(
    loc='lower left', child=legend_box, pad=0., frameon=False,
    bbox_to_anchor=(0, -0.55), bbox_transform=ax.transAxes, borderpad=0.
)
ax.add_artist(anchored_legend)

plt.subplots_adjust(bottom=0.3)

# Save as PDF
fig.savefig('Figura1Ba.png', dpi=300, bbox_inches='tight')

plt.show()

# %% Figure 1Bb - Pie charts
# Dictionary with columns of interest
columns_to_plot = [
    'STATUS',
    'Time_Point_Ranges',
    'CAR_Type',
    'Max_Response',
    'ICANS_Grade_Range',
    'Age_Range',
    'Sex'
]

# Plot each column with its defined color
for col in columns_to_plot:
    value_counts = Df_heatmap[col].value_counts()
    print(value_counts)
    colors = [color_maps[col].get(cat, COLOR_NA) for cat in value_counts.index]  # grey color if NA

    plt.figure(figsize=(6,6))
    plt.pie(value_counts, autopct='%1.1f%%', startangle=140, colors=colors)
    plt.axis('equal')
    plt.tight_layout()

    # Save image
    filename = f"Figura1Bb_{col}_pie_chart.pdf"
    plt.savefig(filename, dpi=300)

    plt.show()

# %% End of script
