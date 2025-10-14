###############################################################################
###############################################################################

# Program: 2.5_Subclustering_V4_Reunite.py
# Author: Sergio Cámara Peña
# Date: 22/02/2024
# Version: FINAL

###############################################################################
###############################################################################

# %% Load all the needed libraries
import scanpy as sc
import os
import random
import pandas as pd
import seaborn as sns
from scanpy_cluster_proportions import get_cluster_proportions, plot_cluster_proportions

# %% Set a random seed
random.seed(2504)

# %% Read files
os.chdir("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/scVI/V4")
adata = sc.read_h5ad("Python_scVI_adata_V4_state3.h5ad")

os.chdir("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/Subclustering_V4/Data/CD4")
CD4_annotation_high = pd.read_csv("adata_CD4_manual_celltype_annotation_high.csv", index_col=0)

os.chdir("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/Subclustering_V4/Data/CD8")
CD8_annotation_high = pd.read_csv("adata_CD8_manual_celltype_annotation_high.csv", index_col=0)

os.chdir("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/Subclustering_V4/Data/GATA3")
GATA3_annotation_high = pd.read_csv("adata_GATA3_manual_celltype_annotation_high.csv", index_col=0)

# %% Add annotation info to general table
metadata = adata.obs.copy()
anno_info = pd.concat([CD4_annotation_high, CD8_annotation_high, GATA3_annotation_high])

merged_df = metadata.join(anno_info)

merged_df["manual_celltype_annotation_high"] = merged_df["manual_celltype_annotation_high"].fillna(merged_df["manual_celltype_annotation_low"])

adata.obs = merged_df.copy()

# %% Plot result
os.chdir("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/Plots/V4/Subclustering/Final")
sc.pl.umap(adata, color = ["manual_celltype_annotation_high"], frameon=False, palette=['#d21820', '#1869ff', '#008a00', '#f36dff', '#710079', '#aafb00', '#00bec2', '#ffa235', '#5d3d04', '#08008a', '#005d5d', '#9a7d82'], save="Atlas_V4_manual_celltype_annotation_high.png")
sc.pl.umap(adata, color = ["manual_celltype_annotation_high"], frameon=False, palette=['#d21820', '#1869ff', '#008a00', '#f36dff', '#710079', '#aafb00', '#00bec2', '#ffa235', '#5d3d04', '#08008a', '#005d5d', '#9a7d82'], save="Atlas_V4_manual_celltype_annotation_high.pdf")
sc.pl.umap(adata, color = ["manual_celltype_annotation_high"], frameon=False, palette=['#d21820', '#1869ff', '#008a00', '#f36dff', '#710079', '#aafb00', '#00bec2', '#ffa235', '#5d3d04', '#08008a', '#005d5d', '#9a7d82'], save="Atlas_V4_manual_celltype_annotation_high.svg")

# %% Atlas colorized based on metadata -- Generation of datasets
os.chdir("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/scVI/V4")
sc_Metadata_to_python_v33 = pd.read_csv("sc_Metadata_to_python_v33_V4.csv", sep=";") # Check is the latest version of the metadata
sc_Metadata_to_python_v33.rename(columns={'Norm_Sample_Name': 'Product_norm'}, inplace=True)
metadata_bis = adata.obs.copy()
metadata_bis['row_number'] = range(1, len(metadata_bis) + 1)
metadata_bis['row_number'].is_monotonic_increasing
metadata_bis['Product_norm'] = metadata_bis['Product_norm'].astype('object')
merged_df = pd.merge(metadata_bis, sc_Metadata_to_python_v33, on="Product_norm", how="right")
merged_df.index = adata.obs.index

# %% Check everything is passed correctly and delete the column
merged_df['row_number'].is_monotonic_increasing
del merged_df['row_number']
del metadata_bis['row_number']

# %% Check everything is passed correctly 2
Check_1 = metadata_bis
Check_2 = merged_df.iloc[:, :39]

print(Check_2.equals(Check_1))

from pandas.testing import assert_frame_equal

try:
    assert_frame_equal(Check_1, Check_2)
    print("DataFrames are the same.")
except AssertionError as e:
    print("DataFrames are not the same.")
    print(e)

# %% Check everything is passed correctly 3
Check_3 = sc_Metadata_to_python_v33
Check_4 = merged_df.iloc[:, [9] + list(range(39, merged_df.shape[1]))].drop_duplicates().reset_index(drop=True)

print(Check_4.equals(Check_3))

from pandas.testing import assert_frame_equal

try:
    assert_frame_equal(Check_3, Check_4)
    print("DataFrames are the same.")
except AssertionError as e:
    print("DataFrames are not the same.")
    print(e)

# %% Check names of rownames and metadata is correct
cell_names = merged_df.index.copy()

import re

# Remove the trailing number (if present) and then remove the sequence part
def clean_name(name):
    # Step 1: Remove the trailing number (if it exists)
    name = re.sub(r'_\d+$', '', name)
    
    # Step 2: Remove everything after the last underscore (the sequence part)
    name = "_".join(name.split("_")[:-1])
    
    return name

processed_names = pd.Series([clean_name(name) for name in cell_names])
Check_5 = processed_names.copy()

Check_6 = merged_df["Product_norm"].copy()

Check_5_reset = Check_5.reset_index(drop=True)
Check_6_reset = Check_6.reset_index(drop=True)

if (Check_5_reset.values != Check_6_reset.values).any():
    # Find the discrepancies
    discrepancies = Check_5_reset[Check_5_reset != Check_6_reset]
    
    # Raise an error with the discrepancies
    raise ValueError(f"Discrepancies found in the following rows:\n{discrepancies}")
else:
    print("The two Series are the same.")

# %% Create the new object
adata2 = adata.copy()
adata2.obs = merged_df.copy()

# %% Save results to state 4
os.chdir("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/scVI/V4")
adata2.write("Python_scVI_adata_V4_state4.h5ad")

# %% End of script
