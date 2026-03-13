###############################################################################
###############################################################################
###############################################################################
###############################################################################

# Program: Supplementary_Table_10.py
# Author: Sergio Cámara Peña
# Date: 30/07/2025
# Version: FINAL

###############################################################################
###############################################################################
###############################################################################
###############################################################################

# %% Load all the needed libraries
import scanpy as sc
import os
import random
import seaborn as sns
import pandas as pd
from scipy.stats import fisher_exact

# %% Set a random seed
random.seed(2504)

# %% Read data
os.chdir("/home/scamara/data_a/scamara/Atlas/Input")
adata = sc.read("Python_scVI_adata_big_V4_state4.h5ad")

# %% Filter and mark IACs
adata2 = adata.copy()
adata2.obs['IACs'] = adata2.obs.groupby('Norm_Patient_Name')['manual_celltype_annotation_high'].transform(lambda x: 'Yes' if 'Monocyte-like T cells' in x.values else 'No')
print(adata2.shape)
adata2 = adata2[adata2.obs["Time_Point_Ranges"] == "Infusion_Product"]
adata2 = adata2[adata2.obs["STATUS"] == "DISEASE"]

print(adata2.shape)
del adata

# %% Contigency table
# Select the relevant columns
df = adata2.obs[["IACs", "ICANS_Grade_Range", "Norm_Patient_Name"]]

# Remove rows with NaN in 'ICANS_Grade_Range'
df = df.dropna(subset=['ICANS_Grade_Range'])

# Remove rows where 'ICANS_Grade_Range' is 0
df = df[~df['ICANS_Grade_Range'].isin([0, '0'])]

# Drop row names and keep only unique occurrences
df = df.drop_duplicates().reset_index(drop=True)

# Set 'Norm_Patient_Name' as the new index
df.set_index('Norm_Patient_Name', inplace=True)

# Create the contingency table
contingency_table = pd.crosstab(df["ICANS_Grade_Range"], df["IACs"])
print(contingency_table)

# %% Perform Fisher’s exact test
oddsratio, p_value = fisher_exact(contingency_table)

print("Odds ratio:", oddsratio)
print("p-value:", p_value)

# %% End of script