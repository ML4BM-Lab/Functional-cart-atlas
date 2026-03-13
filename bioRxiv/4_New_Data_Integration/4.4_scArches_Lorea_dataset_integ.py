###############################################################################
###############################################################################

# Program: 4.4_scArches_Lorea_dataset_integ.py
# Author: Sergio Cámara Peña
# Date: 23/01/2025
# Version: FINAL

###############################################################################
###############################################################################

# %% Link to tutorial followed:
# https://docs.scarches.org/en/latest/scanvi_surgery_pipeline.html

# %% Load all the needed libraries
import os
import scanpy as sc
import torch
import scarches as sca
from scarches.dataset.trvae.data_handling import remove_sparsity
import matplotlib.pyplot as plt
import numpy as np
import random
import pandas as pd
import anndata as ad

# %% Set a random seed
random.seed(2504)

# %% Load needed files
# Load ATLAS
os.chdir("/home/scamara/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Raw_Atlas")
adata_Reference = sc.read_h5ad("Python_scVI_adata_big_V4_state4.h5ad")

# Load Jordanas (Query)
os.chdir("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/Python-Celltypist/V5")
adata_V5 = sc.read_h5ad("Seurat_merged_With_Celltypist.h5ad")
adata_Query = adata_V5.copy()
adata_Query = adata_Query[adata_Query.obs["orig.ident"] == "Jordana_et_al"]

# %% Switches
Train = False
Train_2 = False
Save_h5ad_1 = False
Save_h5ad_2 = False
Save_h5ad_3 = False
Save_h5ad_4 = True

# %% Create a layer with the counts - Is what scVI uses
adata_Reference.layers["counts"] = adata_Reference.X.copy()

# %% Normalize and log scale
sc.pp.normalize_total(adata_Reference, target_sum=1e4)
sc.pp.log1p(adata_Reference)
adata_Reference.raw = adata_Reference  # Freeze the state in ".raw"

# %% Identify and only keep the 2000 most variable genes ----------- GIVES ERROR loess fit (I had modified span argument to 0.6, default is 0.3)
sc.pp.filter_genes(adata_Reference, min_cells=10)
sc.pp.highly_variable_genes(
    adata_Reference,
    flavor="seurat_v3",
    n_top_genes=2000,
    layer="counts",
    batch_key="Product_norm",
    span=0.6,
    subset=True,
)

# %% scVI setup
adata_Reference.obs.rename(columns={"orig.ident": "orig_ident"}, inplace=True)

if Train:
    sca.models.SCVI.setup_anndata(adata_Reference,
        layer="counts",
        batch_key="Product_norm",
        labels_key="manual_celltype_annotation_high")
    
    vae = sca.models.SCVI(
    adata_Reference,
    n_layers=2,
    #n_latent=30, # Go to default, better results
    encode_covariates=True,
    deeply_inject_covariates=False,
    use_layer_norm="both",
    use_batch_norm="none",
    gene_likelihood="nb"
)

    sca.models.SCVI.view_anndata_setup(vae)

# %% Train the model
if Train:
    vae.train()

# %% Save/Load trained model
os.chdir(
    "/home/scamara/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/scArches/V5"
)

if Train:
    vae.save("Atlas_integ_scArches_Reference_V5", overwrite=True, save_anndata=True)
    print("Saved successfully")
else:
    vae = sca.models.SCVI.load("Atlas_integ_scArches_Reference_V5")
    print("Loaded successfully")

# %% Create the SCANVI model instance
if Train:
    scanvae = sca.models.SCANVI.from_scvi_model(vae, unlabeled_category = "Unknown", gene_likelihood="nb")
    print("Labelled Indices: ", len(scanvae._labeled_indices))
    print("Unlabelled Indices: ", len(scanvae._unlabeled_indices))

# %% Train the model
if Train:
    scanvae.train(max_epochs=20, check_val_every_n_epoch=1, batch_size=256)

# %% Save/Load trained model
os.chdir(
    "/home/scamara/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/scArches/V5"
)

if Train:
    scanvae.save("Atlas_integ_scArches_Reference_V5_scanvae", overwrite=True, save_anndata=True)
    print("Saved successfully")
else:
    scanvae = sca.models.SCANVI.load("Atlas_integ_scArches_Reference_V5_scanvae")
    print("Loaded successfully")

# %% Create anndata file of latent representation and compute UMAP
if Save_h5ad_1:
    reference_latent = sc.AnnData(scanvae.get_latent_representation())
    reference_latent.obs["cell_type"] = adata_Reference.obs["manual_celltype_annotation_high"].tolist()
    reference_latent.obs["batch"] = adata_Reference.obs["Product_norm"].tolist()

    sc.pp.neighbors(reference_latent, n_neighbors=8)
    sc.tl.leiden(reference_latent)
    sc.tl.umap(reference_latent)
    sc.pl.umap(reference_latent,
               color=['cell_type'],
               frameon=False,
               wspace=0.6,
               )

# %% Save/Load adata object to state 1 h5ad - Reference
########################################################################
#
# State 1 Load/Save
# 
######################################################################## 
os.chdir("/home/scamara/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/scArches/V5")

if Save_h5ad_1:
    reference_latent.write("Atlas_integ_scArches_Reference_V5_state1.h5ad")
    print("Correctly saved")
else:
    reference_latent = sc.read_h5ad("Atlas_integ_scArches_Reference_V5_state1.h5ad")
    print("Correctly loaded")

# %% Compute the accuracy of the learned classifier
reference_latent.obs['predictions'] = scanvae.predict()
print("Acc: {}".format(np.mean(reference_latent.obs.predictions == reference_latent.obs.cell_type)))

#%% Prepare query
########################################################################
########################################################################
#                                                                      #
########################################################################
########################################################################

# %% Create a layer with the counts - Is what scVI uses
if Train_2:
    adata_Query.layers["counts"] = adata_Query.X.copy()

# %% Normalize and log scale
if Train_2:
    sc.pp.normalize_total(adata_Query, target_sum=1e4)
    sc.pp.log1p(adata_Query)
    adata_Query.raw = adata_Query  # Freeze the state in ".raw"

# %% Rename column of metadata
if Train_2:
    adata_Query.obs.rename(columns={"orig.ident": "orig_ident"}, inplace=True)

# %% Keep the same 2000 genes as the reference in the query dataset
if Train_2:
    adata_Query_copy = adata_Query.copy()
    genes_to_keep = adata_Reference.var_names.tolist()

    # Ensure genes_to_keep contains the 2000 differential genes, to avoid bugs, if not, error.
    if len(genes_to_keep) != 2000:
        raise ValueError(f"Expected 2000 genes, but found {len(genes_to_keep)}.")
    else:
        print("2000 genes, OK!")

    adata_Query = adata_Query[:, genes_to_keep]

# %% Perform surgery on reference model and train on query dataset without cell type labels
if Train_2:
    model = sca.models.SCANVI.load_query_data(
        adata_Query.copy(),
        reference_model="Atlas_integ_scArches_Reference_V5_scanvae",
        freeze_dropout = True,
    )

    model._unlabeled_indices = np.arange(adata_Query.n_obs)
    model._labeled_indices = []
    print("Labelled Indices: ", len(model._labeled_indices))
    print("Unlabelled Indices: ", len(model._unlabeled_indices))

# %% Train model
if Train_2:
    model.train(
        max_epochs=100,
        plan_kwargs=dict(weight_decay=0.0),
        check_val_every_n_epoch=10,
        batch_size=256
)

# %% Save/Load trained model
os.chdir(
    "/home/scamara/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/scArches/V5"
)

if Train_2:
    model.save("Atlas_integ_scArches_Query_V5_scanvae", overwrite=True, save_anndata=True)
    print("Saved successfully")
else:
    model = sca.models.SCANVI.load("Atlas_integ_scArches_Query_V5_scanvae")
    print("Loaded successfully")


# %% Create anndata file of latent representation and compute UMAP
query_latent = sc.AnnData(model.get_latent_representation())
query_latent.obs['cell_type'] = model.predict()
query_latent.obs['batch'] = adata_Query.obs["Product_norm"].tolist()
query_latent.obs["Jordanas_Original_Annotation"] = adata_Query.obs["Lorea_Annotation"].tolist()

sc.pp.neighbors(query_latent)
sc.tl.leiden(query_latent)
sc.tl.umap(query_latent)
plt.figure()
sc.pl.umap(
    query_latent,
    color=["batch", "cell_type"],
    frameon=False,
    wspace=0.6,
)

# %% Equivalence with Original Annotation
name_mapping = {
    "Proliferative": "Proliferative T cells",
    "CD8 effector memory": "CD8 effector memory",
    "CD4 activated": "CD4 effector memory",
    "CD4 memory": "CD4 central memory",
    "CD8 effector 1": "CD8 cytotoxic",
    "CD8 effector 2": "CD8 cytotoxic",
    "CD8 memory": "CD8 memory",
    "CD8 effectory GNLY-": "CD8 cytotoxic",
    "CD4 early memory": "CD4 central memory",
    "CD4 proliferative": "Proliferative T cells",
    "CD8 resting": "Unknown",
    "CD8 effector proliferative": "CD8 cytotoxic"
}

query_latent.obs["Jordanas_Equivalence"] = query_latent.obs["Jordanas_Original_Annotation"].replace(name_mapping)

# %% Save/Load adata object to state 1 h5ad - Query
########################################################################
#
# State 1 Load/Save
# 
######################################################################## 
os.chdir("/home/scamara/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/scArches/V5")

if Save_h5ad_2:
    query_latent.write("Atlas_integ_scArches_Query_V5_state1.h5ad")
    print("Correctly saved")
else:
    query_latent = sc.read_h5ad("Atlas_integ_scArches_Query_V5_state1.h5ad")
    print("Correctly loaded")


# %% Calculate accuracy
# Find common categories in both columns
common_labels = set(query_latent.obs["Jordanas_Equivalence"]) & set(query_latent.obs["cell_type"])

# Filter the dataframe to keep only rows where both columns have common values
filtered_obs = query_latent.obs[
    query_latent.obs["Jordanas_Equivalence"].isin(common_labels) &
    query_latent.obs["cell_type"].isin(common_labels)
]

# Compute accuracy on the filtered data
acc = np.mean(filtered_obs["Jordanas_Equivalence"].astype(str) == filtered_obs["cell_type"].astype(str))
print(f"Acc: {acc}")

# %% Calculate confusion matrix
# Create contingency table
df = filtered_obs.groupby(["Jordanas_Equivalence", "cell_type"]).size().unstack(fill_value=0)

# Drop empty rows (observed categories with no matching predictions)
df = df.loc[(df.sum(axis=1) > 0), :]

# Drop empty columns (predicted categories with no matching observations)
df = df.loc[:, (df.sum(axis=0) > 0)]

# Normalize (row-wise)
norm_df = df.div(df.sum(axis=1), axis=0).fillna(0)

# Plot
plt.figure(figsize=(8, 8))
plt.pcolor(norm_df, cmap="viridis")
plt.xticks(np.arange(0.5, len(df.columns), 1), df.columns, rotation=90)
plt.yticks(np.arange(0.5, len(df.index), 1), df.index)
plt.xlabel("Predicted")
plt.ylabel("Observed")
plt.colorbar(label="Proportion")

plt.show()

# %% Get latent representation of reference + query dataset - Fix missing values for the .obs
# Fix to have a similar/common metadata structure
# List of columns to keep
cols_to_keep = [
    'orig_ident', 'nCount_RNA', 'nFeature_RNA', 'Product',
    'log10GenesPerUMI', 'mitoRatio', 'S.Score', 'G2M.Score', 'Phase',
    'Product_norm', 'Marked_cells', 'barcode', 'timepoint',
    'barcode_timepoint', 'celltypist_predicted_labels',
    'celltypist_over_clustering', 'celltypist_majority_voting',
    'celltypist_conf_score', 'manual_celltype_annotation_high'
]

# Identify columns to drop
cols_to_drop = adata_Reference.obs.loc[:, 'celltypist_conf_score':].columns.difference(cols_to_keep)

# Drop the identified columns
adata_Reference.obs.drop(columns=cols_to_drop, inplace=True)

# Fix adata_Query
adata_Query.obs.rename(columns={"orig.ident": "orig_ident"}, inplace=True)
adata_Query.obs['manual_celltype_annotation_high'] = model.predict()
genes_to_keep = adata_Reference.var_names.tolist()
adata_Query = adata_Query[:, genes_to_keep]

for col in adata_Query.obs.select_dtypes(include=["object"]).columns:
    adata_Query.obs[col] = adata_Query.obs[col].astype("category")

# Create concantenated adata
adata_full = ad.concat( [adata_Reference,adata_Query])
adata_full.layers["counts"] = adata_full.X.copy()
full_latent = sc.AnnData(model.get_latent_representation(adata=adata_full))

# %% Compute UMAP and plot it
full_latent.obs['cell_type'] = adata_full.obs['manual_celltype_annotation_high'].tolist()
full_latent.obs['batch'] = adata_full.obs["Product_norm"].tolist()

sc.pp.neighbors(full_latent)
sc.tl.leiden(full_latent)
sc.tl.umap(full_latent)
plt.figure()
sc.pl.umap(
    full_latent,
    color=["cell_type"],
    frameon=False,
    wspace=0.6,
)

# %% Save/Load adata object to state 1 h5ad - Full
########################################################################
#
# State 1 Load/Save
# 
######################################################################## 
os.chdir("/home/scamara/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/scArches/V5")

if Save_h5ad_3:
    full_latent.write("Atlas_integ_scArches_Full_V5_state1.h5ad")
    print("Correctly saved")
else:
    full_latent = sc.read_h5ad("Atlas_integ_scArches_Full_V5_state1.h5ad")
    print("Correctly loaded")

# %% Load needed files
# Load ATLAS
os.chdir("/home/scamara/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Raw_Atlas")
adata_Reference = sc.read_h5ad("Python_scVI_adata_big_V4_state4.h5ad")

# Load Jordanas (Query)
os.chdir("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/Python-Celltypist/V5")
adata_V5 = sc.read_h5ad("Seurat_merged_With_Celltypist.h5ad")
adata_Query = adata_V5.copy()
adata_Query = adata_Query[adata_Query.obs["orig.ident"] == "Jordana_et_al"]

# Create concantenated adata
adata_big = ad.concat( [adata_Reference,adata_Query])

# %% Pass reductions to big object
adata_big.obs = adata_full.obs.copy()
adata_big.obs["leiden"] = full_latent.obs["leiden"].copy()
adata_big.uns = full_latent.uns.copy()
adata_big.obsm = full_latent.obsm.copy()
adata_big.obsp = full_latent.obsp.copy()

# %% Add final metadata
os.chdir("/home/scamara/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/scVI/V5")
sc_Metadata_to_python_v33 = pd.read_csv("sc_Metadata_to_python_v33_V5.csv", sep=";") # Check that this is the latest version of the metadata
sc_Metadata_to_python_v33.rename(columns={'Norm_Sample_Name': 'Product_norm'}, inplace=True)

# Merge metadata
metadata_bis = adata_big.obs.copy()
metadata_bis['row_number'] = range(1, len(metadata_bis) + 1)
print(metadata_bis['row_number'].is_monotonic_increasing)
metadata_bis['Product_norm'] = metadata_bis['Product_norm'].astype('object')
merged_df = pd.merge(metadata_bis, sc_Metadata_to_python_v33, on="Product_norm", how="right")
merged_df.index = adata_big.obs.index

# Check everything is passed correctly and delete the column
print(merged_df['row_number'].is_monotonic_increasing)
del merged_df['row_number']
del metadata_bis['row_number']

# Create the new object
adata_big_2 = adata_big.copy()
adata_big_2.obs = merged_df.copy()

# %% Save/Load final object
########################################################################
#
# FINAL object Load/Save
# 
######################################################################## 
os.chdir("/home/scamara/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/scArches/V5")

if Save_h5ad_4:
    adata_big_2.raw = None


    ##### Changes in metadata to up-to-date to latest version #####
    ## Load CSV with changes
    os.chdir("/home/scamara/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Raw_Atlas")
    time_point_changes = pd.read_csv("Time_Point_Changes_V5.csv", sep=";", index_col="Norm_Sample_Name")

    # Join directly to update Time_Point and Time_Point_Ranges
    adata_big_2.obs = adata_big_2.obs.drop(columns=["Time_Point", "Time_Point_Ranges"], errors="ignore")
    adata_big_2.obs = adata_big_2.obs.join(
        time_point_changes,
        on="Product_norm"   # match Product_norm in obs with Norm_Sample_Name index
    )

    # Quick check
    pd.set_option("display.max_rows", None)
    Check = adata_big_2.copy()
    print(Check.obs[["Product_norm", "Time_Point", "Time_Point_Ranges"]].drop_duplicates().to_string(index=False))


    adata_big_2.write("Atlas_integ_scArches_FINAL_V5.h5ad")
    print("Correctly saved")
else:
    adata_big_2 = sc.read_h5ad("Atlas_integ_scArches_FINAL_V5.h5ad")
    print("Correctly loaded")

# %% Check umap
sc.pl.umap(
    adata_big_2,
    color=["manual_celltype_annotation_high"],
    frameon=False,
    wspace=0.6,
)

# %% End of script