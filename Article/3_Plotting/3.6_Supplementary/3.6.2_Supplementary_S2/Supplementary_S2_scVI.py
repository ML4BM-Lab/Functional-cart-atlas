###############################################################################
###############################################################################

# Program: Supplementary_S2_scVI.py
# Author: Sergio Cámara Peña
# Date: 09/06/2025
# Version: FINAL

###############################################################################
###############################################################################

##### Use of scVI in its Python version #####
## More info: https://ccbskillssem.github.io/assets/scvi_notebook.html y https://docs.scvi-tools.org/en/stable/tutorials/notebooks/harmonization.html

#%% Load all the needed libraries
import scanpy as sc
import scvi
from scvi.model.utils import mde
import scanpy.external as sce
import os
import random
import pandas as pd
import seaborn as sns

#%% Set a random seed
random.seed(2504)

# %% Define a 12 color palette
palette = [
"#F8766D",
"#00BA38",
"#FF9900",
"#619CFF",
"#F564E3",
"#B79F00"
]

# %% Switches
Primera_vez = False

#%% Read the data
adata = sc.read("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration_methods_lab/scVI/Sin_GT_With_Python/Seurat_merged_RAW_for_Py.h5ad") # Obtained in Supplementary_S2_Methods.R script
adata.raw = adata  # keep full dimension safe

# %% Obtain the cellular annotation
os.chdir("/home/scamara/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Raw_Atlas")
adata_orig = sc.read_h5ad("Python_scVI_adata_big_V4_state4.h5ad")

adata_orig.obs.index = adata_orig.obs.index.str.replace(r"_\d+$", "", regex=True)

common_cells = adata.obs_names.intersection(adata_orig.obs_names)
print(f"Common cells: {len(common_cells)} of {adata.n_obs}")

adata_subset = adata[common_cells].copy()

adata_subset.obs["manual_celltype_annotation_high"] = adata_orig.obs.loc[common_cells, "manual_celltype_annotation_high"]

os.chdir("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration_methods_lab")
df_annotation = adata_subset.obs[["manual_celltype_annotation_high"]].copy()
df_annotation["cell_name"] = df_annotation.index
df_annotation = df_annotation.reset_index(drop=True)
df_annotation.to_csv("cell_annotation_manual.csv", index=False)

adata_subset = adata[common_cells].copy()
adata_subset.obs["manual_celltype_annotation_high"] = adata_orig.obs.loc[common_cells, "manual_celltype_annotation_high"]
del adata

adata = adata_subset.copy()

#%% Read the metadata, process and adapt it
df = pd.read_csv("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Datasets_Metadata/scRNAseq_metadata_CARTs_v10.csv", skipfooter=6, sep=";")
Data_rows = adata.obs.Product.unique()
filtered_df = df[df["Sample name"].isin(Data_rows)]
filtered_df = filtered_df[["Sample name", "Age", "Age_Range", "Sex", "CAR_Construct", "CAR_Gen", "ScFv", "Costim_Domain_1", "Costim_Domain_2", "STATUS", "Stimulated"]]

#%% Normalize the data
adata.layers["counts"] = adata.X.copy() # preserve counts
sc.pp.normalize_total(adata, target_sum=1e4) # scale each cell to a common library size
sc.pp.log1p(adata) # log(expression + 1)

#%% Identify and only keep the 2000 most variable genes
sc.pp.highly_variable_genes(
    adata,
    flavor="seurat_v3",
    n_top_genes=2000,
    batch_key="Product_norm",
    subset=True,
    layer="counts"
)

#%% Standard Workflow
sc.pp.scale(adata) # Normalize the columns (genes)
sc.tl.pca(adata)

adata.obsm["X_pca"]

sc.pp.neighbors(adata) # Compute nearest neighbors
sc.tl.umap(adata)

sc.pl.umap(adata, color="Product")

sce.pp.bbknn(adata, batch_key="Product")
sc.tl.umap(adata)

#%% scVI setup
if Primera_vez:
    scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="Product")
    model = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
    scvi.model.SCVI.view_anndata_setup(model)

#%% Train the model
if Primera_vez:
    model.train()

#%% Save/Load trained model
os.chdir("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration_methods_lab/scVI/Sin_GT_With_Python")

if Primera_vez:
    model.save("Healty_donors_scvi_v1", overwrite=True, save_anndata=True)
else:
    model = scvi.model.SCVI.load('Healty_donors_scvi_v1', adata)

#%% Obtain latent representation from scVI to evaluate it
adata.obsm["X_scVI"] = model.get_latent_representation()

# %% Neighbours, leiden (to cluster cells) and UMAP calculation
sc.pp.neighbors(adata, use_rep="X_scVI")
sc.tl.leiden(adata, key_added="leiden_scvi", resolution=1.2)
sc.tl.umap(adata, min_dist=0.4)

# %% Save or read .h5ad
Editar_Figura = False
os.chdir("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration_methods_lab/scVI/Sin_GT_With_Python")
if Editar_Figura:
    adata = sc.read_h5ad("Suplementaria_S2_scVI_state1.h5ad")
    print("Correctly loaded")
else:
    adata.write("Suplementaria_S2_scVI_state1.h5ad")
    print("Correctly saved")

#%% UMAP representation and save
os.chdir("/home/scamara/data/scamara/Atlas_Mieloma_Multiple/Resultados_Figuras/Suplementarias/")

sns.palplot(palette)

sc.pl.umap(
    adata,
    color=["orig.ident"],
    frameon=False,
    size=4,
    palette=palette,
    save="Suplementaria_S2_2_scVI.pdf"
)

# %% End of script
