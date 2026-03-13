###############################################################################
###############################################################################

# Program: 2.3_Integrated_datasets_Workflow.py
# Author: Sergio Cámara Peña
# Date: 08/08/2023
# Version: FINAL

###############################################################################
###############################################################################

##### Use scVI in its Python version #####
## Performed following these two tutorials: https://ccbskillssem.github.io/assets/scvi_notebook.html y https://docs.scvi-tools.org/en/stable/tutorials/notebooks/harmonization.html
## Youtube video with more info: https://www.youtube.com/watch?v=YT9qTuF6YFk&ab_channel=Sanbomics
## Colab notebook: https://colab.research.google.com/github/yoseflab/scvi_tutorials/blob/0.14.5/api_overview.ipynb#scrollTo=QXIxM_VgqovL

##### More info: #####
## https://www.kaggle.com/code/hiramcho/scrna-seq-differential-expression-with-scvi
## https://docs.scvi-tools.org/en/1.0.1/api/reference/scvi.model.SCVI.html#scvi.model.SCVI
## https://github.com/scverse/scanpy/issues/670 -- Clustree for python


# %% Load all the needed libraries
import scanpy as sc
import scvi
import os
import random
from scanpy_cluster_proportions import get_cluster_proportions, plot_cluster_proportions
import seaborn as sns
import pandas as pd

# %% Switches
Train = False
Save_h5ad = False
Save_h5ad_2 = False
Save_h5ad_2_5 = False
Save_h5ad_3 = False

# %% Set a random seed
random.seed(2504)

# %% Read the data
if Train:
    adata = sc.read(
        "/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/Python-Celltypist/V4/Seurat_merged_With_Celltypist.h5ad"
    )

# %% Create a layer with the counts - Is what scVI uses
if Train:
    adata.layers["counts"] = adata.X.copy()

# %% Normalize and log scale
if Train:
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata  # Freeze the state in ".raw"

# %% Remove samples with fewer than 45 CD3⁺ and CAR⁺ cells — since they cause issues in the next step
if Train:
    sample_column = 'Product_norm'

    # Count the number of cells per sample
    cell_counts = adata.obs[sample_column].value_counts(sort=False)

    # Filter samples with less than 45 cells
    filtered_samples = cell_counts[cell_counts >= 45].index

    # Subset the AnnData object to keep only samples with at least 45 cells
    adata_filtered = adata.copy()
    adata_filtered = adata_filtered[adata_filtered.obs[sample_column].isin(filtered_samples)]

    del adata
    adata = adata_filtered.copy()
    del adata_filtered

# %% Identify and only keep the 2000 most variable genes ----------- GIVES ERROR loess fit (I had modified span argument to 0.6, default is 0.3)
if Train:
    sc.pp.filter_genes(adata, min_cells=10)

    sc.pp.highly_variable_genes(
        adata,
        flavor="seurat_v3",
        n_top_genes=2000,
        layer="counts",
        batch_key="Product_norm",
        span=0.6,
        subset=True,
    )

# %% scVI setup
if Train:
    adata.obs.rename(columns={"orig.ident": "orig_ident"}, inplace=True)
    scvi.model.SCVI.setup_anndata(
        adata,
        layer="counts",
        batch_key="Product_norm",
        categorical_covariate_keys=["orig_ident"],
    )
    model = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
    scvi.model.SCVI.view_anndata_setup(model)

# %% Train the model
if Train:
    model.train()

# %% Save/Load trained model
os.chdir(
    "/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/scVI/V4"
)

if Train:
    model.save("Atlas_integ_scVI_V4", overwrite=True, save_anndata=False)
    print("Saved successfully")
else:
    model = scvi.model.SCVI.load("Atlas_integ_scVI_V4")
    print("Loaded successfully")

# %% Obtain latent representation from scVI to evaluate it
adata.obsm["X_scVI"] = model.get_latent_representation()

# %% Neighbours, leiden (to cluster cells) and UMAP calculation
sc.pp.neighbors(adata, use_rep="X_scVI")
sc.tl.umap(adata)
sc.tl.leiden(adata, key_added="leiden_Res_0.3", resolution=0.3)
sc.tl.leiden(adata, key_added="leiden_Res_0.4", resolution=0.4)
sc.tl.leiden(adata, key_added="leiden_Res_0.5", resolution=0.5)
sc.tl.leiden(adata, key_added="leiden_Res_0.6", resolution=0.6)

# %% Save/Load adata object to state 1 h5ad
########################################################################
#
# State 1 Load/Save
# 
######################################################################## 
if Save_h5ad:
    os.chdir(
        "/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/scVI/V4"
    )
    adata.write("Python_scVI_adata_V4_state1.h5ad")
    print("Correctly saved")
else:
    os.chdir("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/scVI/V4")
    adata = sc.read_h5ad("Python_scVI_adata_V4_state1.h5ad")
    print("Correctly loaded")

# %% I paused to check the resolutions using Clustree
if False:
    import invoke
    import os
    os.chdir("/home/scamara/data/scamara/Atlas_Mieloma_Multiple/Codigo/Codigo_datasets_atlas/Datasets_Integration/Tests_labs")
    CLUSTREE_SCRIPT = "Leiden_resol_lab_clustree_function.R"
    ADATA_PATH = "/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/scVI/V4/Python_scVI_adata_V4_state1.h5ad"
    CLUSTREE_OUT = "/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/scVI/V4/Clustree_resolutions/clustree.png"
    Carrera = invoke.run(f"Rscript {CLUSTREE_SCRIPT} {ADATA_PATH} {CLUSTREE_OUT}")
# After examining the resulting graph, I decided to use resolution 0.3

# %% UMAP plots save
os.chdir(
    "/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/Plots/V4/scVI_integ"
)
Features_to_explore = [
    "celltypist_majority_voting",
    "leiden_Res_0.3",
    "orig_ident",
]

for i in range(len(Features_to_explore)):
    feature_name = Features_to_explore[i]
    filename = f"_{feature_name}_scVI_V4.pdf"
    filename2 = f"_{feature_name}_scVI_V4.png"

    sc.pl.umap(adata, color=feature_name, frameon=False, save=filename)

    sc.pl.umap(adata, color=feature_name, frameon=False, save=filename2)

# %% Clustering QC Plots --- 4 QC parameters, dataset contribution to cluster, cluster composition by dataset, cell cycle
os.chdir("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/Plots/V4/scVI_integ/figures")

cluster_props = get_cluster_proportions(adata, cluster_key="leiden_Res_0.3", sample_key="orig_ident")
f1 = plot_cluster_proportions(cluster_props, xlabel_rotation=90, cluster_palette=sns.color_palette("hls", adata.obs["leiden_Res_0.3"].nunique()))
f1.savefig("QC_Datasets_to_Cluster_distrib.pdf")
f1.savefig("QC_Datasets_to_Cluster_distrib.png")

cluster_props_2 = get_cluster_proportions(adata, cluster_key="orig_ident", sample_key="leiden_Res_0.3")
f2 = plot_cluster_proportions(cluster_props_2, xlabel_rotation=90, cluster_palette=sns.color_palette("hls", adata.obs["orig_ident"].nunique()))
f2.savefig("QC_Clusters_to_Dataset_distrib.pdf")
f2.savefig("QC_Clusters_to_Dataset_distrib.png")

os.chdir("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/Plots/V4/scVI_integ")
sc.pl.umap(adata, color="Phase", frameon=False, save="_QC_Cell_Cycle.pdf")
sc.pl.umap(adata, color="Phase", frameon=False, save="_QC_Cell_Cycle.png")

sc.pl.umap(adata, color=["nCount_RNA", "nFeature_RNA", "log10GenesPerUMI", "mitoRatio"], frameon=False, save="_QC_4params.pdf")
sc.pl.umap(adata, color=["nCount_RNA", "nFeature_RNA", "log10GenesPerUMI", "mitoRatio"], frameon=False, save="_QC_4params.png")

# %% Clustering QC Data --- How many cells are in each cluster? How many cell are in every cell type asigned by celltypist?
metadata = adata.obs.copy()
cluster_counts = metadata["leiden_Res_0.3"].value_counts()
cluster_counts_df = pd.DataFrame({'Counts': cluster_counts.values})
print(cluster_counts_df)

celltypist_count = pd.DataFrame(metadata["celltypist_majority_voting"].value_counts())

os.chdir("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/scVI/V4/Data")
cluster_counts_df.to_csv("Counts_per_cluster.csv")
celltypist_count.to_csv("Celltypist_counts.csv")

# %% Finding marker genes by cluster
sc.pp.log1p(adata) # Calling this function before ranking genes avoids the need to use the raw data and fixes the dictionary issue
sc.tl.rank_genes_groups(adata, "leiden_Res_0.3", method="wilcoxon") # This function expects the data to be log-transformed

os.chdir("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/Plots/V4/scVI_integ")
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, save=".pdf")
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, save=".png")

sc.pl.rank_genes_groups_stacked_violin(adata, n_genes=4, cmap='RdBu_r', groupby="leiden_Res_0.3", save="scvi.pdf")

# Example to compare specific cluster (https://nbisweden.github.io/workshop-archive/workshop-scRNAseq/2020-01-27/labs/compiled/scanpy/scanpy_05_dge.html)
# sc.tl.rank_genes_groups(adata, 'louvain_0.6', groups=['2'], reference='6', method='wilcoxon')
# sc.pl.rank_genes_groups(adata, groups=['2'], n_genes=20)

# %% Marker genes by cluster --- Table
os.chdir("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/scVI/V4/Data")
result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names

Marker_genes_x_clus_df = pd.DataFrame(
    {group + '_' + key[:1]: result[key][group]
    for group in groups for key in ['names', 'pvals']})

Marker_genes_x_clus_df
Marker_genes_x_clus_df.to_csv("Marker_genes_per_cluster.csv")

#%% Load/Save state 2
########################################################################
#
# State 2 Load/Save
# 
######################################################################## 
if Save_h5ad_2:
    os.chdir(
        "/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/scVI/V4"
    )
    adata.write("Python_scVI_adata_V4_state2.h5ad")
    print("Correctly saved")
else:
    os.chdir("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/scVI/V4")
    adata = sc.read_h5ad("Python_scVI_adata_V4_state2.h5ad")
    print("Correctly loaded")

# %% Differential marker genes between clusters --- Plots
os.chdir("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/Plots/V4/scVI_integ")

sc.pl.umap(adata, color=["CD4", "CD3D", "CD8A"], frameon=False, ncols=3, save="_marker_genes1.pdf") #ncols = nº of panels per row
sc.pl.umap(adata, color=["CD4", "CD3D", "CD8A"], frameon=False, ncols=3, save="_marker_genes1.png")
sc.pl.umap(adata, color=["celltypist_majority_voting", "celltypist_conf_score"], frameon=False, ncols=2, save="_celltypist_confidence.pdf")
sc.pl.umap(adata, color=["celltypist_majority_voting", "celltypist_conf_score"], frameon=False, ncols=2, save="_celltypist_confidence.png")

#%% Cell type UMAPs
os.chdir("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/Plots/V4/scVI_markers_exploration")

sc.pl.umap(adata, color=["CD4", "CD3D", "CD8A"], frameon=False, ncols=3, vmin=0, vmax="p99", save="0_Basic_Markers.png") # Basic Markers
sc.pl.umap(adata, color=["CD8A", "CD8B"], frameon=False, ncols=3, vmin=0, vmax="p99", save="0.1_Basic_Markers_CD8.png") # Basic Markers (CD8+)
sc.pl.umap(adata, color=["CD19", "MS4A1", "CD79A"], frameon=False, ncols=2, vmin=0, vmax="p99", save="1_B_Cell.png") # B cells
sc.pl.umap(adata, color=["CD4", "IL7R", "CCR7"], frameon=False, ncols=2, vmin=0, vmax="p99", save="2_T_CD4_Cell.png") # T CD4 cells
sc.pl.umap(adata, color=["FOXP3", "IL2RA", "CTLA4", "CCR4"], frameon=False, ncols=2, vmin=0, vmax="p99", save="3_Treg.png") # T reg cells (CD4+)
sc.pl.umap(adata, color=["TBX21", "IFNG", "CXCR3"], frameon=False, ncols=2, vmin=0, vmax="p99", save="4_Th1.png") # Th1 cells (CD4+)
sc.pl.umap(adata, color=["GATA3", "IL4"], frameon=False, ncols=2, vmin=0, vmax="p99", save="5_Th2.png") # Th2 cells (CD4+)
sc.pl.umap(adata, color=["RORC", "IL17A", "CCR6"], frameon=False, ncols=2, vmin=0, vmax="p99", save="6_Th17.png") # Th17 cells (CD4+)
sc.pl.umap(adata, color=["IL9", "CCR3"], frameon=False, ncols=2, vmin=0, vmax="p99", save="7_Th9.png") # Th9 cells (CD4+)
sc.pl.umap(adata, color=["BCL6", "CXCR5", "IL21", "CD28"], frameon=False, ncols=2, vmin=0, vmax="p99", save="8_Tfh.png") # Tfh cells (CD4+)
sc.pl.umap(adata, color=["AHR", "IL22"], frameon=False, ncols=2, vmin=0, vmax="p99", save="9_Th22.png") # Th22 cells (CD4+)
sc.pl.umap(adata, color=["GNLY", "NKG7", "NCAM1", "KLRD1"], frameon=False, ncols=2, vmin=0, vmax="p99", save="10_celula_NK.png") # NK cells
sc.pl.umap(adata, color=["CD14", "HLA-DRA", "FCGR3A"], frameon=False, ncols=2, vmin=0, vmax="p99", save="11_Monocitos.png") # Monocytes
sc.pl.umap(adata, color=["CD14", "LYZ"], frameon=False, ncols=2, vmin=0, vmax="p99", save="12_Monocitos_CD14.png") # CD14+ monocytes
sc.pl.umap(adata, color=["FCGR3A", "MS4A7"], frameon=False, ncols=2, vmin=0, vmax="p99", save="13_Monocitos_FCGR3A.png") # FCGR3A+ monocytes
sc.pl.umap(adata, color=["FCER1A", "CST3", "ITGAX", "HLA-DRA"], frameon=False, ncols=2, vmin=0, vmax="p99", save="14_Conventional_DC.png") # Conventional dendritic cell
sc.pl.umap(adata, color=["IL3RA", "HLA-DRA", "GZMB", "SERPINF1"], frameon=False, ncols=2, vmin=0, vmax="p99", save="15_Plasmacytoid_DC.png") # Plasmacytoid dendritic cells
sc.pl.umap(adata, color=["PPBP", "ITGA2B"], frameon=False, ncols=2, vmin=0, vmax="p99", save="16_Megak_and_platelts.png") # Megakaryocytes and Platelets
sc.pl.umap(adata, color=["HBB", "HBA2", "GYPA"], frameon=False, ncols=2, vmin=0, vmax="p99", save="17_Eritrocytes.png") # Erythrocytes

##### Functional classification
sc.pl.umap(adata, color=["PRF1", "IFNG", "GZMK", "GZMB", "CD8A", "NKG7"], frameon=False, ncols=2, vmin=0, vmax="p99", save="18_Cytotoxic_function.png") # Cytotoxic (effector) function
sc.pl.umap(adata, color=["GNLY", "GZMA", "GZMM", "CCL5", "CCL3", "CCL4"], frameon=False, ncols=2, vmin=0, vmax="p99", save="18.2_Cytotoxic_function.png") # Cytotoxic (effector) function 2
sc.pl.umap(adata, color=["CD69", "CD27", "CD28", "IL2RA"], frameon=False, ncols=2, vmin=0, vmax="p99", save="19_Activation_function.png") # Activation markers
sc.pl.umap(adata, color=["IL10", "IL32"], frameon=False, ncols=2, vmin=0, vmax="p99", save="20_Cytokines_secretion.png") # Cytokine secretion
sc.pl.umap(adata, color=["TCF7", "CCR7", "SELL", "IL7R"], frameon=False, ncols=2, vmin=0, vmax="p99", save="21_Memory.png") # Memory / Memory-like
sc.pl.umap(adata, color=["ZWINT", "MKI67", "TUBA1B", "TUBB"], frameon=False, ncols=2, vmin=0, vmax="p99", save="22_Proliferation.png") # Proliferation
sc.pl.umap(adata, color=["HAVCR2", "TIGIT", "LAG3", "CTLA4"], frameon=False, ncols=2, vmin=0, vmax="p99", save="23_Exhaustion.png") # Exhaustion
sc.pl.umap(adata, color=["KLRG1", "B3GAT1", "TIGIT", "CD160"], frameon=False, ncols=2, vmin=0, vmax="p99", save="24_Senescense.png") # Senescence
sc.pl.umap(adata, color=["SELL", "CXCR3", "CCR7", "CD44"], frameon=False, ncols=2, vmin=0, vmax="p99", save="25_Other_Memory_markers.png") # Other memory markers
sc.pl.umap(adata, color=["TBX21", "B3GAT1", "SELL", "CD44"], frameon=False, ncols=2, vmin=0, vmax="p99", save="25.2_Other_Memory_markers.png") # Other memory markers 2
sc.pl.umap(adata, color=["CD8A", "CD8B"], frameon=False, ncols=2, vmin=0, vmax="p99", save="26_CD8_Complex.png") # CD8 complex marker genes
sc.pl.umap(adata, color=["CD3D", "CD3G", "CD3E"], frameon=False, ncols=2, vmin=0, vmax="p99", save="27_CD3_Complex.png") # CD3 complex marker genes

## Extra
sc.pl.umap(adata, color=["TNFRSF17", "SDC1", "CD19"], frameon=False, ncols=2) # Study of residual cancer cells — It appears there are no residual cancer cells (TNFRSF17 = BCMA)

#%% Cell type violins
os.chdir("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/Plots/V4/scVI_markers_exploration")

sc.pl.violin(adata, keys=["CD4", "CD3D", "CD8A"], groupby='leiden_Res_0.3', stripplot=False, save="0_Basic_Markers.png") # Basic Markers
sc.pl.violin(adata, keys=["CD19", "MS4A1", "CD79A"], groupby='leiden_Res_0.3', stripplot=False, save="1_B_Cell.png") # B cells
sc.pl.violin(adata, keys=["CD4", "IL7R", "CCR7"], groupby='leiden_Res_0.3', stripplot=False, save="2_T_CD4_Cell.png") # T CD4 cells
sc.pl.violin(adata, keys=["FOXP3", "IL2RA", "CTLA4", "CCR4"], groupby='leiden_Res_0.3', stripplot=False, save="3_Treg.png") # T reg cells (CD4+)
sc.pl.violin(adata, keys=["TBX21", "IFNG", "CXCR3"], groupby='leiden_Res_0.3', stripplot=False, save="4_Th1.png") # Th1 cells (CD4+)
sc.pl.violin(adata, keys=["GATA3", "IL4"], groupby='leiden_Res_0.3', stripplot=False, save="5_Th2.png") # Th2 cells (CD4+)
sc.pl.violin(adata, keys=["RORC", "IL17A", "CCR6"], groupby='leiden_Res_0.3', stripplot=False, save="6_Th17.png") # Th17 cells (CD4+)
sc.pl.violin(adata, keys=["IL9", "CCR3"], groupby='leiden_Res_0.3', stripplot=False, save="7_Th9.png") # Th9 cells (CD4+)
sc.pl.violin(adata, keys=["BCL6", "CXCR5", "IL21", "CD28"], groupby='leiden_Res_0.3', stripplot=False, save="8_Tfh.png") # Tfh cells (CD4+)
sc.pl.violin(adata, keys=["AHR", "IL22"], groupby='leiden_Res_0.3', stripplot=False, save="9_Th22.png") # Th22 cells (CD4+)
sc.pl.violin(adata, keys=["GNLY", "NKG7", "NCAM1", "KLRD1"], groupby='leiden_Res_0.3', stripplot=False, save="10_celula_NK.png") # NK cells
sc.pl.violin(adata, keys=["CD14", "HLA-DRA", "FCGR3A"], groupby='leiden_Res_0.3', stripplot=False, save="11_Monocitos.png") # Monocytes
sc.pl.violin(adata, keys=["CD14", "LYZ"], groupby='leiden_Res_0.3', stripplot=False, save="12_Monocitos_CD14.png") # CD14+ monocytes
sc.pl.violin(adata, keys=["FCGR3A", "MS4A7"], groupby='leiden_Res_0.3', stripplot=False, save="13_Monocitos_FCGR3A.png") # FCGR3A+ monocytes
sc.pl.violin(adata, keys=["FCER1A", "CST3", "ITGAX", "HLA-DRA"], groupby='leiden_Res_0.3', stripplot=False, save="14_Conventional_DC.png") # Conventional dendritic cell
sc.pl.violin(adata, keys=["IL3RA", "HLA-DRA", "GZMB", "SERPINF1"], groupby='leiden_Res_0.3', stripplot=False, save="15_Plasmacytoid_DC.png") # Plasmacytoid dendritic cells
sc.pl.violin(adata, keys=["PPBP", "ITGA2B"], groupby='leiden_Res_0.3', stripplot=False, save="16_Megak_and_platelts.png") # Megakaryocytes and Platelets
sc.pl.violin(adata, keys=["HBB", "HBA2", "GYPA"], groupby='leiden_Res_0.3', stripplot=False, save="17_Eritrocytes.png") # Erythrocytes

##### Functional classification
sc.pl.violin(adata, keys=["PRF1", "IFNG", "GZMK", "GZMB", "CD8A", "NKG7"], groupby='leiden_Res_0.3', stripplot=False, save="18_Cytotoxic_function.png") # Cytotoxic (effector) function
sc.pl.violin(adata, keys=["GNLY", "GZMA", "GZMM", "CCL5", "CCL3", "CCL4"], groupby='leiden_Res_0.3', stripplot=False, save="18.2_Cytotoxic_function.png") # Cytotoxic (effector) function 2
sc.pl.violin(adata, keys=["CD69", "CD27", "CD28", "IL2RA"], groupby='leiden_Res_0.3', stripplot=False, save="19_Activation_function.png") # Activation markers
sc.pl.violin(adata, keys=["IL10", "IL32"], groupby='leiden_Res_0.3', stripplot=False, save="20_Cytokines_secretion.png") # Cytokine secretion
sc.pl.violin(adata, keys=["TCF7", "CCR7", "SELL", "IL7R"], groupby='leiden_Res_0.3', stripplot=False, save="21_Memory.png") # Memory / Memory-like
sc.pl.violin(adata, keys=["ZWINT", "MKI67", "TUBA1B", "TUBB"], groupby='leiden_Res_0.3', stripplot=False, save="22_Proliferation.png") # Proliferation
sc.pl.violin(adata, keys=["HAVCR2", "TIGIT", "LAG3", "CTLA4"], groupby='leiden_Res_0.3', stripplot=False, save="23_Exhaustion.png") # Exhaustion
sc.pl.violin(adata, keys=["KLRG1", "B3GAT1", "TIGIT", "CD160"], groupby='leiden_Res_0.3', stripplot=False, save="24_Senescense.png") # Senescence
sc.pl.violin(adata, keys=["SELL", "CXCR3", "CCR7", "CD44"], groupby='leiden_Res_0.3', stripplot=False, save="25_Other_Memory_markers.png") # Other memory markers
sc.pl.violin(adata, keys=["TBX21", "B3GAT1", "SELL", "CD44"], groupby='leiden_Res_0.3', stripplot=False, save="25.2_Other_Memory_markers.png") # Other memory markers 2
sc.pl.violin(adata, keys=["CD8A", "CD8B"], groupby='leiden_Res_0.3', stripplot=False, save="26_CD8_Complex.png") # CD8 complex marker genes
sc.pl.violin(adata, keys=["CD3D", "CD3G", "CD3E"], groupby='leiden_Res_0.3', stripplot=False, save="27_CD3_Complex.png") # CD3 complex marker genes

# %% Gen exploration -- Paula's choice
os.chdir("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/Plots/V4/scVI_markers_exploration_2")

sc.pl.stacked_violin(adata, var_names=["CD3D", "CD4", "CD8A"], groupby='leiden_Res_0.3', title="Marcadores basicos de celulas T", save="1_Basic_T_markers.pdf") # 1. Basic markers of T cells
sc.pl.stacked_violin(adata, var_names=["FOXP3", "IL2RA", "IL17A"], groupby='leiden_Res_0.3', title="Otros marcadores", save="2_Other_markers.pdf") # 2. Other markers
sc.pl.stacked_violin(adata, var_names=["NKG7", "GNLY", "IFNG", "PRF1", "GZMK", "GZMB"], groupby='leiden_Res_0.3', title="Citotoxicidad", save="3_Cytotoxicity.pdf") # 3. Cytotoxicity markers
sc.pl.stacked_violin(adata, var_names=["CD69", "CD27", "CD28", "HLA-DRA", "TNFRSF9"], groupby='leiden_Res_0.3', title="Activacion", save="4_Activation.pdf") # 4. Activation markers
sc.pl.stacked_violin(adata, var_names=["TCF7", "CCR7", "SELL"], groupby='leiden_Res_0.3', title="Memoria", save="5_Memory.pdf") # 5. Memory markers
sc.pl.stacked_violin(adata, var_names=["PTPRC"], groupby='leiden_Res_0.3', title="CD45", save="6_CD45.pdf") # 6. CD45

sc.pl.violin(adata, keys=["CD3D", "CD4", "CD8A"], groupby='leiden_Res_0.3', stripplot=False, save="1_Basic_T_markers.pdf") # 1. Basic markers of T cells
sc.pl.violin(adata, keys=["FOXP3", "IL2RA", "IL17A"], groupby='leiden_Res_0.3', stripplot=False, save="2_Other_markers.pdf") # 2. Other markers
sc.pl.violin(adata, keys=["NKG7", "GNLY", "IFNG", "PRF1", "GZMK", "GZMB"], groupby='leiden_Res_0.3', stripplot=False, save="3_Cytotoxicity.pdf") # 3. Cytotoxicity markers
sc.pl.violin(adata, keys=["CD69", "CD27", "CD28", "HLA-DRA", "TNFRSF9"], groupby='leiden_Res_0.3', stripplot=False, save="4_Activation.pdf") # 4. Activation markers
sc.pl.violin(adata, keys=["TCF7", "CCR7", "SELL"], groupby='leiden_Res_0.3', stripplot=False, save="5_Memory.pdf") # 5. Memory markers
sc.pl.violin(adata, keys=["PTPRC"], groupby='leiden_Res_0.3', stripplot=False, save="6_CD45.pdf") # 6. CD45

sc.pl.umap(adata, color=["CD3D", "CD4", "CD8A"], frameon=False, ncols=3, vmin=0, vmax="p99", save="1_Basic_T_markers.png") # 1. Basic markers of T cells
sc.pl.umap(adata, color=["FOXP3", "IL2RA", "IL17A"], frameon=False, ncols=2, vmin=0, vmax="p99", save="2_Other_markers.png") # 2. Other markers
sc.pl.umap(adata, color=["NKG7", "GNLY", "IFNG", "PRF1", "GZMK", "GZMB"], frameon=False, ncols=2, vmin=0, vmax="p99", save="3_Cytotoxicity.png") # 3. Cytotoxicity markers
sc.pl.umap(adata, color=["CD69", "CD27", "CD28", "HLA-DRA", "TNFRSF9"], frameon=False, ncols=2, vmin=0, vmax="p99", save="4_Activation.png") # 4. Activation markers
sc.pl.umap(adata, color=["TCF7", "CCR7", "SELL"], frameon=False, ncols=2, vmin=0, vmax="p99", save="5_Memory.png") # 5. Memory markers
sc.pl.umap(adata, color=["PTPRC"], frameon=False, ncols=2, vmin=0, vmax="p99", save="6_CD45.png") # 6. CD45

# %% Extra QC --- Ribosomal and hemoglobin scoring
os.chdir("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/Plots/V4/scVI_markers_exploration_3")

# Ribosomal genes
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))

# Hemoglobin genes
adata.var["hb"] = adata.var_names.str.contains(("^HB[^(P)]"))

# Mitocondrial genes
adata.var["mt"] = adata.var_names.str.startswith(("MT-"))

# Apoptosis genes (gsea: KEGG_APOPTOSIS + REACTOME_APOPTOSIS)
Apoptosis_genes = ['CASP3', 'CASP8', 'BAD', 'FAS', 'PPP3CC', 'FASLG', 'CASP9', 'AKT1', 'TRAF2', 'BAK1', 'CSE1L', 'DFFB', 'PPP3R1', 'TP53', 'BIRC2']

sc.pl.stacked_violin(adata, var_names=Apoptosis_genes, groupby='leiden_Res_0.3', title="Apoptosis markers", save="Apoptosis_markers.pdf") # Apoptosis markers
sc.pl.umap(adata, color=Apoptosis_genes, frameon=False, ncols=3, vmin=0, vmax="p99", save="Apoptosis_markers.png") # Apoptosis markers

adata.var["apop"] = adata.var_names.isin(Apoptosis_genes)

sc.pp.calculate_qc_metrics(
    adata, qc_vars=["ribo", "hb", "mt", "apop"], inplace=True, percent_top=[20]
)
adata

sc.pl.umap(adata, color=["pct_counts_ribo", "pct_counts_hb", "pct_counts_mt"], frameon=False, vmin=0, vmax="p99", save="_QC_extra_params.pdf")
sc.pl.umap(adata, color=["pct_counts_ribo", "pct_counts_hb", "pct_counts_mt"], frameon=False, vmin=0, vmax="p99", save="_QC_extra_params.png")

# Necrosis genes (gsea: REACTOME_REGULATED_NECROSIS)
Necrosis_genes = ["GZMB", "IL18", "UBC", "UBA52", "IL1A", "BAX", "CFLAR", "BIRC3", "TP63", "BAK1"]

sc.pl.stacked_violin(adata, var_names=Necrosis_genes, groupby='leiden_Res_0.3', title="Necrosis markers", save="Necrosis_markers.pdf") # Necrosis markers
sc.pl.umap(adata, color=Necrosis_genes, frameon=False, ncols=2, vmin=0, vmax="p99", save="Necrosis_markers.png") # Necrosis markers

# %% Get deeper in the possible erytroid contamination.
os.chdir("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/Plots/V4/scVI_markers_exploration_3")
sc.pl.umap(adata, color=("HBB", "HBA1", "HBA2", "HBD", "HBM", "AHSP", "ALAS2", "CA1", "SLC4A1", "IFIT1B", "TRIM58", "SELENBP1", "TMCC2"), frameon=False, vmin=0, vmax="p99", save="Erythroid_contamination.png")

sc.tl.score_genes(adata, gene_list=["HBB", "HBA1", "HBA2", "HBD", "HBM", "AHSP", "ALAS2", "CA1", "SLC4A1", "IFIT1B", "TRIM58", "SELENBP1", "TMCC2"], score_name="score_ery")
sc.pl.umap(adata, color="score_ery", frameon=False, ncols=2, vmin=0, vmax="p99", save="Erythroid_contamination_2.png")
# There are erythocytes that should be removed. But I don't see contamination in the rest, they are only expressing HBB in a much lower manner.

sc.tl.score_genes(adata, gene_list=['CASP3', 'CASP8', 'BAD', 'FAS', 'PPP3CC', 'AKT1', 'TRAF2', 'BAK1', 'CSE1L', 'PPP3R1', 'TP53', 'BIRC2'], score_name="score_apop")
sc.pl.umap(adata, color="score_apop", frameon=False, ncols=2, vmin=0, vmax="p99", save="Apoptosis_score.png")
sc.pl.umap(adata, color="score_apop", frameon=False, ncols=2, vmin=0, vmax="p99", save="Apoptosis_score.pdf")
sc.pl.umap(adata, color="score_apop", frameon=False, ncols=2, save="Apoptosis_score_2.pdf")

# %% Extra QC --- Extra monocytes and IACs markers
os.chdir("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/Plots/V4/scVI_markers_exploration_3")

sc.pl.stacked_violin(adata, var_names=["C1QA", "PSAP", "TGFBI", "BCL2A1", "JAKMIP2", "OTOA", "PLPP3", "CXCL8"], groupby='leiden_Res_0.3', title="Extra monocytes markers", save="Extra_monocytes_markers.pdf") # Extra monocytes markers
sc.pl.umap(adata, color=["C1QA", "PSAP", "TGFBI", "BCL2A1", "JAKMIP2", "OTOA", "PLPP3", "CXCL8"], frameon=False, ncols=2, vmin=0, vmax="p99", save="Extra_monocytes_markers.png") # Extra monocytes markers

sc.pl.stacked_violin(adata, var_names=["CD68", "LYZ", "SPI1", "LILRB4", "SIRPA"], groupby='leiden_Res_0.3', save="1_IACs.pdf") # First IACs graph (look at Deng's paper)
sc.pl.umap(adata, color=["CD68", "LYZ", "SPI1", "LILRB4", "SIRPA"], frameon=False, ncols=2, vmin=0, vmax="p99", save="1_IACs.png") # First IACs graph (look at Deng's paper)

sc.pl.stacked_violin(adata, var_names=["IL1B", "CXCL8"], groupby='leiden_Res_0.3', save="2_IACs.pdf") # Second IACs graph (look at Deng's paper)
sc.pl.umap(adata, color=["IL1B", "CXCL8"], frameon=False, ncols=2, vmin=0, vmax="p99", save="2_IACs.png") # Second IACs graph (look at Deng's paper)

# %% Remove erytrocytes detected in Atlas
Hemo_list1 = adata[adata[: , 'HBB'].X > 0.5, :].obs_names
Hemo_list2 = adata[adata[: , 'HBA1'].X > 0.5, :].obs_names
Hemo_list3 = adata[adata[: , 'HBA2'].X > 0.5, :].obs_names

Hemo_to_remove = list(set(Hemo_list1) & set(Hemo_list2) & set(Hemo_list3))

adata_filtered = adata[~adata.obs_names.isin(Hemo_to_remove)]

adata = adata_filtered.copy()

# %% Graphs to detect if contamination has been removed
os.chdir("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/Plots/V4/scVI_markers_exploration_3")
sc.pl.umap(adata, color=("HBB", "HBA1", "HBA2", "HBD", "HBM", "AHSP", "ALAS2", "CA1", "SLC4A1", "IFIT1B", "TRIM58", "SELENBP1", "TMCC2"), frameon=False, vmin=0, vmax="p99", save="Erythroid_contamination_bis.png")

sc.tl.score_genes(adata, gene_list=["HBB", "HBA1", "HBA2", "HBD", "HBM", "AHSP", "ALAS2", "CA1", "SLC4A1", "IFIT1B", "TRIM58", "SELENBP1", "TMCC2"], score_name="score_ery")

sc.pl.umap(adata, color="score_ery", frameon=False, ncols=2, vmin=0, vmax="p99", save="Erythroid_contamination_2_bis.png")

# %% Atlas contribution (in nº of cells)

## Set display options to show all rows and columns
# pd.set_option('display.max_rows', None)
# pd.set_option('display.max_columns', None)

os.chdir("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/scVI/V4/Data")
adata.obs[["orig_ident"]].value_counts(sort=False).to_csv("Orig_ident_final_values.csv")
adata.obs[["Product_norm"]].value_counts(sort=False).to_csv("Product_norm_final_values.csv")

## Reset display options to default if needed
# pd.reset_option('display.max_rows')
# pd.reset_option('display.max_columns')

# %% UMAP plots save
os.chdir(
    "/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/Plots/V4/scVI_integ"
)
Features_to_explore = [
    "celltypist_majority_voting",
    "leiden_Res_0.3",
    "orig_ident",
]

for i in range(len(Features_to_explore)):
    feature_name = Features_to_explore[i]
    filename = f"_{feature_name}_scVI_V4_WO_ery.pdf"
    filename2 = f"_{feature_name}_scVI_V4_WO_ery.png"

    sc.pl.umap(adata, color=feature_name, frameon=False, save=filename)

    sc.pl.umap(adata, color=feature_name, frameon=False, save=filename2)

# %% Load/Save state 2.5
########################################################################
#
# State 2.5 Load/Save
# 
######################################################################## 
if Save_h5ad_2_5:
    os.chdir(
        "/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/scVI/V4"
    )
    adata.write("Python_scVI_adata_V4_state2_5.h5ad")
    print("Correctly saved")
else:
    os.chdir("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/scVI/V4")
    adata = sc.read_h5ad("Python_scVI_adata_V4_state2_5.h5ad")
    print("Correctly loaded")

# %% Code to visualize where a specific cluster is located
### Search for clusters at 0.4 resolution
adata.obs['cluster_dummy9'] = adata.obs['leiden_Res_0.4'] == adata.obs['leiden_Res_0.4'].cat.categories[9]
adata.obs['cluster_dummy9'] = adata.obs['cluster_dummy9']*1
sc.pl.umap(adata, color='cluster_dummy9', frameon=False)

adata.obs['cluster_dummy10'] = adata.obs['leiden_Res_0.4'] == adata.obs['leiden_Res_0.4'].cat.categories[10]
adata.obs['cluster_dummy10'] = adata.obs['cluster_dummy10']*1
sc.pl.umap(adata, color='cluster_dummy10', frameon=False)

### Search for clusters at 0.3 resolution
adata.obs['cluster_dummy6'] = adata.obs['leiden_Res_0.3'] == adata.obs['leiden_Res_0.3'].cat.categories[6]
adata.obs['cluster_dummy6'] = adata.obs['cluster_dummy6']*1
sc.pl.umap(adata, color='cluster_dummy6', frameon=False)

# %% Load state 2.8
########################################################################
#
# State 2.8 Load
# 
######################################################################## 
os.chdir("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/scVI/V4")
adata = sc.read_h5ad("Python_scVI_adata_V4_state2_8.h5ad")

# %% Annotation 1 --- Version 4
########################################################################
#
# Annotation - Resolution "low" - low number of clusters, more general
# 
######################################################################## 

os.chdir("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/Plots/V4/Annotation/Low_Res")

cl_annotation_low = {
    "0": "CD4 memory", # Subclusterize (Missing for example T regs)
    "1": "CD8 cytotoxic", # Subclusterize a bit to check subpopulations
    "2": "Proliferative T cells",
    "3": "CD4+ CD8+ high-GATA3", # Subclusterize mix between CD4 and CD8 cells
    "4": "Ribosomal enriched",
    "5": "Monocyte-like T cells",
    "6": "Apoptotic T cells"
}

adata.obs["manual_celltype_annotation_low"] = adata.obs["leiden_Res_0.3"].map(cl_annotation_low)

sc.pl.umap(adata, color = ["manual_celltype_annotation_low"], frameon=False, save="manual_celltype_annotation_low.png")
sc.pl.umap(adata, color = ["manual_celltype_annotation_low"], frameon=False, save="manual_celltype_annotation_low.pdf")

# %% Load/Save state 3
########################################################################
#
# State 3 Load/Save
# 
######################################################################## 
if Save_h5ad_3:
    os.chdir(
        "/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/scVI/V4"
    )
    adata.write("Python_scVI_adata_V4_state3.h5ad")
    print("Correctly saved")
else:
    os.chdir("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/scVI/V4")
    adata = sc.read_h5ad("Python_scVI_adata_V4_state3.h5ad")
    print("Correctly loaded")

# %% Atlas colorized based on metadata -- Generation of datasets
os.chdir("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/scVI/V4")
sc_Metadata_to_python_v33 = pd.read_csv("sc_Metadata_to_python_v33_V4.csv", sep=";") # Check it is the latest version
sc_Metadata_to_python_v33.rename(columns={'Norm_Sample_Name': 'Product_norm'}, inplace=True)
metadata_bis = adata.obs.copy()

metadata_bis['row_number'] = range(1, len(metadata_bis) + 1)
print(metadata_bis['row_number'].is_monotonic_increasing)

merged_df = pd.merge(metadata_bis, sc_Metadata_to_python_v33, on="Product_norm", how="right")
print(merged_df['row_number'].is_monotonic_increasing)

merged_df.index = adata.obs.index

adata2 = adata.copy()
adata2.obs = merged_df.copy()

adata3 = adata.copy()
adata3.obs = merged_df.copy()
adata3 = adata3[adata3.obs["CAR_Gen"] == "Second"]

adata4 = adata2.copy()
adata4.obs = merged_df.copy()
dummy_column = merged_df["Time_Point_Ranges"] == "Patient_Infusion_Product"
adata4.obs['Dummy_Column'] = dummy_column.astype(int)
adata4.obs.loc[adata4.obs['Dummy_Column'] == 0, 'Therapy_Toxicities'] = pd.NA
adata4.obs.loc[adata4.obs['Dummy_Column'] == 0, 'CRS'] = pd.NA
adata4.obs.loc[adata4.obs['Dummy_Column'] == 0, 'CRS_Grade_Range'] = pd.NA
adata4.obs.loc[adata4.obs['Dummy_Column'] == 0, 'ICANS'] = pd.NA
adata4.obs.loc[adata4.obs['Dummy_Column'] == 0, 'ICANS_Grade_Range'] = pd.NA

# %% Atlas colorized based on metadata -- Plotting
os.chdir("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/Plots/V4/Colorized_Atlas_Metadata")

sc.pl.umap(adata2, color = ["Time_Point_Ranges"], frameon=False, save="1_Time_Point_Ranges.png")
sc.pl.umap(adata2, color = ["Time_Point_Ranges"], frameon=False, save="1_Time_Point_Ranges.pdf")

sc.pl.umap(adata2, color = ["Sex"], frameon=False, save="2_Sex.png")
sc.pl.umap(adata2, color = ["Sex"], frameon=False, save="2_Sex.pdf")

sc.pl.umap(adata2, color = ["Age_Range"], frameon=False, save="3_Age_Range.png")
sc.pl.umap(adata2, color = ["Age_Range"], frameon=False, save="3_Age_Range.pdf")

sc.pl.umap(adata2, color = ["ScFv"], frameon=False, save="5_ScFv.png")
sc.pl.umap(adata2, color = ["ScFv"], frameon=False, save="5_ScFv.pdf")

sc.pl.umap(adata2, color = ["CAR_Construct"], frameon=False, save="6_CAR_Construct.png")
sc.pl.umap(adata2, color = ["CAR_Construct"], frameon=False, save="6_CAR_Construct.pdf")

sc.pl.umap(adata2, color = ["CAR_Gen"], frameon=False, save="7_CAR_Gen.png")
sc.pl.umap(adata2, color = ["CAR_Gen"], frameon=False, save="7_CAR_Gen.pdf")

sc.pl.umap(adata3, color = ["Costim_Domain_1"], frameon=False, save="8_Costim_Domain_1.png")
sc.pl.umap(adata3, color = ["Costim_Domain_1"], frameon=False, save="8_Costim_Domain_1.pdf")

sc.pl.umap(adata2, color = ["Person_source_General"], frameon=False, save="9.1_Person_source_General.png")
sc.pl.umap(adata2, color = ["Person_source_General"], frameon=False, save="9.1_Person_source_General.pdf")

sc.pl.umap(adata2, color = ["Person_source_Specific"], frameon=False, save="9.2_Person_source_Specific.png")
sc.pl.umap(adata2, color = ["Person_source_Specific"], frameon=False, save="9.2_Person_source_Specific.pdf")

sc.pl.umap(adata4, color = ["Therapy_Toxicities"], frameon=False, save="10_Therapy_Toxicities_only_IP.png")
sc.pl.umap(adata4, color = ["Therapy_Toxicities"], frameon=False, save="10_Therapy_Toxicities_only_IP.pdf")

sc.pl.umap(adata2, color = ["Therapy_Toxicities"], frameon=False, save="10_Therapy_Toxicities.png")
sc.pl.umap(adata2, color = ["Therapy_Toxicities"], frameon=False, save="10_Therapy_Toxicities.pdf")

sc.pl.umap(adata2, color = ["Stimulated"], frameon=False, save="11_Stimulated.png")
sc.pl.umap(adata2, color = ["Stimulated"], frameon=False, save="11_Stimulated.pdf")

sc.pl.umap(adata2, color = ["Starting_Point"], frameon=False, save="12_Starting_Point.png")
sc.pl.umap(adata2, color = ["Starting_Point"], frameon=False, save="12_Starting_Point.pdf")

sc.pl.umap(adata2, color = ["Technology"], frameon=False, save="13_Technology.png")
sc.pl.umap(adata2, color = ["Technology"], frameon=False, save="13_Technology.pdf")

sc.pl.umap(adata2, color = ["Antigen"], frameon=False, save="14_Antigen.png")
sc.pl.umap(adata2, color = ["Antigen"], frameon=False, save="14_Antigen.pdf")

sc.pl.umap(adata2, color = ["Source"], frameon=False, save="15_Source.png")
sc.pl.umap(adata2, color = ["Source"], frameon=False, save="15_Source.pdf")

sc.pl.umap(adata2, color = ["STATUS"], frameon=False, save="16_STATUS.png")
sc.pl.umap(adata2, color = ["STATUS"], frameon=False, save="16_STATUS.pdf")

sc.pl.umap(adata4, color = ["CRS"], frameon=False, save="17_CRS_only_IP.png")
sc.pl.umap(adata4, color = ["CRS"], frameon=False, save="17_CRS_only_IP.pdf")

sc.pl.umap(adata2, color = ["CRS"], frameon=False, save="17_CRS.png")
sc.pl.umap(adata2, color = ["CRS"], frameon=False, save="17_CRS.pdf")

sc.pl.umap(adata4, color = ["CRS_Grade_Range"], frameon=False, save="18_CRS_Grade_Range_only_IP.png")
sc.pl.umap(adata4, color = ["CRS_Grade_Range"], frameon=False, save="18_CRS_Grade_Range_only_IP.pdf")

sc.pl.umap(adata2, color = ["CRS_Grade_Range"], frameon=False, save="18_CRS_Grade_Range.png")
sc.pl.umap(adata2, color = ["CRS_Grade_Range"], frameon=False, save="18_CRS_Grade_Range.pdf")

sc.pl.umap(adata4, color = ["ICANS"], frameon=False, save="19_ICANS_only_IP.png")
sc.pl.umap(adata4, color = ["ICANS"], frameon=False, save="19_ICANS_only_IP.pdf")

sc.pl.umap(adata2, color = ["ICANS"], frameon=False, save="19_ICANS.png")
sc.pl.umap(adata2, color = ["ICANS"], frameon=False, save="19_ICANS.pdf")

sc.pl.umap(adata4, color = ["ICANS_Grade_Range"], frameon=False, save="20_ICANS_Grade_Range_only_IP.png")
sc.pl.umap(adata4, color = ["ICANS_Grade_Range"], frameon=False, save="20_ICANS_Grade_Range_only_IP.pdf")

sc.pl.umap(adata2, color = ["ICANS_Grade_Range"], frameon=False, save="20_ICANS_Grade_Range.png")
sc.pl.umap(adata2, color = ["ICANS_Grade_Range"], frameon=False, save="20_ICANS_Grade_Range.pdf")

sc.pl.umap(adata2, color = ["CAR_Detection_Method"], frameon=False, save="21_CAR_Detection_Method.png")
sc.pl.umap(adata2, color = ["CAR_Detection_Method"], frameon=False, save="21_CAR_Detection_Method.pdf")

# %% Subclustering -- To get more resolution
os.chdir("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/Subclustering_V4")

adata_CD4 = adata2[adata2.obs["manual_celltype_annotation_low"] == "CD4 memory"].copy()
adata_CD8 = adata2[adata2.obs["manual_celltype_annotation_low"] == "CD8 cytotoxic"].copy()
adata_GATA3 = adata2[adata2.obs["manual_celltype_annotation_low"] == "CD4+ CD8+ high-GATA3"].copy()

# Delete previous leiden information - create a list of columns to delete
columns_to_delete_CD4 = [col for col in adata_CD4.obs.columns if 'leiden' in col]
columns_to_delete_CD8 = [col for col in adata_CD8.obs.columns if 'leiden' in col]
columns_to_delete_GATA3 = [col for col in adata_GATA3.obs.columns if 'leiden' in col]

# Delete the selected columns from adata.obs
adata_CD4.obs = adata_CD4.obs.drop(columns=columns_to_delete_CD4)
adata_CD8.obs = adata_CD8.obs.drop(columns=columns_to_delete_CD8)
adata_GATA3.obs = adata_GATA3.obs.drop(columns=columns_to_delete_GATA3)

# Save objects
adata_CD4.write("adata_CD4.h5ad")
adata_CD8.write("adata_CD8.h5ad")
adata_GATA3.write("adata_GATA3.h5ad")

# %% End of script
