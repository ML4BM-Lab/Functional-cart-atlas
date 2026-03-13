###############################################################################
###############################################################################

# Program: 2.4_Subclustering_V4_GATA3.py
# Author: Sergio Cámara Peña
# Date: 16/02/2024
# Version: FINAL

###############################################################################
###############################################################################

# %% Load all the needed libraries
import scanpy as sc
import os
import random
from scanpy_cluster_proportions import get_cluster_proportions, plot_cluster_proportions
import seaborn as sns
import pandas as pd
import gseapy

# %% Switches
Save_h5ad = False

# %% Set a random seed
random.seed(2504)

# %% Load data
if Save_h5ad:
    os.chdir("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/Subclustering_V4")
    adata_GATA3 = sc.read_h5ad("adata_GATA3.h5ad")

# %% Calculate neighbours and leiden resolutions
if Save_h5ad:
    sc.pp.neighbors(adata_GATA3, use_rep="X_scVI")
    sc.tl.leiden(adata_GATA3, key_added="leiden_Sub_Res_0.2", resolution=0.2)
    sc.tl.leiden(adata_GATA3, key_added="leiden_Sub_Res_0.3", resolution=0.3)
    sc.tl.leiden(adata_GATA3, key_added="leiden_Sub_Res_0.4", resolution=0.4)
    sc.tl.leiden(adata_GATA3, key_added="leiden_Sub_Res_0.5", resolution=0.5)
    sc.tl.leiden(adata_GATA3, key_added="leiden_Sub_Res_0.6", resolution=0.6)

# %% Save/Load adata object to state 1 h5ad
########################################################################
#
# State 1 Load/Save
# 
######################################################################## 
if Save_h5ad:
    os.chdir("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/Subclustering_V4")
    adata_GATA3.write("adata_GATA3_state1.h5ad")
    print("Correctly saved")
else:
    os.chdir("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/Subclustering_V4")
    adata_GATA3 = sc.read_h5ad("adata_GATA3_state1.h5ad")
    print("Correctly loaded")

# %% Check which resolution to choose for subclustering using Clustree
if False:
    import invoke
    os.chdir("/home/scamara/data/scamara/Atlas_Mieloma_Multiple/Codigo/Codigo_datasets_atlas/Datasets_Integration/Tests_labs")
    CLUSTREE_SCRIPT = "Leiden_resol_lab_clustree_function.R"
    ADATA_PATH = "/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/Subclustering_V4/adata_GATA3_state1.h5ad"
    CLUSTREE_OUT = "/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/Subclustering_V4/Clustree_res/clustree_GATA3.png"
    Carrera = invoke.run(f"Rscript {CLUSTREE_SCRIPT} {ADATA_PATH} {CLUSTREE_OUT}")

# After examining the resulting graph, I decided to use 0.3 resolution

# %% Clustering Plots --- GATA3 subclustering UMAP
os.chdir("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/Plots/V4/Subclustering/GATA3")
sc.pl.umap(adata_GATA3, color="leiden_Sub_Res_0.3", frameon=False, save = "0_GATA3_subclustering.png")
sc.pl.umap(adata_GATA3, color="leiden_Sub_Res_0.3", frameon=False, save = "0_GATA3_subclustering.pdf")

# %% Clustering Plots --- Dataset contribution to cluster & cluster composition by dataset
cluster_props = get_cluster_proportions(adata_GATA3, cluster_key="leiden_Sub_Res_0.3", sample_key="orig_ident")
f1 = plot_cluster_proportions(cluster_props, xlabel_rotation=90, cluster_palette=sns.color_palette("hls", adata_GATA3.obs["leiden_Sub_Res_0.3"].nunique()))
f1.savefig("QC_Datasets_to_Cluster_distrib.pdf")
f1.savefig("QC_Datasets_to_Cluster_distrib.png")

cluster_props_2 = get_cluster_proportions(adata_GATA3, cluster_key="orig_ident", sample_key="leiden_Sub_Res_0.3")
f2 = plot_cluster_proportions(cluster_props_2, xlabel_rotation=90, cluster_palette=sns.color_palette("hls", adata_GATA3.obs["orig_ident"].nunique()))
f2.savefig("QC_Clusters_to_Dataset_distrib.pdf")
f2.savefig("QC_Clusters_to_Dataset_distrib.png")

# %% Finding marker genes by subcluster
sc.pp.log1p(adata_GATA3) # Calling this function before ranking genes avoids the need to use the raw data and fixes the dictionary issue
sc.tl.rank_genes_groups(adata_GATA3, "leiden_Sub_Res_0.3", method="wilcoxon") # This function expects the data to be log-transformed

os.chdir("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/Plots/V4/Subclustering/GATA3")
sc.pl.rank_genes_groups(adata_GATA3, n_genes=25, sharey=False, save=".pdf")
sc.pl.rank_genes_groups(adata_GATA3, n_genes=25, sharey=False, save=".png")

sc.pl.rank_genes_groups_stacked_violin(adata_GATA3, n_genes=4, cmap='RdBu_r', groupby="leiden_Sub_Res_0.3", save="scvi.pdf")

# %% Marker genes by cluster --- Table
os.chdir("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/Subclustering_V4/Data/GATA3")
result = adata_GATA3.uns['rank_genes_groups']
groups = result['names'].dtype.names

Marker_genes_x_clus_df = pd.DataFrame(
    {group + '_' + key[:1]: result[key][group]
    for group in groups for key in ['names', 'pvals']})

Marker_genes_x_clus_df
Marker_genes_x_clus_df.to_csv("Marker_genes_per_cluster.csv")

# %% Gen exploration -- Paula's choice
os.chdir("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/Plots/V4/Subclustering/GATA3")

sc.pl.stacked_violin(adata_GATA3, var_names=["CD3D", "CD4", "CD8A"], groupby='leiden_Sub_Res_0.3', title="Marcadores basicos de celulas T", save="1_Basic_T_markers.pdf") # 1. Basic T cell markers
sc.pl.stacked_violin(adata_GATA3, var_names=["FOXP3", "IL2RA", "IL17A"], groupby='leiden_Sub_Res_0.3', title="Otros marcadores", save="2_Other_markers.pdf") # 2. Other markers
sc.pl.stacked_violin(adata_GATA3, var_names=["NKG7", "GNLY", "IFNG", "PRF1", "GZMK", "GZMB"], groupby='leiden_Sub_Res_0.3', title="Citotoxicidad", save="3_Cytotoxicity.pdf") # 3. Cytotoxicity markers
sc.pl.stacked_violin(adata_GATA3, var_names=["CD69", "CD27", "CD28", "HLA-DRA", "TNFRSF9"], groupby='leiden_Sub_Res_0.3', title="Activacion", save="4_Activation.pdf") # 4. Activation markers
sc.pl.stacked_violin(adata_GATA3, var_names=["TCF7", "CCR7", "SELL"], groupby='leiden_Sub_Res_0.3', title="Memoria", save="5_Memory.pdf") # 5. Memory markers
sc.pl.stacked_violin(adata_GATA3, var_names=["PTPRC"], groupby='leiden_Sub_Res_0.3', title="CD45", save="6_CD45.pdf") # 6. CD45

sc.pl.violin(adata_GATA3, keys=["CD3D", "CD4", "CD8A"], groupby='leiden_Sub_Res_0.3', stripplot=False, save="1_Basic_T_markers.pdf") # 1. Basic T cell markers
sc.pl.violin(adata_GATA3, keys=["FOXP3", "IL2RA", "IL17A"], groupby='leiden_Sub_Res_0.3', stripplot=False, save="2_Other_markers.pdf") # 2. Other markers
sc.pl.violin(adata_GATA3, keys=["NKG7", "GNLY", "IFNG", "PRF1", "GZMK", "GZMB"], groupby='leiden_Sub_Res_0.3', stripplot=False, save="3_Cytotoxicity.pdf") # 3. Cytotoxicity markers
sc.pl.violin(adata_GATA3, keys=["CD69", "CD27", "CD28", "HLA-DRA", "TNFRSF9"], groupby='leiden_Sub_Res_0.3', stripplot=False, save="4_Activation.pdf") # 4. Activation markers
sc.pl.violin(adata_GATA3, keys=["TCF7", "CCR7", "SELL"], groupby='leiden_Sub_Res_0.3', stripplot=False, save="5_Memory.pdf") # 5. Memory markers
sc.pl.violin(adata_GATA3, keys=["PTPRC"], groupby='leiden_Sub_Res_0.3', stripplot=False, save="6_CD45.pdf") # 6. CD45

sc.pl.umap(adata_GATA3, color=["CD3D", "CD4", "CD8A"], frameon=False, ncols=2, vmin=0, vmax="p99", save="1_Basic_T_markers.png") # 1. Basic T cell markers
sc.pl.umap(adata_GATA3, color=["FOXP3", "IL2RA", "IL17A"], frameon=False, ncols=2, vmin=0, vmax="p99", save="2_Other_markers.png") # 2. Other markers
sc.pl.umap(adata_GATA3, color=["NKG7", "GNLY", "IFNG", "PRF1", "GZMK", "GZMB"], frameon=False, ncols=2, vmin=0, vmax="p99", save="3_Cytotoxicity.png") # 3. Cytotoxicity markers
sc.pl.umap(adata_GATA3, color=["CD69", "CD27", "CD28", "HLA-DRA", "TNFRSF9"], frameon=False, ncols=2, vmin=0, vmax="p99", save="4_Activation.png") # 4. Activation markers
sc.pl.umap(adata_GATA3, color=["TCF7", "CCR7", "SELL"], frameon=False, ncols=2, vmin=0, vmax="p99", save="5_Memory.png") # 5. Memory markers
sc.pl.umap(adata_GATA3, color=["PTPRC"], frameon=False, ncols=2, vmin=0, vmax="p99", save="6_CD45.png") # 6. CD45

# Extra 
sc.pl.umap(adata_GATA3, color=["NCAM1"], frameon=False, ncols=2, vmin=0, vmax="p99", save="7_NCAM1.png") # 7. NCAM1
sc.pl.stacked_violin(adata_GATA3, var_names=["GATA3", "IL4"], groupby='leiden_Sub_Res_0.3', title="Th2", save="8_Th2.pdf") # 8. Th2

# %% Gene Set Analysis -- Preprocessing steps
gene_set_names = gseapy.get_library_name(organism='Human')
print(gene_set_names)

os.chdir("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/Subclustering_V4/Data/GATA3")
Marker_genes_x_clus_df = pd.read_csv("Marker_genes_per_cluster.csv", index_col=0)

# %% Gene Set Analysis -- Analysis for clusters
os.chdir("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/Plots/V4/Subclustering/GATA3/GSA")

# Define a list of clusters and their corresponding gene lists
clusters = [
    {"name": "Cluster_0", "genes": Marker_genes_x_clus_df.iloc[:100, 0]},
    {"name": "Cluster_1", "genes": Marker_genes_x_clus_df.iloc[:100, 2]},
    {"name": "Cluster_2", "genes": Marker_genes_x_clus_df.iloc[:100, 4]},
    {"name": "Cluster_3", "genes": Marker_genes_x_clus_df.iloc[:100, 6]},
    {"name": "Cluster_4", "genes": Marker_genes_x_clus_df.iloc[:100, 8]},
]

# Define the gene set databases
gene_set_databases = {
    "GO_BP": "GO_Biological_Process_2023",
    "GO_MF": "GO_Molecular_Function_2023",
    "GO_CC": "GO_Cellular_Component_2023",
    "KEGG": "KEGG_2021_Human"
}

# Set the output directory
output_directory = "/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/Plots/V4/Subclustering/GATA3/GSA"

# Loop through clusters and gene set databases
for cluster in clusters:
    cluster_name = cluster["name"]
    cluster_genes = cluster["genes"]

    for database_name, gene_set in gene_set_databases.items():
        os.chdir(output_directory)

        # Perform gene set analysis
        enr_res = gseapy.enrichr(
            gene_list=cluster_genes,
            organism='Human',
            gene_sets=gene_set,
            cutoff=0.5
        )

        # Generate a bar plot for the results
        plot_title = f"{cluster_name}_{database_name}_2023"
        plot_filename = f"enr_res_{cluster_name}_{database_name}.png"
        gseapy.barplot(enr_res.res2d, title=plot_title, ofname=plot_filename)

# %% Subclustering annotation
cl_annotation_high = {
    "0": "CD8 cytotoxic",
    "1": "CD4 central memory",
    "2": "CD8 cytotoxic",
    "3": "CD8 cytotoxic",
    "4": "CD8 effector memory",
}

adata_GATA3.obs["manual_celltype_annotation_high"] = adata_GATA3.obs["leiden_Sub_Res_0.3"].map(cl_annotation_high)

os.chdir("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/Plots/V4/Subclustering/GATA3")
sc.pl.umap(adata_GATA3, color = ["manual_celltype_annotation_high"], frameon=False, save="manual_celltype_annotation_high_GATA3.png")

# %% Save the subclustering information
os.chdir("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/Subclustering_V4/Data/GATA3")
adata_GATA3.obs["manual_celltype_annotation_high"].to_csv("adata_GATA3_manual_celltype_annotation_high.csv")

# %% End of script