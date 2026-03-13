###############################################################################
###############################################################################

# Program: 2.4_Subclustering_V4_CD4.py
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
    adata_CD4 = sc.read_h5ad("adata_CD4.h5ad")

# %% Calculate neighbours and leiden resolutions
if Save_h5ad:
    sc.pp.neighbors(adata_CD4, use_rep="X_scVI")
    sc.tl.leiden(adata_CD4, key_added="leiden_Sub_Res_0.3", resolution=0.3)
    sc.tl.leiden(adata_CD4, key_added="leiden_Sub_Res_0.4", resolution=0.4)
    sc.tl.leiden(adata_CD4, key_added="leiden_Sub_Res_0.5", resolution=0.5)
    sc.tl.leiden(adata_CD4, key_added="leiden_Sub_Res_0.6", resolution=0.6)
    sc.tl.leiden(adata_CD4, key_added="leiden_Sub_Res_0.7", resolution=0.7)
    sc.tl.leiden(adata_CD4, key_added="leiden_Sub_Res_0.8", resolution=0.8)
    sc.tl.leiden(adata_CD4, key_added="leiden_Sub_Res_0.9", resolution=0.9)

# %% Save/Load adata object to state 1 h5ad
########################################################################
#
# State 1 Load/Save
# 
######################################################################## 
if Save_h5ad:
    os.chdir("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/Subclustering_V4")
    adata_CD4.write("adata_CD4_state1.h5ad")
    print("Correctly saved")
else:
    os.chdir("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/Subclustering_V4")
    adata_CD4 = sc.read_h5ad("adata_CD4_state1.h5ad")
    print("Correctly loaded")

# %% Check which resolution to choose for subclustering using Clustree
if False:
    import invoke
    os.chdir("/home/scamara/data/scamara/Atlas_Mieloma_Multiple/Codigo/Codigo_datasets_atlas/Datasets_Integration/Tests_labs")
    CLUSTREE_SCRIPT = "Leiden_resol_lab_clustree_function.R"
    ADATA_PATH = "/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/Subclustering_V4/adata_CD4_state1.h5ad"
    CLUSTREE_OUT = "/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/Subclustering_V4/Clustree_res/clustree_CD4.png"
    Carrera = invoke.run(f"Rscript {CLUSTREE_SCRIPT} {ADATA_PATH} {CLUSTREE_OUT}")

# After examining the resulting graph, I decided to use the 0.7 resolution

# %% Clustering Plots --- CD4 subclustering UMAP
os.chdir("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/Plots/V4/Subclustering/CD4")
sc.pl.umap(adata_CD4, color="leiden_Sub_Res_0.7", frameon=False, save = "0_CD4_subclustering.png")
sc.pl.umap(adata_CD4, color="leiden_Sub_Res_0.7", frameon=False, save = "0_CD4_subclustering.pdf")

# %% Clustering Plots --- Dataset contribution to cluster & cluster composition by dataset
cluster_props = get_cluster_proportions(adata_CD4, cluster_key="leiden_Sub_Res_0.7", sample_key="orig_ident")
f1 = plot_cluster_proportions(cluster_props, xlabel_rotation=90, cluster_palette=sns.color_palette("hls", adata_CD4.obs["leiden_Sub_Res_0.7"].nunique()))
f1.savefig("QC_Datasets_to_Cluster_distrib.pdf")
f1.savefig("QC_Datasets_to_Cluster_distrib.png")

cluster_props_2 = get_cluster_proportions(adata_CD4, cluster_key="orig_ident", sample_key="leiden_Sub_Res_0.7")
f2 = plot_cluster_proportions(cluster_props_2, xlabel_rotation=90, cluster_palette=sns.color_palette("hls", adata_CD4.obs["orig_ident"].nunique()))
f2.savefig("QC_Clusters_to_Dataset_distrib.pdf")
f2.savefig("QC_Clusters_to_Dataset_distrib.png")

# %% Finding marker genes by subcluster
sc.pp.log1p(adata_CD4) # Calling this function before ranking genes avoids the need to use the raw data and fixes the dictionary issue
sc.tl.rank_genes_groups(adata_CD4, "leiden_Sub_Res_0.7", method="wilcoxon") # This function expects the data to be log-transformed

os.chdir("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/Plots/V4/Subclustering/CD4")
sc.pl.rank_genes_groups(adata_CD4, n_genes=25, sharey=False, save=".pdf")
sc.pl.rank_genes_groups(adata_CD4, n_genes=25, sharey=False, save=".png")

sc.pl.rank_genes_groups_stacked_violin(adata_CD4, n_genes=4, cmap='RdBu_r', groupby="leiden_Sub_Res_0.7", save="scvi.pdf")

# %% Marker genes by cluster --- Table
os.chdir("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/Subclustering_V4/Data/CD4")
result = adata_CD4.uns['rank_genes_groups']
groups = result['names'].dtype.names

Marker_genes_x_clus_df = pd.DataFrame(
    {group + '_' + key[:1]: result[key][group]
    for group in groups for key in ['names', 'pvals']})

Marker_genes_x_clus_df
Marker_genes_x_clus_df.to_csv("Marker_genes_per_cluster.csv")

# %% Gen exploration -- Paula's choice
os.chdir("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/Plots/V4/Subclustering/CD4")

sc.pl.stacked_violin(adata_CD4, var_names=["CD3D", "CD4", "CD8A"], groupby='leiden_Sub_Res_0.7', title="Marcadores basicos de celulas T", save="1_Basic_T_markers.pdf") # 1. Basic T cell markers
sc.pl.stacked_violin(adata_CD4, var_names=["FOXP3", "IL2RA", "IL17A"], groupby='leiden_Sub_Res_0.7', title="Otros marcadores", save="2_Other_markers.pdf") # 2. Other markers
sc.pl.stacked_violin(adata_CD4, var_names=["NKG7", "GNLY", "IFNG", "PRF1", "GZMK", "GZMB"], groupby='leiden_Sub_Res_0.7', title="Citotoxicidad", save="3_Cytotoxicity.pdf") # 3. Cytotoxicity markers
sc.pl.stacked_violin(adata_CD4, var_names=["CD69", "CD27", "CD28", "HLA-DRA", "TNFRSF9"], groupby='leiden_Sub_Res_0.7', title="Activacion", save="4_Activation.pdf") # 4. Activation markers
sc.pl.stacked_violin(adata_CD4, var_names=["TCF7", "CCR7", "SELL"], groupby='leiden_Sub_Res_0.7', title="Memoria", save="5_Memory.pdf") # 5. Memory markers
sc.pl.stacked_violin(adata_CD4, var_names=["PTPRC"], groupby='leiden_Sub_Res_0.7', title="CD45", save="6_CD45.pdf") # 6. CD45

sc.pl.violin(adata_CD4, keys=["CD3D", "CD4", "CD8A"], groupby='leiden_Sub_Res_0.7', stripplot=False, save="1_Basic_T_markers.pdf") # 1. Basic T cell markers
sc.pl.violin(adata_CD4, keys=["FOXP3", "IL2RA", "IL17A"], groupby='leiden_Sub_Res_0.7', stripplot=False, save="2_Other_markers.pdf") # 2. Other markers
sc.pl.violin(adata_CD4, keys=["NKG7", "GNLY", "IFNG", "PRF1", "GZMK", "GZMB"], groupby='leiden_Sub_Res_0.7', stripplot=False, save="3_Cytotoxicity.pdf") # 3. Cytotoxicity markers
sc.pl.violin(adata_CD4, keys=["CD69", "CD27", "CD28", "HLA-DRA", "TNFRSF9"], groupby='leiden_Sub_Res_0.7', stripplot=False, save="4_Activation.pdf") # 4. Activation markers
sc.pl.violin(adata_CD4, keys=["TCF7", "CCR7", "SELL"], groupby='leiden_Sub_Res_0.7', stripplot=False, save="5_Memory.pdf") # 5. Memory markers
sc.pl.violin(adata_CD4, keys=["PTPRC"], groupby='leiden_Sub_Res_0.7', stripplot=False, save="6_CD45.pdf") # 6. CD45

sc.pl.umap(adata_CD4, color=["CD3D", "CD4", "CD8A"], frameon=False, ncols=2, vmin=0, vmax="p99", save="1_Basic_T_markers.png") # 1. Basic T cell markers
sc.pl.umap(adata_CD4, color=["FOXP3", "IL2RA", "IL17A"], frameon=False, ncols=2, vmin=0, vmax="p99", save="2_Other_markers.png") # 2. Other markers
sc.pl.umap(adata_CD4, color=["NKG7", "GNLY", "IFNG", "PRF1", "GZMK", "GZMB"], frameon=False, ncols=2, vmin=0, vmax="p99", save="3_Cytotoxicity.png") # 3. Cytotoxicity markers
sc.pl.umap(adata_CD4, color=["CD69", "CD27", "CD28", "HLA-DRA", "TNFRSF9"], frameon=False, ncols=2, vmin=0, vmax="p99", save="4_Activation.png") # 4. Activation markers
sc.pl.umap(adata_CD4, color=["TCF7", "CCR7", "SELL"], frameon=False, ncols=2, vmin=0, vmax="p99", save="5_Memory.png") # 5. Memory markers
sc.pl.umap(adata_CD4, color=["PTPRC"], frameon=False, ncols=2, vmin=0, vmax="p99", save="6_CD45.png") # 6. CD45

# Extra 
sc.pl.umap(adata_CD4, color=["NCAM1"], frameon=False, ncols=2, vmin=0, vmax="p99", save="7_NCAM1.png") # 7. NCAM1

# %% Gen exploration -- CD4 subtypes exploration
os.chdir("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/Plots/V4/Subclustering/CD4")
sc.pl.stacked_violin(adata_CD4, var_names=["TBX21", "IFNG", "CXCR3"], groupby='leiden_Sub_Res_0.7', title="Th1", save="7_Th1.pdf") # 7. Th1
sc.pl.stacked_violin(adata_CD4, var_names=["GATA3", "IL4"], groupby='leiden_Sub_Res_0.7', title="Th2", save="8_Th2.pdf") # 8. Th2
sc.pl.stacked_violin(adata_CD4, var_names=["RORC", "IL17A", "CCR6"], groupby='leiden_Sub_Res_0.7', title="Th17", save="9_Th17.pdf") # 9. Th17
sc.pl.stacked_violin(adata_CD4, var_names=["IL9", "CCR3"], groupby='leiden_Sub_Res_0.7', title="Th9", save="10_Th9.pdf") # 10. Th9
sc.pl.stacked_violin(adata_CD4, var_names=["BCL6", "CXCR5", "IL21", "CD28"], groupby='leiden_Sub_Res_0.7', title="Tfh", save="11_Tfh.pdf") # 11. Tfh
sc.pl.stacked_violin(adata_CD4, var_names=["AHR", "IL22"], groupby='leiden_Sub_Res_0.7', title="Th22", save="12_Th22.pdf") # 12. Th22

sc.pl.violin(adata_CD4, keys=["TBX21", "IFNG", "CXCR3"], groupby='leiden_Sub_Res_0.7', stripplot=False, save="7_Th1.pdf") # 7. Th1
sc.pl.violin(adata_CD4, keys=["GATA3", "IL4"], groupby='leiden_Sub_Res_0.7', stripplot=False, save="8_Th2.pdf") # 8. Th2
sc.pl.violin(adata_CD4, keys=["RORC", "IL17A", "CCR6"], groupby='leiden_Sub_Res_0.7', stripplot=False, save="9_Th17.pdf") # 9. Th17
sc.pl.violin(adata_CD4, keys=["IL9", "CCR3"], groupby='leiden_Sub_Res_0.7', stripplot=False, save="10_Th9.pdf") # 10. Th9
sc.pl.violin(adata_CD4, keys=["BCL6", "CXCR5", "IL21", "CD28"], groupby='leiden_Sub_Res_0.7', stripplot=False, save="11_Tfh.pdf") # 11. Tfh
sc.pl.violin(adata_CD4, keys=["AHR", "IL22"], groupby='leiden_Sub_Res_0.7', stripplot=False, save="12_Th22.pdf") # 12. Th22

sc.pl.umap(adata_CD4, color=["TBX21", "IFNG", "CXCR3"], frameon=False, ncols=2, vmin=0, vmax="p99", save="7_Th1.png") # 7. Th1
sc.pl.umap(adata_CD4, color=["GATA3", "IL4"], frameon=False, ncols=2, vmin=0, vmax="p99", save="8_Th2.png") # 8. Th2
sc.pl.umap(adata_CD4, color=["RORC", "IL17A", "CCR6"], frameon=False, ncols=2, vmin=0, vmax="p99", save="9_Th17.png") # 9. Th17
sc.pl.umap(adata_CD4, color=["IL9", "CCR3"], frameon=False, ncols=2, vmin=0, vmax="p99", save="10_Th9.png") # 10. Th9
sc.pl.umap(adata_CD4, color=["BCL6", "CXCR5", "IL21", "CD28"], frameon=False, ncols=2, vmin=0, vmax="p99", save="11_Tfh.png") # 11. Tfh
sc.pl.umap(adata_CD4, color=["AHR", "IL22"], frameon=False, ncols=2, vmin=0, vmax="p99", save="12_Th22.png") # 12. Th22

# %% Code to visualize where a specific cluster is located
adata_CD4.obs['cluster_dummy9'] = adata_CD4.obs['leiden_Sub_Res_0.7'] == adata_CD4.obs['leiden_Sub_Res_0.7'].cat.categories[9]
adata_CD4.obs['cluster_dummy9'] = adata_CD4.obs['cluster_dummy9']*1
sc.pl.umap(adata_CD4, color='cluster_dummy9', frameon=False)

# %% More markers doubts
sc.pl.stacked_violin(adata_CD4, var_names=["PRF1", "IFNG", "GZMK", "GZMB", "CD8A", "NKG7"], groupby='leiden_Sub_Res_0.7', title="Marcadores dudas") # Cytotoxic (effector) function
sc.pl.stacked_violin(adata_CD4, var_names=["GNLY", "GZMA", "GZMM", "CCL5", "CCL3", "CCL4"], groupby='leiden_Sub_Res_0.7', title="Marcadores dudas") # Cytotoxic (effector) function 2
sc.pl.stacked_violin(adata_CD4, var_names=["CD69", "CD27", "CD28", "IL2RA"], groupby='leiden_Sub_Res_0.7', title="Marcadores dudas") # Activation markers
sc.pl.stacked_violin(adata_CD4, var_names=["IL10", "IL32"], groupby='leiden_Sub_Res_0.7', title="Marcadores dudas") # Cytokine secretion
sc.pl.stacked_violin(adata_CD4, var_names=["TCF7", "CCR7", "SELL", "IL7R"], groupby='leiden_Sub_Res_0.7', title="Marcadores dudas") # Memory / Memory-like
sc.pl.stacked_violin(adata_CD4, var_names=["ZWINT", "MKI67", "TUBA1B", "TUBB"], groupby='leiden_Sub_Res_0.7', title="Marcadores dudas") # Proliferation
sc.pl.stacked_violin(adata_CD4, var_names=["HAVCR2", "TIGIT", "LAG3", "CTLA4"], groupby='leiden_Sub_Res_0.7', title="Marcadores dudas") # Exhaustion
sc.pl.stacked_violin(adata_CD4, var_names=["KLRG1", "B3GAT1", "TIGIT", "CD160"], groupby='leiden_Sub_Res_0.7', title="Marcadores dudas") # Senescence
sc.pl.stacked_violin(adata_CD4, var_names=["SELL", "CXCR3", "CCR7", "CD44"], groupby='leiden_Sub_Res_0.7', title="Marcadores dudas") # Other memory markers
sc.pl.stacked_violin(adata_CD4, var_names=["TBX21", "B3GAT1", "SELL", "CD44"], groupby='leiden_Sub_Res_0.7', title="Marcadores dudas") # Other memory markers 2

# %% Gene Set Analysis -- Preprocessing steps
gene_set_names = gseapy.get_library_name(organism='Human')
print(gene_set_names)

os.chdir("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/Subclustering_V4/Data/CD4")
Marker_genes_x_clus_df = pd.read_csv("Marker_genes_per_cluster.csv", index_col=0)

# %% Gene Set Analysis -- Analysis for clusters
os.chdir("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/Plots/V4/Subclustering/CD4/GSA")

# Define a list of clusters and their corresponding gene lists
clusters = [
    {"name": "Cluster_0", "genes": Marker_genes_x_clus_df.iloc[:100, 0]},
    {"name": "Cluster_1", "genes": Marker_genes_x_clus_df.iloc[:100, 2]},
    {"name": "Cluster_2", "genes": Marker_genes_x_clus_df.iloc[:100, 4]},
    {"name": "Cluster_3", "genes": Marker_genes_x_clus_df.iloc[:100, 6]},
    {"name": "Cluster_4", "genes": Marker_genes_x_clus_df.iloc[:100, 8]},
    {"name": "Cluster_5", "genes": Marker_genes_x_clus_df.iloc[:100, 10]},
    {"name": "Cluster_6", "genes": Marker_genes_x_clus_df.iloc[:100, 12]},
    {"name": "Cluster_7", "genes": Marker_genes_x_clus_df.iloc[:100, 14]},
    {"name": "Cluster_8", "genes": Marker_genes_x_clus_df.iloc[:100, 16]},
    {"name": "Cluster_9", "genes": Marker_genes_x_clus_df.iloc[:100, 18]},
    {"name": "Cluster_10", "genes": Marker_genes_x_clus_df.iloc[:100, 20]}
]

# Define the gene set databases
gene_set_databases = {
    "GO_BP": "GO_Biological_Process_2023",
    "GO_MF": "GO_Molecular_Function_2023",
    "GO_CC": "GO_Cellular_Component_2023",
    "KEGG": "KEGG_2021_Human"
}

# Set the output directory
output_directory = "/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/Plots/V4/Subclustering/CD4/GSA"

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
    "0": "CD4 effector memory",
    "1": "CD4 central memory",
    "2": "CD8 effector memory",
    "3": "CD4 central memory",
    "4": "CD4 central memory",
    "5": "Regulatory T cells",
    "6": "CD4 central memory",
    "7": "CD8 memory",
    "8": "CD4 central memory",
    "9": "CD4 central memory",
    "10": "CD4 central memory",
}

adata_CD4.obs["manual_celltype_annotation_high"] = adata_CD4.obs["leiden_Sub_Res_0.7"].map(cl_annotation_high)

os.chdir("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/Plots/V4/Subclustering/CD4")
sc.pl.umap(adata_CD4, color = ["manual_celltype_annotation_high"], frameon=False, save="manual_celltype_annotation_high_CD4.png")

# %% Save the subclustering information
os.chdir("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/Subclustering_V4/Data/CD4")
adata_CD4.obs["manual_celltype_annotation_high"].to_csv("adata_CD4_manual_celltype_annotation_high.csv")

# %% End of script
