###############################################################################
###############################################################################

# Program: Supplementary_Table_3.py
# Author: Sergio Cámara Peña
# Date: 09/06/2025
# Version: FINAL

###############################################################################
###############################################################################

#%% Load all the needed libraries
import scanpy as sc
import scvi
from scvi.model.utils import mde
import scanpy.external as sce
import os
import random
import pandas as pd
import seaborn as sns
from scgraph import scGraph

#%% Set a random seed
random.seed(2504)

#%% Set PATH to save
os.chdir("/home/scamara/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration_methods_lab/scGraph")

###############################################################################
###############################################################################

# %% Merged
# Initialize the graph analyzer
scgraph_merged = scGraph(
    adata_path="/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration_methods_lab/Merged_WO_integration/Sin_GT/Seurat_merged.h5ad",   # Path to AnnData object
    batch_key="Product_norm",                     # Column name for batch information
    label_key="manual_celltype_annotation_high",                 # Column name for cell type labels
    trim_rate=0.05,                        # Trim rate for robust mean calculation
    thres_batch=100,                       # Minimum number of cells per batch
    thres_celltype=10,                      # Minimum number of cells per cell type
    only_umap=True                        # Only evaluate 2D embeddings (mostly umaps)
)

# Run the analysis, return a pandas dataframe
results_merged = scgraph_merged.main()

# Save the results
results_merged.to_csv("embedding_evaluation_results_merged.csv")

# %% scVI
# Initialize the graph analyzer
scgraph_scVI = scGraph(
    adata_path="/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration_methods_lab/scVI/Sin_GT_With_Python/Suplementaria_S2_scVI_state1.h5ad",   # Path to AnnData object
    batch_key="Product_norm",                     # Column name for batch information
    label_key="manual_celltype_annotation_high",                 # Column name for cell type labels
    trim_rate=0.05,                        # Trim rate for robust mean calculation
    thres_batch=100,                       # Minimum number of cells per batch
    thres_celltype=10,                      # Minimum number of cells per cell type
    only_umap=True                        # Only evaluate 2D embeddings (mostly umaps)
)

# Run the analysis, return a pandas dataframe
results_scVI = scgraph_scVI.main()

# Save the results
results_scVI.to_csv("embedding_evaluation_results_scVI.csv")

###############################################################################
###############################################################################

# %% Harmony
# Initialize the graph analyzer
scgraph_Harmony = scGraph(
    adata_path="/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration_methods_lab/Harmony/Sin_GT/Seurat_harmony.h5ad",   # Path to AnnData object
    batch_key="Product_norm",                     # Column name for batch information
    label_key="manual_celltype_annotation_high",                 # Column name for cell type labels
    trim_rate=0.05,                        # Trim rate for robust mean calculation
    thres_batch=100,                       # Minimum number of cells per batch
    thres_celltype=10,                      # Minimum number of cells per cell type
    only_umap=True                        # Only evaluate 2D embeddings (mostly umaps)
)

# Run the analysis, return a pandas dataframe
results_Harmony = scgraph_Harmony.main()

# Save the results
results_Harmony.to_csv("embedding_evaluation_results_Harmony.csv")

###############################################################################
###############################################################################

# %% LIGER
# Initialize the graph analyzer
scgraph_LIGER = scGraph(
    adata_path="/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration_methods_lab/LIGER/Sin_GT/Seurat_liger.h5ad",   # Path to AnnData object
    batch_key="Product_norm",                     # Column name for batch information
    label_key="manual_celltype_annotation_high",                 # Column name for cell type labels
    trim_rate=0.05,                        # Trim rate for robust mean calculation
    thres_batch=100,                       # Minimum number of cells per batch
    thres_celltype=10,                      # Minimum number of cells per cell type
    only_umap=True                        # Only evaluate 2D embeddings (mostly umaps)
)

# Run the analysis, return a pandas dataframe
results_LIGER = scgraph_LIGER.main()

# Save the results
results_LIGER.to_csv("embedding_evaluation_results_LIGER.csv")

###############################################################################
###############################################################################

# %% STACAS
# Initialize the graph analyzer
scgraph_STACAS = scGraph(
    adata_path="/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration_methods_lab/STACAS/Sin_GT/Seurat_STACAS_integ2.h5ad",   # Path to AnnData object
    batch_key="Product_norm",                     # Column name for batch information
    label_key="manual_celltype_annotation_high",                 # Column name for cell type labels
    trim_rate=0.05,                        # Trim rate for robust mean calculation
    thres_batch=100,                       # Minimum number of cells per batch
    thres_celltype=10,                      # Minimum number of cells per cell type
    only_umap=True                        # Only evaluate 2D embeddings (mostly umaps)
)

# Run the analysis, return a pandas dataframe
results_STACAS = scgraph_STACAS.main()

# Save the results
results_STACAS.to_csv("embedding_evaluation_results_STACAS.csv")

###############################################################################
###############################################################################

# %% Seurat (RPCA)
# Initialize the graph analyzer
scgraph_RPCA = scGraph(
    adata_path="/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration_methods_lab/Seurat_RPCA/Sin_GT/Seurat_RPCA_integ.h5ad",   # Path to AnnData object
    batch_key="Product_norm",                     # Column name for batch information
    label_key="manual_celltype_annotation_high",                 # Column name for cell type labels
    trim_rate=0.05,                        # Trim rate for robust mean calculation
    thres_batch=100,                       # Minimum number of cells per batch
    thres_celltype=10,                      # Minimum number of cells per cell type
    only_umap=True                        # Only evaluate 2D embeddings (mostly umaps)
)

# Run the analysis, return a pandas dataframe
results_RPCA = scgraph_RPCA.main()

# Save the results
results_RPCA.to_csv("embedding_evaluation_results_RPCA.csv")

###############################################################################
###############################################################################

# %% fastMNN
# Initialize the graph analyzer
scgraph_fastMNN = scGraph(
    adata_path="/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration_methods_lab/fastMNN/Sin_GT/Seurat_fastMNN2.h5ad",   # Path to AnnData object
    batch_key="Product_norm",                     # Column name for batch information
    label_key="manual_celltype_annotation_high",                 # Column name for cell type labels
    trim_rate=0.05,                        # Trim rate for robust mean calculation
    thres_batch=100,                       # Minimum number of cells per batch
    thres_celltype=10,                      # Minimum number of cells per cell type
    only_umap=True                        # Only evaluate 2D embeddings (mostly umaps)
)

# Run the analysis, return a pandas dataframe
results_fastMNN = scgraph_fastMNN.main()

# Save the results
results_fastMNN.to_csv("embedding_evaluation_results_fastMNN.csv")

###############################################################################
###############################################################################

# %% End of script