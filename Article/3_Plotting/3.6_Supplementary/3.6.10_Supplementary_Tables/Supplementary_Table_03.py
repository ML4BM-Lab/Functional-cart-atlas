###############################################################################
###############################################################################

# Program: Supplementary_Table_03.py
# Author: Sergio Cámara Peña
# Date: 09/06/2025
# Version: FINAL
# Machine: Margaret

###############################################################################
###############################################################################

# %% Command-line and environment path configuration

import argparse
import os
from pathlib import Path

_SCRIPT_DIR = Path(__file__).resolve().parent if "__file__" in globals() else Path.cwd()
_ARTICLE_DIR = next(
    (candidate for candidate in (_SCRIPT_DIR, *_SCRIPT_DIR.parents) if candidate.name == "Article"),
    _SCRIPT_DIR / "Article" if (_SCRIPT_DIR / "Article").is_dir() else _SCRIPT_DIR,
)
_path_parser = argparse.ArgumentParser()
_path_parser.add_argument(
    "--project-dir",
    default=os.environ.get("CART_ATLAS_PROJECT_DIR", str(_ARTICLE_DIR)),
    help="Root containing the repository-relative data, results, figure, model, and metadata subdirectories.",
)
_path_args, _unknown_path_args = _path_parser.parse_known_args()
project_dir = Path(_path_args.project_dir).expanduser().resolve()
if not project_dir.is_dir():
    raise FileNotFoundError(
        f"Project directory does not exist: {project_dir}. "
        "Set --project-dir or CART_ATLAS_PROJECT_DIR."
    )

_generated_input_paths = {
    project_dir / "Resultados" / "Joined_datasets" / "Integration_methods_lab" / "Merged_WO_integration" / "Sin_GT" / "Seurat_merged.h5ad",
    project_dir / "Resultados" / "Joined_datasets" / "Integration_methods_lab" / "scVI" / "Sin_GT_With_Python" / "Suplementaria_S2_scVI_state1.h5ad",
    project_dir / "Resultados" / "Joined_datasets" / "Integration_methods_lab" / "Harmony" / "Sin_GT" / "Seurat_harmony.h5ad",
    project_dir / "Resultados" / "Joined_datasets" / "Integration_methods_lab" / "LIGER" / "Sin_GT" / "Seurat_liger.h5ad",
    project_dir / "Resultados" / "Joined_datasets" / "Integration_methods_lab" / "STACAS" / "Sin_GT" / "Seurat_STACAS_integ2.h5ad",
    project_dir / "Resultados" / "Joined_datasets" / "Integration_methods_lab" / "Seurat_RPCA" / "Sin_GT" / "Seurat_RPCA_integ.h5ad",
    project_dir / "Resultados" / "Joined_datasets" / "Integration_methods_lab" / "fastMNN" / "Sin_GT" / "Seurat_fastMNN2.h5ad",
}

def _require_path(path):
    if path.exists():
        return path
    if path in _generated_input_paths:
        input_path = project_dir / "Input" / path.name
        if input_path.exists():
            return input_path
        raise FileNotFoundError(
            "Required generated input file does not exist. Checked: "
            + ", ".join(str(candidate) for candidate in (path, input_path))
        )
    raise FileNotFoundError(f"Required input path does not exist: {path}")


def _input_path(directory, filename):
    path = directory / filename
    if path.exists():
        return path
    if path in _generated_input_paths:
        input_path = project_dir / "Input" / path.name
        if input_path.exists():
            return input_path
        raise FileNotFoundError(
            "Required generated input file does not exist. Checked: "
            + ", ".join(str(candidate) for candidate in (path, input_path))
        )
    raise FileNotFoundError(f"Required input path does not exist: {path}")


def _output_path(directory, filename):
    directory.mkdir(parents=True, exist_ok=True)
    return directory / filename

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
results_dir = project_dir / 'Resultados' / 'Joined_datasets' / 'Integration_methods_lab' / 'scGraph'

###############################################################################
###############################################################################

# %% Merged
# Initialize the graph analyzer
scgraph_merged = scGraph(
    adata_path=_require_path(project_dir / 'Resultados' / 'Joined_datasets' / 'Integration_methods_lab' / 'Merged_WO_integration' / 'Sin_GT' / 'Seurat_merged.h5ad'),   # Path to AnnData object
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
results_merged.to_csv(_output_path(results_dir, "embedding_evaluation_results_merged.csv"))

# %% scVI
# Initialize the graph analyzer
scgraph_scVI = scGraph(
    adata_path=_require_path(project_dir / 'Resultados' / 'Joined_datasets' / 'Integration_methods_lab' / 'scVI' / 'Sin_GT_With_Python' / 'Suplementaria_S2_scVI_state1.h5ad'),   # Path to AnnData object
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
results_scVI.to_csv(_output_path(results_dir, "embedding_evaluation_results_scVI.csv"))

###############################################################################
###############################################################################

# %% Harmony
# Initialize the graph analyzer
scgraph_Harmony = scGraph(
    adata_path=_require_path(project_dir / 'Resultados' / 'Joined_datasets' / 'Integration_methods_lab' / 'Harmony' / 'Sin_GT' / 'Seurat_harmony.h5ad'),   # Path to AnnData object
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
results_Harmony.to_csv(_output_path(results_dir, "embedding_evaluation_results_Harmony.csv"))

###############################################################################
###############################################################################

# %% LIGER
# Initialize the graph analyzer
scgraph_LIGER = scGraph(
    adata_path=_require_path(project_dir / 'Resultados' / 'Joined_datasets' / 'Integration_methods_lab' / 'LIGER' / 'Sin_GT' / 'Seurat_liger.h5ad'),   # Path to AnnData object
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
results_LIGER.to_csv(_output_path(results_dir, "embedding_evaluation_results_LIGER.csv"))

###############################################################################
###############################################################################

# %% STACAS
# Initialize the graph analyzer
scgraph_STACAS = scGraph(
    adata_path=_require_path(project_dir / 'Resultados' / 'Joined_datasets' / 'Integration_methods_lab' / 'STACAS' / 'Sin_GT' / 'Seurat_STACAS_integ2.h5ad'),   # Path to AnnData object
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
results_STACAS.to_csv(_output_path(results_dir, "embedding_evaluation_results_STACAS.csv"))

###############################################################################
###############################################################################

# %% Seurat (RPCA)
# Initialize the graph analyzer
scgraph_RPCA = scGraph(
    adata_path=_require_path(project_dir / 'Resultados' / 'Joined_datasets' / 'Integration_methods_lab' / 'Seurat_RPCA' / 'Sin_GT' / 'Seurat_RPCA_integ.h5ad'),   # Path to AnnData object
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
results_RPCA.to_csv(_output_path(results_dir, "embedding_evaluation_results_RPCA.csv"))

###############################################################################
###############################################################################

# %% fastMNN
# Initialize the graph analyzer
scgraph_fastMNN = scGraph(
    adata_path=_require_path(project_dir / 'Resultados' / 'Joined_datasets' / 'Integration_methods_lab' / 'fastMNN' / 'Sin_GT' / 'Seurat_fastMNN2.h5ad'),   # Path to AnnData object
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
results_fastMNN.to_csv(_output_path(results_dir, "embedding_evaluation_results_fastMNN.csv"))

###############################################################################
###############################################################################

# %% End of script