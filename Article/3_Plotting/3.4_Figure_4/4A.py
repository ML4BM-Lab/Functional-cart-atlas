###############################################################################
###############################################################################

# Program: 4A.py
# Author: Sergio Cámara Peña
# Date: 27/01/2025
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

def _require_path(path):
    if not path.exists():
        raise FileNotFoundError(f"Required input path does not exist: {path}")
    return path


def _input_path(directory, filename):
    path = directory / filename
    if not path.exists():
        raise FileNotFoundError(f"Required input path does not exist: {path}")
    return path


def _output_path(directory, filename):
    directory.mkdir(parents=True, exist_ok=True)
    return directory / filename

###############################################################################
###############################################################################

# %% Load all the needed libraries
import scanpy as sc
import os
import random
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import matplotlib.patches as mpatches

# %% Set a random seed
random.seed(2504)

# %% Read needed files
data_dir = project_dir / 'Resultados' / 'Joined_datasets' / 'Raw_Atlas'
if not data_dir.is_dir():
    raise FileNotFoundError(f"Required input directory does not exist: {data_dir}")
adata = sc.read_h5ad(_input_path(data_dir, "Python_scVI_adata_big_V4_state4.h5ad"))
print(adata.shape[0])  # Print the number of cells

# %% Filter object
filtered_adata = adata[adata.obs["Antigen"] == "Blood"].copy()
print(filtered_adata.shape[0])  # Print the number of cells

# Display the distribution of 'manual_celltype_annotation_high'
print(filtered_adata.obs["manual_celltype_annotation_high"].value_counts())

# %% Add a new column 'IACs_dummy' based on 'manual_celltype_annotation_high'
filtered_adata.obs["IACs_dummy"] = np.where(
    filtered_adata.obs["manual_celltype_annotation_high"] == "Monocyte-like T cells",
    "IACs",
    "Rest_of_Cells"
)

# Display the distribution of 'IACs_dummy'
print(filtered_adata.obs["IACs_dummy"].value_counts())

# %% Plot UMAP - Figure 4Aa
# Define colors for the categories
figures_dir = project_dir / 'Resultados_Figuras' / 'Figura_4'
figures_dir.mkdir(parents=True, exist_ok=True)
sc.settings.figdir = figures_dir
custom_colors = {
    "Rest_of_Cells": "lightgrey",
    "IACs": "#9D00FF"
}

# Map colors to the `IACs_dummy` column
filtered_adata.obs["IACs_color"] = filtered_adata.obs["IACs_dummy"].map(custom_colors)

# Figure
fig, ax = plt.subplots(figsize=(8, 8))
ax.scatter(
    filtered_adata.obsm["X_umap"][:, 0],  # UMAP 1
    filtered_adata.obsm["X_umap"][:, 1],  # UMAP 2
    c=filtered_adata.obs["IACs_color"],  # Color by IACs_dummy mapping
    s=10,  # Size of the points
    alpha=0.8  # Transparency
)
ax.axis("off")  # Remove axes
plt.tight_layout()
plt.savefig(_output_path(figures_dir, "Figura_4Aa.pdf"))
plt.show()

# %% Plot UMAP - Figure 4Ab
# Filter and only keep IACs
adata2 = adata.copy()
adata2 = adata2[adata2.obs["manual_celltype_annotation_high"] == "Monocyte-like T cells"]
adata2 = adata2[adata2.obs["STATUS"] == "DISEASE"]
filtered_adata_2 = adata2.copy()

# %% Add a new column 'IACs_dummy' based on 'Time_Point_Ranges'
filtered_adata_2.obs["IACs_dummy"] = np.where(
    filtered_adata_2.obs["Time_Point_Ranges"] == "Infusion_Product",
    "IP",
    "Post"
)

# Set working directory
figures_dir = project_dir / 'Resultados_Figuras' / 'Figura_4'
figures_dir.mkdir(parents=True, exist_ok=True)
sc.settings.figdir = figures_dir

# Define colors for the categories
custom_colors = {
    "IP": "#fdae61",
    "Post": "#2c7bb6"
}

# Map colors to the `IACs_dummy` column
filtered_adata_2.obs["IACs_color"] = filtered_adata_2.obs["IACs_dummy"].map(custom_colors)

# Figure
fig, ax = plt.subplots(figsize=(8, 8))
ax.scatter(
    filtered_adata_2.obsm["X_umap"][:, 0],  # UMAP 1
    filtered_adata_2.obsm["X_umap"][:, 1],  # UMAP 2
    c=filtered_adata_2.obs["IACs_color"],   # Color by IACs_dummy mapping
    s=10,
    alpha=0.8
)
ax.axis("off")  # Remove axes
plt.tight_layout()

# Add legend
legend_handles = [
    mpatches.Patch(color=color, label=label)
    for label, color in custom_colors.items()
]
ax.legend(handles=legend_handles, title="IACs", loc="lower left", frameon=False)

# Save and show
plt.savefig(_output_path(figures_dir, "Figura_4Ab.pdf"))
plt.show()

# %% End of script