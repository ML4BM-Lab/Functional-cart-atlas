###############################################################################
###############################################################################

# Program: 3F.py
# Author: Sergio Cámara Peña
# Date: 04/08/2025
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
    if path.exists():
        return path
    atlas_filenames = {
        "Python_scVI_adata_big_V4_state4.h5ad",
        "Python_scVI_adata_big_V4_state4_Normalized.h5ad",
        "Atlas_integ_scArches_FINAL_V5.h5ad",
    }
    atlas_directories = (
        project_dir / "Input",
        project_dir / "Resultados" / "Joined_datasets" / "Raw_Atlas",
    )
    if filename in atlas_filenames and directory in atlas_directories:
        alternate_paths = [
            candidate / filename for candidate in atlas_directories if candidate != directory
        ]
        for alternate_path in alternate_paths:
            if alternate_path.exists():
                return alternate_path
        checked_paths = [path, *alternate_paths]
        raise FileNotFoundError(
            "Required atlas input file does not exist. Checked: "
            + ", ".join(str(candidate) for candidate in checked_paths)
        )
    raise FileNotFoundError(f"Required input path does not exist: {path}")


def _output_path(directory, filename):
    directory.mkdir(parents=True, exist_ok=True)
    return directory / filename

# %% Import libraries
import os
import scanpy as sc
import anndata
import numpy as np
import random
import pandas as pd

# %% Set random seed
random.seed(2504)

# %% Load data
path = project_dir / 'Resultados' / 'Joined_datasets' / 'Raw_Atlas'
file = _input_path(path, "Python_scVI_adata_big_V4_state4.h5ad")
adata = anndata.read_h5ad(file)
print(adata.shape[0])

# %% Normalize
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Optionally, set raw to keep a copy of the normalized data
adata.raw = adata

# %% Filter object
# Antigen == "Blood"
adata_filtered = adata[adata.obs['Antigen'] == "Blood"].copy()
print(adata_filtered.shape[0])

# manual_celltype_annotation_high == "CD8 cytotoxic"
adata_filtered = adata_filtered[adata_filtered.obs['manual_celltype_annotation_high'] == "CD8 cytotoxic"].copy()
print(adata_filtered.shape[0])

# STATUS is not "Healthy"
adata_filtered = adata_filtered[adata_filtered.obs['STATUS'] != "HEALTHY"].copy()
print(adata_filtered.shape[0])

# Remove rows where Max_Response is NA
adata_filtered = adata_filtered[~adata_filtered.obs['Max_Response'].isna()].copy()

# Remove cells where Time_Point_Ranges == "Infusion_Product" and Stimulated == "YES"
mask = ~((adata_filtered.obs['Time_Point_Ranges'] == "Infusion_Product") & (adata_filtered.obs['Stimulated'] == "YES"))
adata_filtered = adata_filtered[mask].copy()
print(adata_filtered.shape[0])

adata_filtered_bis = adata_filtered.copy()

# %% Plot figure 3F
color_maps = {
    'Max_Response': {
        'CR': '#006400',
        'PR': '#66c2a5',
        'NR': '#fc8d62'
    }
}

# Convert to an ordered categorical factor
adata_filtered_bis.obs['Max_Response'] = pd.Categorical(
    adata_filtered_bis.obs['Max_Response'],
    categories=['CR', 'PR', 'NR'],
    ordered=True
)

# Visualize using the predefined color mapping
figures_dir = project_dir / 'Resultados_Figuras' / 'Figura_3'
figures_dir.mkdir(parents=True, exist_ok=True)
sc.settings.figdir = figures_dir
sc.pl.umap(
    adata_filtered_bis,
    color="Max_Response",
    palette=[color_maps['Max_Response'][cat] for cat in adata_filtered_bis.obs['Max_Response'].cat.categories],
    frameon=False,
    save="Figura_3F.pdf"
)

# %% End of script
