###############################################################################
###############################################################################

# Program: 4B.py
# Author: Sergio Cámara Peña
# Date: 29/05/2025
# Version: FINAL
# Machine: Rocinante

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

# %% Set a random seed
random.seed(2504)

# %% Read data
data_dir = project_dir / 'Input'
adata = sc.read(_input_path(data_dir, "Python_scVI_adata_big_V4_state4.h5ad"))

# %% PATH to save figs
figures_dir = project_dir / 'Figuras' / 'Figura_4'
figures_dir.mkdir(parents=True, exist_ok=True)
sc.settings.figdir = figures_dir

# %% Normalize and log data
# 1. Normalize total counts per cell
sc.pp.normalize_total(adata, target_sum=1e4)

# 2. Log-transform the data
sc.pp.log1p(adata)

# 3. Optionally, set raw to keep a copy of the normalized data
adata.raw = adata
adata_orig = adata.copy()

# %% Set genes to plot
Genes_IACs_vs_Rest = ["CD3D", "CD8A", "CD4", "CD14", "LYZ", "CXCL16", "CFD", "GRN", "TYROBP"]

# %% Figure 4Ba - Create the dummy column to separate IACs from the rest
adata.obs["IACs_dummy"] = adata.obs["manual_celltype_annotation_high"].apply(
    lambda x: "IACs" if x == "Monocyte-like T cells" else "Rest_of_Cells"
)

print(adata.obs["IACs_dummy"].value_counts())
print(adata.obs["manual_celltype_annotation_high"].value_counts())

# %% Figure 4Ba - Create and save plots
sc.pl.stacked_violin(adata, var_names=Genes_IACs_vs_Rest, groupby='IACs_dummy', title="", use_raw=True, save="4Ba.pdf", vmax=3)

# %% Figure 4Bb - Filter and only keep IACs
adata2 = adata_orig.copy()
adata2 = adata2[adata2.obs["manual_celltype_annotation_high"] == "Monocyte-like T cells"]
adata2

# %% Figure 4Bb - Create the dummy column to separate pre- and post-Infusion
adata2.obs["IACs_dummy"] = adata2.obs["Time_Point_Ranges"].apply(
    lambda x: "Infusion Product" if x == "Infusion_Product" else "Post Infusion"
)

print(adata2.obs["IACs_dummy"].value_counts())
print(adata2.obs["Time_Point_Ranges"].value_counts())

# %% Figure 4Bb - Create and save plots
sc.pl.stacked_violin(adata2, var_names=Genes_IACs_vs_Rest, groupby='IACs_dummy', title="", use_raw=True, save="4Bb.pdf", vmax=3)

# %% End of script
