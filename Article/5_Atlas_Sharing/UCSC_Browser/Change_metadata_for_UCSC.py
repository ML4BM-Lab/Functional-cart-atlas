###############################################################################
###############################################################################

# Program: Change_metadata_for_UCSC.py
# Author: Sergio Cámara Peña
# Date: 27/01/2026
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

# %% Import libraries
import os
import scanpy as sc

# %% Read files
_current_dir = Path(project_dir / 'Resultados' / 'Joined_datasets' / 'Raw_Atlas')
adata_V5 = sc.read_h5ad(_input_path(_current_dir, "Atlas_integ_scArches_FINAL_V5_Normalized.h5ad"))
adata_V5_orig = adata_V5.copy()

# %% Change metadata columns
adata_V5.obs.rename(columns={'ScFv': 'CAR_Type'}, inplace=True)

# %% Remove raw and genes that have few counts
adata_V5.raw = None
sc.pp.filter_genes(adata_V5, min_counts=20)
sc.pp.filter_genes(adata_V5, min_cells=20)

# %% Remove extra columns
cols_to_remove = ["Marked_cells", "barcode", "timepoint", "barcode_timepoint", "celltypist_predicted_labels", "celltypist_over_clustering", "celltypist_majority_voting", "celltypist_conf_score", "Product", "_scvi_batch", "_scvi_labels", "leiden"]
adata_V5.obs.drop(columns=cols_to_remove, inplace=True)

keep_uns = ["log1p", "umap"]

for key in list(adata_V5.uns.keys()):
    if key not in keep_uns:
        del adata_V5.uns[key]

for key in ["_scvi_extra_categorical_covs"]:
    if key in adata_V5.obsm:
        del adata_V5.obsm[key]

for key in ["X_pca"]:
    if key in adata_V5.obsm:
        del adata_V5.obsm[key]

adata_V5.var.drop(columns=["n_counts", "n_cells"], inplace=True)

for key in ["distances", "connectivities"]:
    if key in adata_V5.obsp:
        del adata_V5.obsp[key]

# %% Save files
_current_dir = Path(project_dir / 'UCSC_Browser')
adata_V5.write(_output_path(_current_dir, "Atlas_integ_scArches_FINAL_V5_Normalized_for_UCSC.h5ad"))

# %% End of script