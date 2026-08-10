###############################################################################
###############################################################################

# Program: Celltypist_integ_plus_Jordanas.py
# Author: Sergio Cámara Peña
# Date: 30/10/2024
# Version: V5
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
_path_parser.add_argument(
    "--models-dir",
    default=None,
    help="CellTypist model directory (env: CART_ATLAS_MODELS_DIR; default: PROJECT_DIR/Models/CellTypist).",
)
_path_args, _unknown_path_args = _path_parser.parse_known_args()
project_dir = Path(_path_args.project_dir).expanduser().resolve()
if not project_dir.is_dir():
    raise FileNotFoundError(
        f"Project directory does not exist: {project_dir}. "
        "Set --project-dir or CART_ATLAS_PROJECT_DIR."
    )

models_dir = Path(
    _path_args.models_dir
    or os.environ.get("CART_ATLAS_MODELS_DIR", project_dir / "Models" / "CellTypist")
).expanduser().resolve()

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

# %% Import needed packages + Download and select celltypist model
import celltypist
from celltypist import models
import scanpy as sc
import os

Online=True
models.models_path = str(models_dir)
if Online:
    models_dir.mkdir(parents=True, exist_ok=True)
    models.download_models(force_update=True)

model_to_use = models.Model.load(model=str(_input_path(models_dir, "Immune_All_Low.pkl")))

####################################################################################################################
####################################################################################################################
####################################################################################################################

# %% Select (.csv) file on where to do the predictions
dir_merged = _require_path(project_dir / 'Resultados' / 'Joined_datasets' / 'Integration' / 'Python-Celltypist' / 'V5' / 'Seurat_merged.h5ad')
adata_Seurat_merged = sc.read(dir_merged)
sc.pp.normalize_per_cell(adata_Seurat_merged, counts_per_cell_after=10000)
sc.pp.log1p(adata_Seurat_merged)

# %% Predictions generation
predictions_merged = celltypist.annotate(
    filename=adata_Seurat_merged,
    model=model_to_use,
    majority_voting=True,
    transpose_input=True,
)
print(predictions_merged.predicted_labels)

# %% Choose destination path
adata = predictions_merged.to_adata()
_current_dir = Path(project_dir / 'Resultados' / 'Joined_datasets' / 'Integration' / 'Python-Celltypist' / 'V5')

# Save metadata table to a file
adata.obs.to_csv(_output_path(_current_dir, "celltypist_metadata_table.csv"), sep=",", index=True, header=True)

# %% End of script
