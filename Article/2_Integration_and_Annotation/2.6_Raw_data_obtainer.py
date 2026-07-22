###############################################################################
###############################################################################

# Program: 2.6_Raw_data_obtainer.py
# Author: Sergio Cámara Peña
# Date: 12/11/2024
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
import pandas as pd

# %% Switches
Guardar = True

# %% Read files - V4
_current_dir = Path(project_dir / 'Resultados' / 'Joined_datasets' / 'Integration' / 'scVI' / 'V4')
adata_small_V4 = sc.read_h5ad(_input_path(_current_dir, "Python_scVI_adata_V4_state4.h5ad"))

_current_dir = Path(project_dir / 'Resultados' / 'Joined_datasets' / 'Integration' / 'Python-Celltypist' / 'V4')
adata_V4 = sc.read_h5ad(_input_path(_current_dir, "Seurat_merged_With_Celltypist.h5ad"))

# %% Remove irrelevant cells from original object
if Guardar:
    # Find the intersection of the two lists
    list1 = adata_small_V4.obs.index.copy()
    list2 = adata_V4.obs.index.copy()

    intersection = set(list1).intersection(set(list2))
    print(len(list1))
    print(len(list2))
    print(len(intersection))

    # Filter extra cells
    filter_mask = adata_V4.obs_names.isin(intersection)
    filtered_adata_V4 = adata_V4[filter_mask]

    # Find the intersection of the two lists
    list1 = adata_small_V4.obs.index.copy()
    list2 = filtered_adata_V4.obs.index.copy()

    intersection_2 = set(list1).intersection(set(list2))
    print(len(list1))
    print(len(list2))
    print(len(intersection_2))

# %% Create big object, with all relevant layers.
if Guardar:
    adata_big_V4 = filtered_adata_V4.copy()
    adata_big_V4.obs = adata_small_V4.obs.copy()
    adata_big_V4.obsm = adata_small_V4.obsm.copy()
    adata_big_V4.obsp = adata_small_V4.obsp.copy()
    del adata_small_V4.uns["log1p"]
    adata_big_V4.uns = adata_small_V4.uns.copy()
    adata_big_V4.__dict__['_raw'].__dict__['_var'] = adata_big_V4.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'}) # Fix an error while saving
    adata_big_V4.raw = None

# %% Change to adjust to final version of metadata
if Guardar:
    # Load CSV with changes
    _current_dir = Path(project_dir / 'Resultados' / 'Joined_datasets' / 'Raw_Atlas')
    time_point_changes = pd.read_csv(_input_path(_current_dir, "Time_Point_Changes_V4.csv"), sep=";", index_col="Norm_Sample_Name")

    # Join directly to update Time_Point and Time_Point_Ranges
    adata_big_V4.obs = adata_big_V4.obs.drop(columns=["Time_Point", "Time_Point_Ranges"], errors="ignore")
    adata_big_V4.obs = adata_big_V4.obs.join(
        time_point_changes,
        on="Product_norm"   # match Product_norm in obs with Norm_Sample_Name index
    )

    # Quick check
    pd.set_option("display.max_rows", None)
    Check = adata_big_V4.copy()
    print(Check.obs[["Product_norm", "Time_Point", "Time_Point_Ranges"]].drop_duplicates().to_string(index=False))

# %% Save files - V4
_current_dir = Path(project_dir / 'Resultados' / 'Joined_datasets' / 'Raw_Atlas')
if Guardar:
    adata_small_V4.write(_output_path(_current_dir, "Python_scVI_adata_small_V4_state4.h5ad"))
    adata_big_V4.write(_output_path(_current_dir, "Python_scVI_adata_big_V4_state4.h5ad"))
else:
    adata_big_V4 = sc.read_h5ad(_input_path(_current_dir, "Python_scVI_adata_big_V4_state4.h5ad"))

# %% Save files normalized - V4
_current_dir = Path(project_dir / 'Resultados' / 'Joined_datasets' / 'Raw_Atlas')
if Guardar:
    sc.pp.normalize_total(adata_big_V4, target_sum=1e4)
    sc.pp.log1p(adata_big_V4)
    adata_big_V4.raw = adata_big_V4.copy()
    adata_big_V4.write(_output_path(_current_dir, "Python_scVI_adata_big_V4_state4_Normalized.h5ad"))
else:
    adata_big_V4 = sc.read_h5ad(_input_path(_current_dir, "Python_scVI_adata_big_V4_state4_Normalized.h5ad"))

# %% End of script
