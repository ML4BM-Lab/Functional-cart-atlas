###############################################################################
###############################################################################

# Program: 5D.py
# Author: Sergio Cámara Peña
# Date: 04/08/2025
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
adata = sc.read(_input_path(data_dir, "Atlas_integ_scArches_FINAL_V5.h5ad"))

# %% Normalize and log data
# 1. Normalize total counts per cell
sc.pp.normalize_total(adata, target_sum=1e4)

# 2. Log-transform the data
sc.pp.log1p(adata)

# 3. Optionally, set raw to keep a copy of the normalized data
adata.raw = adata

# %% Filters:
# 1. Filter only cells where 'Antigen' is "Blood"
adata_filtered = adata[adata.obs['Antigen'] == "Blood"].copy()
print(adata_filtered.n_obs)

# 2. Filter to remove cells where 'Time_Point_Ranges' == "Infusion_Product" and 'Stimulated' == "YES"
mask = ~((adata_filtered.obs['Time_Point_Ranges'] == "Infusion_Product") & (adata_filtered.obs['Stimulated'] == "YES"))
adata_filtered = adata_filtered[mask].copy()
print(adata_filtered.n_obs)

# 3. Filter to remove cells with NA value in 'Max_Response'
adata_filtered = adata_filtered[~adata_filtered.obs['Max_Response'].isna()].copy()

# 4. Filter to keep only cells with 'Max_Response' == "CR"
adata_filtered = adata_filtered[adata_filtered.obs['Max_Response'] == "CR"].copy()
print(adata_filtered.n_obs)

# 5. Filter to keep only cells with 'STATUS' == "DISEASE" - NOT NEEDED
adata_filtered = adata_filtered[adata_filtered.obs['STATUS'] == "DISEASE"].copy()
print(adata_filtered.n_obs)

# 6. Filter to remove cells with 'Stimulation_Location' == "In_vitro" - NOT NEEDED
adata_filtered = adata_filtered[adata_filtered.obs['Stimulation_Location'] != "In_vitro"].copy()
print(adata_filtered.n_obs)

# 7. Filter to keep only cells with 'Time_Point_Ranges' == "2_weeks-3_months"
adata_filtered = adata_filtered[adata_filtered.obs['Time_Point_Ranges'] == "2_weeks-3_months"].copy()
print(adata_filtered.n_obs)

# %% Create and save plots
figures_dir = project_dir / 'Figuras' / 'Figura_5'
figures_dir.mkdir(parents=True, exist_ok=True)
sc.settings.figdir = figures_dir

Genes_exhaust = ["TIGIT", "CTLA4", "HAVCR2", "LAG3", "SPATA2"] # HAVCR2 = TIM3; SPATA2 = PD1

# Filter scFv == "BCMA"
adata_scFv1 = adata_filtered[adata_filtered.obs['ScFv'] == "BCMA"].copy()
sc.pl.stacked_violin(adata_scFv1, var_names=Genes_exhaust, groupby='Product_norm', use_raw=True, vmax=5, save="5D_BCMA.pdf")

# Filter scFv == "CD19"
adata_scFv2 = adata_filtered[adata_filtered.obs['ScFv'] == "CD19"].copy()
sc.pl.stacked_violin(adata_scFv2, var_names=Genes_exhaust, groupby='Product_norm', use_raw=True, vmax=5, save="5D_CD19.pdf")

# %% End of script
