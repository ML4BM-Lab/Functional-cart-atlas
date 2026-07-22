###############################################################################
###############################################################################

# Program: Supplementary_Table_9.py
# Author: Sergio Cámara Peña
# Date: 15/09/2025
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
    if not path.exists():
        raise FileNotFoundError(f"Required input path does not exist: {path}")
    return path


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
if not data_dir.is_dir():
    raise FileNotFoundError(f"Required input directory does not exist: {data_dir}")
adata = sc.read(_input_path(data_dir, "Python_scVI_adata_big_V4_state4.h5ad"))

# %% Normalize and log data
# 1. Normalize total counts per cell
sc.pp.normalize_total(adata, target_sum=1e4)

# 2. Log-transform the data
sc.pp.log1p(adata)

# 3. Optionally, set raw to keep a copy of the normalized data
adata.raw = adata

# %% Filter and only keep IACs
adata2 = adata.copy()
adata2 = adata2[adata2.obs["manual_celltype_annotation_high"] == "Monocyte-like T cells"]
adata2 = adata2[adata2.obs["STATUS"] == "DISEASE"]
adata2
del adata

# %% Create the dummy column to separate IACs from the rest
adata2.obs["IACs_dummy"] = adata2.obs["Time_Point_Ranges"].apply(
    lambda x: "Infusion Product" if x == "Infusion_Product" else "Post Infusion"
)

print(adata2.obs["IACs_dummy"].value_counts())
print(adata2.obs["Time_Point_Ranges"].value_counts())

Genes_IACs_clusters = ["AIF1", "LTA4H", "LST1"]

# %% Apply statistics
from scipy.stats import mannwhitneyu

def get_significance(p):
    if p < 0.001:
        return '***'
    elif p < 0.01:
        return '**'
    elif p < 0.05:
        return '*'
    else:
        return 'ns'

results = []

for gene in Genes_IACs_clusters:
    expr = adata2.raw[:, gene].X.toarray().flatten()
    groups = adata2.obs["IACs_dummy"]
    
    # Separate the groups
    group1 = expr[groups == "Infusion Product"]
    group2 = expr[groups == "Post Infusion"]
    
    # Mann-Whitney test
    stat, p = mannwhitneyu(group1, group2, alternative='two-sided')
    
    results.append({
        "Gene": gene,
        "U statistic": stat,
        "p-value": p
    })

# Convert into DataFrame for better visualization.
results_df = pd.DataFrame(results)
results_df["Significance"] = results_df["p-value"].apply(get_significance)
print(results_df)

# %% End of script