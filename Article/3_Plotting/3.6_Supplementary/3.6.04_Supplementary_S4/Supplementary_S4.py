###############################################################################
###############################################################################

# Program: Supplementary_S4.py
# Author: Sergio Cámara Peña
# Date: 10/06/2025
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

# %% Load all the needed libraries
import scanpy as sc
import os
import random
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import milopy.core as milo
import milopy.plot as milopl
import milopy.utils
import numpy as np
from plotnine import *

from matplotlib.colors import to_rgb
import matplotlib.patches as mpatches
from collections import OrderedDict

# %% Set a random seed
random.seed(2504)

# %% Load data
data_dir = project_dir / 'Resultados' / 'Joined_datasets' / 'Raw_Atlas'
adata_orig = sc.read_h5ad(_input_path(data_dir, "Python_scVI_adata_big_V4_state4.h5ad"))
adata_orig_normalized = sc.read_h5ad(_input_path(data_dir, "Python_scVI_adata_big_V4_state4_Normalized.h5ad"))

# %% Set path to save figs
figures_dir = project_dir / 'Resultados_Figuras' / 'Suplementarias'
figures_dir.mkdir(parents=True, exist_ok=True)
sc.settings.figdir = figures_dir

# %% Define a 12 color palette
palette = [
"#F8766D",
"#00BA38",
"#B79F00",
"#FF9900",
"#619CFF",
"#F564E3",
"#A81818",
"#9D00FF",
"#006400",
"#00BFC4",
"#00008A"
]

sns.palplot(palette)

categorias = [
    "Apoptotic T cells",
    "CD4 central memory",
    "CD4 cytotoxic",
    "CD4 effector memory",
    "CD8 cytotoxic",
    "CD8 effector memory",
    "CD8 memory",
    "Monocyte-like T cells",
    "Proliferative T cells",
    "Regulatory T cells",
    "Ribosomal enriched"
]

# Fixed colour dictionary
palette_fija = dict(zip(categorias, palette))

########################################################################################################################
########################################################################################################################

# %% S4A - UMAP Red/gray IP analysis
adata2 = adata_orig_normalized.copy()
adata2 = adata2[adata2.obs["Antigen"] == "Blood"]

adata2 = adata2[
    (adata2.obs["Age_Range"] == "40-60") | (adata2.obs["Age_Range"] == ">60")
]
adata2 = adata2[(adata2.obs["STATUS"] == "DISEASE")]
adata2 = adata2[(adata2.obs["Time_Point_Ranges"] == "Infusion_Product")]
adata2 = adata2[
    (adata2.obs["Max_Response"] == "CR") | (adata2.obs["Max_Response"] == "NR")
]
adata2 = adata2[(adata2.obs["Stimulated"] == "NO")]

rownames = adata2.obs.index.tolist()

adata1 = adata_orig_normalized.copy()
adata1.obs["dummy_variable"] = "No"
adata1.obs.loc[adata1.obs.index.isin(rownames), "dummy_variable"] = "Yes"

# Define colors: strong red for "Yes", grey for "No"
colors = ["grey", "red"]

# Make sure the categories are in the correct order
adata1.obs["dummy_variable"] = adata1.obs["dummy_variable"].astype("category")
adata1.obs["dummy_variable"].cat.set_categories(["No", "Yes"], inplace=True)

## Plot UMAP
# Get UMAP coordinates
umap = adata1.obsm["X_umap"]
df = adata1.obs.copy()
df["UMAP1"] = umap[:, 0]
df["UMAP2"] = umap[:, 1]

# Plot grey dots ("No") first
plt.figure(figsize=(6, 6))
df_no = df[df["dummy_variable"] == "No"]
plt.scatter(df_no["UMAP1"], df_no["UMAP2"], c="lightgrey", s=0.3, label="No", alpha=0.6)

# Plot red dots ("Yes") on top
df_yes = df[df["dummy_variable"] == "Yes"]
plt.scatter(df_yes["UMAP1"], df_yes["UMAP2"], c="darkred", s=0.3, label="Yes", alpha=0.9)

plt.tight_layout()
plt.savefig(_output_path(figures_dir, "Suplementaria_S4A.pdf"), dpi=300)
plt.show()

########################################################################################################################
########################################################################################################################

# %% S4B - UMAP of cell type composition of the red zone
sc.pl.umap(
    adata2, 
    color=["manual_celltype_annotation_high"],
    frameon=False,
    palette=palette_fija,
    save="Suplementaria_S4B.pdf"
)

########################################################################################################################
########################################################################################################################

# %% S4C - Neighborhood size distribution for Milo k selection

##### IN Supplementary_S4C.py script #####

########################################################################################################################
########################################################################################################################

# %% S4D - Splitted violins of each cell type impairment produced by IL10

##### IN Supplementary_S4D.R script #####

########################################################################################################################
########################################################################################################################

# %% S4E - scProportion test summary

##### IN Supplementary_S4E.R script #####

# %% End of script