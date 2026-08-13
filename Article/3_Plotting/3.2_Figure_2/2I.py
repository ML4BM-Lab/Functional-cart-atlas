###############################################################################
###############################################################################

# Program: 2I.py
# Author: Sergio Cámara Peña
# Date: 06/03/2025
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
import scProportionTest as pt

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
figures_dir = project_dir / 'Resultados_Figuras' / 'Figura_2'
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

# Fixed color dict
palette_fija = dict(zip(categorias, palette))

########################################################################################################################
########################################################################################################################

# %% Figure 2I - Proportion Test CR Female vs CR Male
adata3 = adata_orig.copy()

## Filter the object
adata3 = adata3[adata3.obs["Antigen"] == "Blood"]
adata_Second = adata3.copy()

# Apply filters directly on adata_Second.obs
adata_Second = adata_Second[
    (adata_Second.obs["Time_Point_Ranges"] == "Infusion_Product") &
    (adata_Second.obs["Stimulated"] == "NO") &
    (adata_Second.obs["STATUS"] == "DISEASE")
].copy()

adata_Second = adata_Second[(adata_Second.obs["Max_Response"] == "CR")].copy()

# Create a mask for the age ranges to exclude
mask = (
    (adata_Second.obs["Age_Range"] == "20-40") |
    (adata_Second.obs["Age_Range"] == "<20")
)

# Apply the mask to filter cells
adata_Second = adata_Second[~mask].copy()

## Make and Save/Load the proportions test
results_dir = project_dir / 'Resultados' / 'Joined_datasets' / 'Post_Integration' / 'Data'

if True:
    results_1 = pt.permutation_test(adata_Second,"F","M",group_col='Sex',cell_type_col='manual_celltype_annotation_high',nperm=10000,alpha=0.05,n_bootstrap=10000,verbose=True)
    results_1.to_pickle(_output_path(results_dir, "F_CR_vs_M_CR.pkl"))
else:
    results_1 = pd.read_pickle(_input_path(results_dir, "F_CR_vs_M_CR.pkl"))

## Create plots
figures_dir = project_dir / 'Resultados_Figuras' / 'Figura_2'
figures_dir.mkdir(parents=True, exist_ok=True)
sc.settings.figdir = figures_dir
pt.point_range_plot(results_1, figsize=(10,8))
plot_1 = plt.gcf()
plot_1.savefig(
    _output_path(figures_dir, "2I_ProportionTest_F_CR_vs_M_CR.pdf"),
    bbox_inches="tight",
    dpi=300
)

# %% End of script
