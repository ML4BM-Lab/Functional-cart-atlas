###############################################################################
###############################################################################

# Program: 1C.py
# Author: Sergio Cámara Peña
# Date: 05/06/2025
# Version: FINAL
# Machine: Margaret

###############################################################################
###############################################################################

# %% Command-line and environment path configuration

import argparse
import os
import sys
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

code_dir = project_dir / 'Codigo' / 'Codigo_datasets_atlas' / 'Datasets_Integration' / 'Initial_Version_Atlas'
if not code_dir.is_dir():
    raise FileNotFoundError(f"Required input directory does not exist: {code_dir}")
sys.path.insert(0, str(code_dir))
from scanpy_cluster_proportions import get_cluster_proportions, plot_cluster_proportions


# %% Set a random seed
random.seed(2504)

# %% Load data
data_dir = project_dir / 'Resultados' / 'Joined_datasets' / 'Raw_Atlas'
if not data_dir.is_dir():
    raise FileNotFoundError(f"Required input directory does not exist: {data_dir}")
adata_orig = sc.read_h5ad(_input_path(data_dir, "Python_scVI_adata_big_V4_state4.h5ad"))
adata_orig_normalized = sc.read_h5ad(_input_path(data_dir, "Python_scVI_adata_big_V4_state4_Normalized.h5ad"))

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

# %% Set path to save figs
figures_dir = project_dir / 'Resultados_Figuras' / 'Figura_1'
figures_dir.mkdir(parents=True, exist_ok=True)
sc.settings.figdir = figures_dir

########################################################################################################################
########################################################################################################################

# %% Figure 1C - Dotplot of curated markers
adata2 = adata_orig_normalized.copy()
markers_V2 = ["CD3D", "CD4", "CD8A", "NKG7", "GNLY", "PRF1", "GZMK", "GZMB", "CD69", "CD27", "CD28", "HLA-DRA", "TNFRSF9", "TCF7", "CCR7", "SELL", "FOXP3", "IL2RA", "CTLA4", "IKZF2", "ZWINT", "MKI67", "TUBA1B", "TUBB", "CD68", "CD14", "LYZ"]
sc.pl.dotplot(adata2, markers_V2, use_raw=False, groupby='manual_celltype_annotation_high', dendrogram=True, figsize=(16, 4), save="1C_Annotation.pdf")

# %% End of script
