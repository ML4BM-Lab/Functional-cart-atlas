###############################################################################
###############################################################################

# Program: 1F.py
# Author: Sergio Cámara Peña
# Date: 07/01/2026
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

code_dir = _ARTICLE_DIR / '3_Plotting' / 'Functions'
if not code_dir.is_dir():
    raise FileNotFoundError(f"Required input directory does not exist: {code_dir}")

sys.path.insert(0, str(code_dir))
from scanpy_cluster_proportions import get_cluster_proportions, plot_cluster_proportions

# %% Set a random seed
random.seed(2504)

# %% Load data
data_dir = project_dir / 'Resultados' / 'Joined_datasets' / 'Raw_Atlas'
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

# %% Set path to save figs
figures_dir = project_dir / 'Resultados_Figuras' / 'Figura_1'
figures_dir.mkdir(parents=True, exist_ok=True)
sc.settings.figdir = figures_dir

########################################################################################################################
########################################################################################################################

# %% 1F - Stacked barplot with cell cycle phases
# Cluster proportion calculation and plotting
adata2 = adata_orig_normalized.copy()
metadata = adata2.obs.copy()

proportions = (
    metadata.groupby(
        ["ScFv", "manual_celltype_annotation_high"]
    )
    .size()
    .reset_index(name="count")
)
proportions["proportion"] = proportions.groupby(["ScFv"])[
    "count"
].transform(lambda x: x / x.sum())

proportions["manual_celltype_annotation_high"] = pd.Categorical(
    proportions["manual_celltype_annotation_high"],
    categories=list(palette_fija.keys()),
    ordered=True
)

scfv_order = ["BCMA", "CD19", "APRIL", "GD2", "HER2"]
proportions["ScFv"] = pd.Categorical(
    proportions["ScFv"],
    categories=scfv_order,
    ordered=True
)

plot = (
    ggplot(
        proportions,
        aes(x="ScFv", y="proportion", fill="manual_celltype_annotation_high"),
    )
    + geom_bar(stat="identity")
    + scale_fill_manual(values=palette_fija)
    + theme_classic()
    + ggtitle("")
    + xlab("")
    + ylab("Proportion")
    + theme(axis_text_x=element_text(rotation=90, ha='right'))
)

# Save it
ggsave(plot, _output_path(figures_dir, "1F_HEALTHY_Stacked_barplot_ScFv.pdf"))

########################################################################################################################
########################################################################################################################

# %% End of script
