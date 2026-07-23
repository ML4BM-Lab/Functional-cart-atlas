###############################################################################
###############################################################################

# Program: Supplementary_S1.py
# Author: Sergio Cámara Peña
# Date: 05/06/2025
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

from matplotlib.colors import to_rgb
import matplotlib.patches as mpatches
from collections import OrderedDict

# %% Set a random seed
random.seed(2504)

# %% Load data
data_dir = project_dir / 'Resultados' / 'Joined_datasets' / 'Raw_Atlas'
adata_orig = sc.read_h5ad(_input_path(data_dir, "Python_scVI_adata_big_V4_state4.h5ad"))

# %% Set path to save figs
figures_dir = project_dir / 'Resultados_Figuras' / 'Suplementarias'
figures_dir.mkdir(parents=True, exist_ok=True)
sc.settings.figdir = figures_dir

# %% Set palette
color_maps = {
    'STATUS': {'HEALTHY': '#a6d854', 'DISEASE': '#fc8d62'},
    'Time_Point_Ranges': {
        'Infusion_Product': '#fdae61',
        '<2_weeks': '#fee08b',
        '2_weeks-3_months': '#f46d43',
        '>3_months': '#d73027'
    },
    'Age_Range': {
    '<20':   '#a6d854',
    '20-40': '#66c2a5',
    '40-60': '#fc8d62',
    '>60':   '#e31a1c'
    },
    'Sex': {'F': '#e5c1ff', 'M': '#9ecae1'},
    'ScFv': {
        'CD19': '#66c2a5',
        'BCMA': '#fc8d62',
        'APRIL': '#8da0cb',
        'GD2': '#e78ac3',
        'HER2': '#a6d854'
    },
    'Max_Response': {
        'CR': '#006400',
        'PR': '#66c2a5',
        'NR': '#fc8d62'
    },
    'ICANS_Grade_Range': {
        '0': '#a6d854',
        '1-2': '#fdae61',
        '3-4': '#d73027'
    },
    'orig_ident': {
        'Bai_et_al': "#F8766D",
        'Boroughs_et_al': "#00BA38",
        'Lynn_et_al': "#FF9900",
        'Rodriguez_Marquez_et_al': "#619CFF",
        'Wang_et_al': "#F564E3",
        'Xhangolli_et_al': "#B79F00",
        'Deng_et_al': "#A81818",
        'Good_et_al': "#9D00FF",
        'Haradvala_et_al': "#006400",
        'Li_X_Cancer_Cell_letter_et_al': "#00BFC4",
        'Li_X_et_al': "#00008A",
        'Melenhorst_et_al': "#753900",
        'Sheih_et_al': "#EBFF0F"

}
}

# Desired order
status_order = ["DISEASE", "HEALTHY"]
time_order = [
    "Infusion_Product",
    "<2_weeks",
    "2_weeks-3_months",
    ">3_months"
]
age_order = [
    "<20", "20-40", "40-60", ">60"
]

response_order = [
    "NR", "PR", "CR"
]

# Transform into ordered categories
adata_orig.obs["STATUS"] = pd.Categorical(adata_orig.obs["STATUS"], categories=status_order, ordered=True)
adata_orig.obs["Time_Point_Ranges"] = pd.Categorical(adata_orig.obs["Time_Point_Ranges"], categories=time_order, ordered=True)
adata_orig.obs["Age_Range"] = pd.Categorical(adata_orig.obs["Age_Range"], categories=age_order, ordered=True)
adata_orig.obs["Max_Response"] = pd.Categorical(adata_orig.obs["Max_Response"], categories=response_order, ordered=True)

########################################################################################################################
########################################################################################################################

# %% Plot UMAPs
sc.pl.umap(
    adata_orig, 
    color=["orig_ident"],
    frameon=False,
    palette=color_maps["orig_ident"],
    save="_Suplementaria_S1_Orig_Ident.pdf")

sc.pl.umap(
    adata_orig, 
    color=["STATUS"],
    frameon=False,
    palette=color_maps["STATUS"],
    save="_Suplementaria_S1_STATUS.pdf")

sc.pl.umap(
    adata_orig, 
    color=["Time_Point_Ranges"],
    frameon=False,
    palette=color_maps["Time_Point_Ranges"],
    save="_Suplementaria_S1_Time_Point_Ranges.pdf")

sc.pl.umap(
    adata_orig, 
    color=["ScFv"],
    frameon=False,
    palette=color_maps["ScFv"],
    title="CAR Type",
    save="_Suplementaria_S1_CAR_Type.pdf")

sc.pl.umap(
    adata_orig, 
    color=["Max_Response"],
    frameon=False,
    palette=color_maps["Max_Response"],
    save="_Suplementaria_S1_Max_Response.pdf")

sc.pl.umap(
    adata_orig, 
    color=["ICANS_Grade_Range"],
    frameon=False,
    palette=color_maps["ICANS_Grade_Range"],
    save="_Suplementaria_S1_ICANS_Grade_Range.pdf")

sc.pl.umap(
    adata_orig, 
    color=["Age_Range"],
    frameon=False,
    palette=color_maps["Age_Range"],
    save="_Suplementaria_S1_Age_Range.pdf")

sc.pl.umap(
    adata_orig, 
    color=["Sex"],
    frameon=False,
    palette=color_maps["Sex"],
    save="_Suplementaria_S1_Sex.pdf")

# %% End of script