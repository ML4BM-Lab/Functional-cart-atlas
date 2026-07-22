###############################################################################
###############################################################################

# Program: scProportiontest_S4E.py
# Author: Sergio Cámara Peña
# Date: 30/09/2025
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
if not data_dir.is_dir():
    raise FileNotFoundError(f"Required input directory does not exist: {data_dir}")
adata_orig = sc.read_h5ad(_input_path(data_dir, "Python_scVI_adata_big_V4_state4.h5ad"))

# %% Set path to save figs
figures_dir = project_dir / 'Resultados_Figuras' / 'Data'
figures_dir.mkdir(parents=True, exist_ok=True)
sc.settings.figdir = figures_dir

# %% Set switches
First_time = True

########################################################################################################################
########################################################################################################################

# %% Test men >60 yo
adata1 = adata_orig.copy()

## Filter the object
adata1 = adata1[adata1.obs["Antigen"] == "Blood"]
adata_First = adata1.copy()

# Apply filters directly on adata_First.obs
adata_First = adata_First[
    (adata_First.obs["Time_Point_Ranges"] == "Infusion_Product") &
    (adata_First.obs["Stimulated"] == "NO") &
    (adata_First.obs["STATUS"] == "DISEASE")
].copy()

adata_First = adata_First[(adata_First.obs["Sex"] == "M")].copy()

# Create a mask for the age ranges to include
mask = (
    (adata_First.obs["Age_Range"] == ">60")
)

# Apply the mask to filter cells
adata_First = adata_First[mask].copy()

## Make and Save/Load the proportions test
if First_time:
    results_1 = pt.permutation_test(adata_First,"NR","CR",group_col='Max_Response',cell_type_col='manual_celltype_annotation_high',nperm=10000,alpha=0.05,n_bootstrap=10000,verbose=True)
    results_1.to_pickle(_output_path(figures_dir, "mo.pkl"))
    results_1.to_csv(_output_path(figures_dir, "mo.csv"))
else:
    results_1 = pd.read_pickle(_input_path(figures_dir, "mo.pkl"))

########################################################################################################################
########################################################################################################################

# %% Test women >60 yo
adata1 = adata_orig.copy()

## Filter the object
adata1 = adata1[adata1.obs["Antigen"] == "Blood"]
adata_Second = adata1.copy()

# Apply filters directly on adata_Second.obs
adata_Second = adata_Second[
    (adata_Second.obs["Time_Point_Ranges"] == "Infusion_Product") &
    (adata_Second.obs["Stimulated"] == "NO") &
    (adata_Second.obs["STATUS"] == "DISEASE")
].copy()

adata_Second = adata_Second[(adata_Second.obs["Sex"] == "F")].copy()

# Create a mask for the age ranges to include
mask = (
    (adata_Second.obs["Age_Range"] == ">60")
)

# Apply the mask to filter cells
adata_Second = adata_Second[mask].copy()

## Make and Save/Load the proportions test
if First_time:
    results_2 = pt.permutation_test(adata_Second,"NR","CR",group_col='Max_Response',cell_type_col='manual_celltype_annotation_high',nperm=10000,alpha=0.05,n_bootstrap=10000,verbose=True)
    results_2.to_pickle(_output_path(figures_dir, "wo.pkl"))
    results_2.to_csv(_output_path(figures_dir, "wo.csv"))
else:
    results_2 = pd.read_pickle(_input_path(figures_dir, "wo.pkl"))

########################################################################################################################
########################################################################################################################

# %% Test men 40-60 yo
adata1 = adata_orig.copy()

## Filter the object
adata1 = adata1[adata1.obs["Antigen"] == "Blood"]
adata_Third = adata1.copy()

# Apply filters directly on adata_Third.obs
adata_Third = adata_Third[
    (adata_Third.obs["Time_Point_Ranges"] == "Infusion_Product") &
    (adata_Third.obs["Stimulated"] == "NO") &
    (adata_Third.obs["STATUS"] == "DISEASE")
].copy()

adata_Third = adata_Third[(adata_Third.obs["Sex"] == "M")].copy()

# Create a mask for the age ranges to include
mask = (
    (adata_Third.obs["Age_Range"] == "40-60")
)

# Apply the mask to filter cells
adata_Third = adata_Third[mask].copy()

## Make and Save/Load the proportions test
if First_time:
    results_3 = pt.permutation_test(adata_Third,"NR","CR",group_col='Max_Response',cell_type_col='manual_celltype_annotation_high',nperm=10000,alpha=0.05,n_bootstrap=10000,verbose=True)
    results_3.to_pickle(_output_path(figures_dir, "my.pkl"))
    results_3.to_csv(_output_path(figures_dir, "my.csv"))
else:
    results_3 = pd.read_pickle(_input_path(figures_dir, "my.pkl"))

########################################################################################################################
########################################################################################################################

# %% Test women 40-60 yo
adata1 = adata_orig.copy()

## Filter the object
adata1 = adata1[adata1.obs["Antigen"] == "Blood"]
adata_Fourth = adata1.copy()

# Apply filters directly on adata_Fourth.obs
adata_Fourth = adata_Fourth[
    (adata_Fourth.obs["Time_Point_Ranges"] == "Infusion_Product") &
    (adata_Fourth.obs["Stimulated"] == "NO") &
    (adata_Fourth.obs["STATUS"] == "DISEASE")
].copy()

adata_Fourth = adata_Fourth[(adata_Fourth.obs["Sex"] == "F")].copy()

# Create a mask for the age ranges to include
mask = (
    (adata_Fourth.obs["Age_Range"] == "40-60")
)

# Apply the mask to filter cells
adata_Fourth = adata_Fourth[mask].copy()

## Make and Save/Load the proportions test
if First_time:
    results_4 = pt.permutation_test(adata_Fourth,"NR","CR",group_col='Max_Response',cell_type_col='manual_celltype_annotation_high',nperm=10000,alpha=0.05,n_bootstrap=10000,verbose=True)
    results_4.to_pickle(_output_path(figures_dir, "wy.pkl"))
    results_4.to_csv(_output_path(figures_dir, "wy.csv"))
else:
    results_4 = pd.read_pickle(_input_path(figures_dir, "wy.pkl"))

# %% End of script