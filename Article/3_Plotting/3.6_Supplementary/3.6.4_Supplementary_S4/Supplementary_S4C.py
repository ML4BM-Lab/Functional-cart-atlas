###############################################################################
###############################################################################

# Program: Supplementary_S4C.py
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
data_dir = project_dir / "Input"
if not data_dir.is_dir():
    raise FileNotFoundError(f"Required input directory does not exist: {data_dir}")
adata_orig = sc.read_h5ad(_input_path(data_dir, "Python_scVI_adata_big_V4_state4.h5ad"))
adata_orig_normalized = sc.read_h5ad(_input_path(data_dir, "Python_scVI_adata_big_V4_state4_Normalized.h5ad"))

# %% Set path to save figs
figures_dir = project_dir / "Figuras" / "Suplementarias"

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

# %% Figure S4C
# https://nbviewer.org/github/emdann/milopy/blob/master/notebooks/milopy_example.ipynb

adata2 = adata_orig_normalized.copy()
adata2 = adata2[adata2.obs["Antigen"] == "Blood"]

adata2 = adata2[
    (adata2.obs["Age_Range"] == "40-60") | (adata2.obs["Age_Range"] == ">60")
]

adata2 = adata2[(adata2.obs["Time_Point_Ranges"] == "Infusion_Product")]

adata2 = adata2[
    (adata2.obs["Max_Response"] == "CR") | (adata2.obs["Max_Response"] == "NR")
]

adata2 = adata2[(adata2.obs["Stimulated"] == "NO")]
print(adata2)

d = 30
k = 60

print(adata2.obs["Age_Range"].value_counts())
print(adata2.obs["Sex"].value_counts())
print(adata2.obs["Max_Response"].value_counts())

sc.pp.neighbors(adata2, use_rep="X_scVI", n_neighbors=k, n_pcs=d)

milo.make_nhoods(adata2)

nhood_size = np.array(adata2.obsm["nhoods"].sum(0)).ravel()

plt.figure(figsize=(8, 6))

plt.hist(nhood_size, bins=100)

plt.axvline(
    len(adata2.obs["Product_norm"].unique()) * 3,
    color="black",
    linestyle="dashed",
    linewidth=1,
)

plt.xlabel("Neighbourhood size")
plt.ylabel("Frequency")
plt.tight_layout()

plt.savefig(
    _output_path(figures_dir, "S4C_nhood_size_distribution.pdf"),
    format="pdf",
    bbox_inches="tight",
)

plt.show()

# %% End of script
