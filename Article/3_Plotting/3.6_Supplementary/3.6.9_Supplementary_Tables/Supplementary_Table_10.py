###############################################################################
###############################################################################

# Program: Supplementary_Table_10.py
# Author: Sergio Cámara Peña
# Date: 30/07/2025
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


# %% Load all the needed libraries
import scanpy as sc
import os
import random
import seaborn as sns
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
from scipy.stats.contingency import odds_ratio

# %% Set a random seed
random.seed(2504)

# %% Read data
data_dir = project_dir / 'Input'
if not data_dir.is_dir():
    raise FileNotFoundError(f"Required input directory does not exist: {data_dir}")
adata = sc.read(_input_path(data_dir, "Python_scVI_adata_big_V4_state4.h5ad"))

# %% Filter and mark IACs
adata2 = adata.copy()
adata2 = adata2[adata2.obs["Time_Point_Ranges"] == "Infusion_Product"]
adata2 = adata2[adata2.obs["STATUS"] == "DISEASE"]
print(adata2.shape)

adata2.obs['IACs'] = adata2.obs.groupby('Norm_Patient_Name')['manual_celltype_annotation_high'].transform(lambda x: 'Yes' if 'Monocyte-like T cells' in x.values else 'No')
print(adata2.shape)

del adata

# %% Contigency table
# Select the relevant columns
df = adata2.obs[["IACs", "ICANS_Grade_Range", "Norm_Patient_Name"]]

# Remove rows with NaN in 'ICANS_Grade_Range'
df = df.dropna(subset=['ICANS_Grade_Range'])

# Remove rows where 'ICANS_Grade_Range' is 0
df = df[~df['ICANS_Grade_Range'].isin([0, '0'])]

# Drop row names and keep only unique occurrences
df = df.drop_duplicates().reset_index(drop=True)

# Set 'Norm_Patient_Name' as the new index
df.set_index('Norm_Patient_Name', inplace=True)

# Create the contingency table
contingency_table = pd.crosstab(df["ICANS_Grade_Range"], df["IACs"])
print(contingency_table)

# %% Perform Fisher’s exact test + 95% CI for odds ratio

# Reorder table explicitly to avoid mistakes
# Rows: ICANS severity
# Columns: IAC status
contingency_table = contingency_table.reindex(
    index=["1-2", "3-4"],
    columns=["No", "Yes"],
    fill_value=0
)

print("\nContingency table:")
print(contingency_table)

# Build 2x2 table so that OR means:
# odds of severe ICANS in IAC+ patients vs IAC- patients
#
#              Severe ICANS   Mild ICANS
# IAC Yes            a             b
# IAC No             c             d

a = contingency_table.loc["3-4", "Yes"]
b = contingency_table.loc["1-2", "Yes"]
c = contingency_table.loc["3-4", "No"]
d = contingency_table.loc["1-2", "No"]

table_2x2 = np.array([
    [a, b],
    [c, d]
])

print("\n2x2 table used for Fisher test:")
print(table_2x2)

# Fisher exact test
oddsratio, p_value = fisher_exact(table_2x2)

# 95% CI for the odds ratio
or_result = odds_ratio(table_2x2, kind="sample")
ci = or_result.confidence_interval(confidence_level=0.95)

print("\nOdds ratio:", oddsratio)
print("95% CI:", ci.low, "-", ci.high)
print("p-value:", p_value)

# %% End of script
