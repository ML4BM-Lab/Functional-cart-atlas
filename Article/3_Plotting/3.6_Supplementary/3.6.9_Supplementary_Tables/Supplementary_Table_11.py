###############################################################################
###############################################################################

# Program: Supplementary_Table_11.py
# Author: Sergio Cámara Peña
# Date: 15/12/2025
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
import os
import random
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import mannwhitneyu

# %% Set a random seed
random.seed(2504)

# %% Read data
data_dir = project_dir / 'Input'
adata = sc.read(_input_path(data_dir, "Python_scVI_adata_big_V4_state4.h5ad"))

# %% Filter
adata = adata[adata.obs["Time_Point_Ranges"] == "Infusion_Product"]
adata = adata[adata.obs["STATUS"] == "DISEASE"]

print(adata.shape)

# %% Mark IACs per patient
# Total number of cells per patient
total_cells = adata.obs.groupby("Norm_Patient_Name").size()

# IACs cells per patient
iac_mask = adata.obs["manual_celltype_annotation_high"] == "Monocyte-like T cells"
iac_cells = adata.obs[iac_mask].groupby("Norm_Patient_Name").size()

# Put 0 where there is no IACs
iac_cells = iac_cells.reindex(total_cells.index, fill_value=0)

# Percentage of IACs per patient
perc_iac = (iac_cells / total_cells) * 100
perc_iac.name = "IAC_percentage"

# %% Build table patient-level with ICANS info
patient_info = (
    adata.obs[["Norm_Patient_Name", "ICANS_Grade_Range"]]
    .drop_duplicates()
    .set_index("Norm_Patient_Name")
    .join(perc_iac)
)

# Remove patients with out data of ICANS (NA) or ICANS grade == 0
patient_info = patient_info.dropna(subset=["ICANS_Grade_Range"])
patient_info = patient_info[~patient_info["ICANS_Grade_Range"].isin([0, "0"])]

print("\nTable patient-level:")
print(patient_info)

# %% Define groups of ICANS grade (1–2 vs 3–4)
icans_str = patient_info["ICANS_Grade_Range"].astype(str)

low = patient_info.loc[icans_str.isin(["1-2"]), "IAC_percentage"].dropna()
high = patient_info.loc[icans_str.isin(["3-4"]), "IAC_percentage"].dropna()

print("\nN patients ICANS 1-2:", len(low))
print("N patients ICANS 3-4:", len(high))

# %% Mann–Whitney tests (Wilcoxon rank-sum)
if len(low) > 0 and len(high) > 0:
    stat, p_wilcoxon = mannwhitneyu(low, high, alternative="two-sided")

    # Medians
    median_low = np.median(low)
    median_high = np.median(high)

    # Means
    mean_low = np.mean(low)
    mean_high = np.mean(high)

    # Value ranges
    range_low = (np.min(low), np.max(low))
    range_high = (np.min(high), np.max(high))

    # Fold-change (median and mean)
    fold_median = (median_high / median_low) if median_low > 0 else np.inf
    fold_mean = (mean_high / mean_low) if mean_low > 0 else np.inf

    print("\nWilcoxon rank-sum (IAC % per patient):")
    print("U-statistic:", stat)
    print("p-value:", p_wilcoxon)

    print("\nICANS 1-2:")
    print("  N =", len(low))
    print("  Median IAC%     =", median_low)
    print("  Mean IAC%       =", mean_low)
    print("  Range IAC%      =", range_low)

    print("\nICANS 3-4:")
    print("  N =", len(high))
    print("  Median IAC%     =", median_high)
    print("  Mean IAC%       =", mean_high)
    print("  Range IAC%      =", range_high)

    print("\nFold-change comparisons:")
    print("  Median fold-change (3-4 vs 1-2):", fold_median)
    print("  Mean fold-change (3-4 vs 1-2):", fold_mean)

# %% End of script
