###############################################################################
###############################################################################

# Program: Supplementary_Tables_10_11_Multiple_Testing_Correction.py
# Author: Sergio Cámara Peña
# Date: 23/07/2026
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


# %% Import libraries
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import fisher_exact, mannwhitneyu
from statsmodels.stats.multitest import multipletests

# %% Set input path
data_dir = project_dir / "Input"
input_file = _input_path(data_dir, "Python_scVI_adata_big_V4_state4.h5ad")

# %% Load and filter the data as in Supplementary Tables 10 and 11
adata = sc.read_h5ad(input_file)
obs = adata.obs.loc[
    (adata.obs["Time_Point_Ranges"] == "Infusion_Product")
    & (adata.obs["STATUS"] == "DISEASE"),
    [
        "Norm_Patient_Name",
        "manual_celltype_annotation_high",
        "ICANS_Grade_Range",
    ],
].copy()
del adata

# %% Build one validated row per patient
obs["IAC_cell"] = (
    obs["manual_celltype_annotation_high"] == "Monocyte-like T cells"
)

icans_per_patient = obs.groupby(
    "Norm_Patient_Name", observed=True
)["ICANS_Grade_Range"].apply(
    lambda values: pd.unique(values.dropna())
)
conflicting_icans = icans_per_patient[icans_per_patient.str.len() > 1]
if not conflicting_icans.empty:
    raise ValueError(
        "Patients with conflicting ICANS categories: "
        f"{conflicting_icans.index.tolist()}"
    )

patient_data = obs.groupby(
    "Norm_Patient_Name", observed=True
).agg(
    IAC_presence=("IAC_cell", lambda values: "Yes" if values.any() else "No"),
    IAC_percentage=("IAC_cell", lambda values: values.mean() * 100),
)
patient_data["ICANS_Grade_Range"] = icans_per_patient.apply(
    lambda values: values[0] if len(values) == 1 else np.nan
)

if len(patient_data) != 88:
    raise ValueError("Expected 88 unique patients.")

# Grade 0 and missing ICANS are excluded in both original analyses.
analysis_data = patient_data[
    patient_data["ICANS_Grade_Range"].isin(["1-2", "3-4"])
].copy()
if (
    len(analysis_data) != 28
    or (analysis_data["ICANS_Grade_Range"] == "1-2").sum() != 10
    or (analysis_data["ICANS_Grade_Range"] == "3-4").sum() != 18
):
    raise ValueError("The complete-case cohort does not match Tables 10 and 11.")

# %% Reproduce Supplementary Table 10: Fisher's exact test
contingency = pd.crosstab(
    analysis_data["ICANS_Grade_Range"],
    analysis_data["IAC_presence"],
).reindex(index=["1-2", "3-4"], columns=["No", "Yes"], fill_value=0)

table_2x2 = np.array(
    [
        [contingency.loc["3-4", "Yes"], contingency.loc["1-2", "Yes"]],
        [contingency.loc["3-4", "No"], contingency.loc["1-2", "No"]],
    ]
)
expected_table = np.array([[13, 3], [5, 7]])
if not np.array_equal(table_2x2, expected_table):
    raise ValueError("The Fisher 2x2 table does not match Supplementary Table 10.")

odds_ratio, fisher_p = fisher_exact(table_2x2, alternative="two-sided")

# %% Reproduce Supplementary Table 11: Mann-Whitney/Wilcoxon rank-sum test
low = analysis_data.loc[
    analysis_data["ICANS_Grade_Range"] == "1-2",
    "IAC_percentage",
]
high = analysis_data.loc[
    analysis_data["ICANS_Grade_Range"] == "3-4",
    "IAC_percentage",
]
u_statistic, wilcoxon_p = mannwhitneyu(
    low,
    high,
    alternative="two-sided",
)

# Check the unrounded results against the values reported in the manuscript.
if round(fisher_p, 4) != 0.0497 or round(wilcoxon_p, 3) != 0.028:
    raise ValueError(
        "The reproduced p-values do not match Supplementary Tables 10 and 11."
    )

# %% Apply Benjamini-Hochberg to the two clinical association tests
raw_p_values = np.array([fisher_p, wilcoxon_p])
reject, adjusted_p_values, _, _ = multipletests(
    raw_p_values,
    alpha=0.05,
    method="fdr_bh",
)

results = pd.DataFrame(
    {
        "Supplementary_Table": ["10", "11"],
        "Analysis": [
            "IAC presence vs ICANS severity",
            "IAC percentage vs ICANS severity",
        ],
        "Statistical_test": [
            "Fisher exact test",
            "Mann-Whitney U / Wilcoxon rank-sum test",
        ],
        "Test_statistic": [odds_ratio, u_statistic],
        "Raw_p_value": raw_p_values,
        "BH_adjusted_p_value": adjusted_p_values,
        "Reject_at_FDR_0.05": reject,
        "Multiple_testing_family": [
            "Two patient-level IAC-ICANS association tests",
            "Two patient-level IAC-ICANS association tests",
        ],
    }
)

# %% Report
print(results.to_string(index=False))

# %% End of script
