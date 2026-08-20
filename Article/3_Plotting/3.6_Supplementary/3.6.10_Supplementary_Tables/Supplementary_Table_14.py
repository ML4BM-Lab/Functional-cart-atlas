###############################################################################
###############################################################################

# Program: Supplementary_Table_14.py
# Author: Sergio Cámara Peña
# Date: 20/08/2026
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

_generated_input_paths = {
    project_dir / "Resultados" / "Joined_datasets" / "Integration" / "Python-Celltypist" / "V5" / "Seurat_merged_With_Celltypist.h5ad",
    project_dir / "Resultados" / "Joined_datasets" / "Raw_Atlas" / "Atlas_integ_scArches_FINAL_V5.h5ad",
}

def _input_path(directory, filename):
    path = directory / filename
    if path.exists():
        return path
    if path in _generated_input_paths:
        input_path = project_dir / "Input" / path.name
        if input_path.exists():
            return input_path
        raise FileNotFoundError(
            "Required generated input file does not exist. Checked: "
            + ", ".join(str(candidate) for candidate in (path, input_path))
        )
    raise FileNotFoundError(f"Required input path does not exist: {path}")


def _output_path(directory, filename):
    directory.mkdir(parents=True, exist_ok=True)
    return directory / filename


# %% Load all the needed libraries
import scanpy as sc
import pandas as pd
import numpy as np

# %% Load Jordana object
jordana_dir = project_dir / "Resultados" / "Joined_datasets" / "Integration" / "Python-Celltypist" / "V5"

adata_Jordana = sc.read_h5ad(
    _input_path(jordana_dir, "Seurat_merged_With_Celltypist.h5ad")
)

adata_Jordana = adata_Jordana[
    adata_Jordana.obs["orig.ident"] == "Jordana_et_al"
].copy()


# %% Load final scArches object
atlas_dir = project_dir / "Resultados" / "Joined_datasets" / "Raw_Atlas"

adata_final = sc.read_h5ad(
    _input_path(atlas_dir, "Atlas_integ_scArches_FINAL_V5.h5ad")
)


# %% Match cells and reconstruct annotation comparison

common_cells = adata_Jordana.obs_names.intersection(
    adata_final.obs_names
)

comparison = pd.DataFrame(index=common_cells)

comparison["Lorea_Annotation"] = (
    adata_Jordana.obs.loc[
        common_cells,
        "Lorea_Annotation"
    ].astype(str)
)

comparison["scArches_Annotation"] = (
    adata_final.obs.loc[
        common_cells,
        "manual_celltype_annotation_high"
    ].astype(str)
)


# %% Convert Lorea annotation to scArches nomenclature

name_mapping = {
    "Proliferative": "Proliferative T cells",
    "CD8 effector memory": "CD8 effector memory",
    "CD4 activated": "CD4 effector memory",
    "CD4 memory": "CD4 central memory",
    "CD8 effector 1": "CD8 cytotoxic",
    "CD8 effector 2": "CD8 cytotoxic",
    "CD8 memory": "CD8 memory",
    "CD8 effectory GNLY-": "CD8 cytotoxic",
    "CD4 early memory": "CD4 central memory",
    "CD4 proliferative": "Proliferative T cells",
    "CD8 resting": "Unknown",
    "CD8 effector proliferative": "CD8 cytotoxic"
}

comparison["Jordanas_Equivalence"] = (
    comparison["Lorea_Annotation"].replace(name_mapping)
)


# %% LOREA EQUIVALENCE vs scArches

common_labels = (
    set(comparison["Jordanas_Equivalence"])
    &
    set(comparison["scArches_Annotation"])
)

comparison_filtered = comparison[
    comparison["Jordanas_Equivalence"].isin(common_labels)
    &
    comparison["scArches_Annotation"].isin(common_labels)
].copy()

accuracy_original = np.mean(
    comparison_filtered["Jordanas_Equivalence"].astype(str)
    ==
    comparison_filtered["scArches_Annotation"].astype(str)
)

print("\n")
print("=" * 70)
print("1. LOREA EQUIVALENCE vs scArches")
print("=" * 70)

print(
    f"\nAccuracy: "
    f"{accuracy_original:.6f}"
)

print(
    f"Accuracy (%): "
    f"{accuracy_original * 100:.2f}%"
)

confusion_counts = pd.crosstab(
    comparison_filtered["Jordanas_Equivalence"],
    comparison_filtered["scArches_Annotation"],
    rownames=["Jordanas_Equivalence"],
    colnames=["scArches_Annotation"]
)

print("\nConfusion matrix - cell counts")
print("-" * 70)
print(confusion_counts)

confusion_percentages = confusion_counts.div(
    confusion_counts.sum(axis=1),
    axis=0
) * 100

print("\nConfusion matrix - percentages")
print("-" * 70)
print(confusion_percentages.round(2))


# %% BROAD ACCURACY: CD4 vs CD8 vs PROLIFERATIVE

broad_comparison = comparison.copy()

broad_comparison["Jordanas_Broad"] = (
    broad_comparison["Lorea_Annotation"].replace({
        "Proliferative": "Proliferative",
        "CD4 proliferative": "Proliferative",

        "CD4 activated": "CD4",
        "CD4 memory": "CD4",
        "CD4 early memory": "CD4",

        "CD8 effector memory": "CD8",
        "CD8 effector 1": "CD8",
        "CD8 effector 2": "CD8",
        "CD8 memory": "CD8",
        "CD8 effectory GNLY-": "CD8",
        "CD8 resting": "CD8",
        "CD8 effector proliferative": "CD8"
    })
)

broad_comparison["scArches_Broad"] = (
    broad_comparison["scArches_Annotation"].replace({
        "CD4 central memory": "CD4",
        "CD8 cytotoxic": "CD8",
        "Proliferative T cells": "Proliferative"
    })
)

broad_accuracy = np.mean(
    broad_comparison["Jordanas_Broad"].astype(str)
    ==
    broad_comparison["scArches_Broad"].astype(str)
)

print("\n\n")
print("=" * 70)
print("2. BROAD ACCURACY: CD4 vs CD8 vs PROLIFERATIVE")
print("=" * 70)

print(
    f"\nAccuracy: "
    f"{broad_accuracy:.6f}"
)

print(
    f"Accuracy (%): "
    f"{broad_accuracy * 100:.2f}%"
)

broad_confusion_counts = pd.crosstab(
    broad_comparison["Jordanas_Broad"],
    broad_comparison["scArches_Broad"],
    rownames=["Jordanas_Broad"],
    colnames=["scArches_Broad"]
)

print("\nConfusion matrix - cell counts")
print("-" * 70)
print(broad_confusion_counts)

broad_confusion_percentages = broad_confusion_counts.div(
    broad_confusion_counts.sum(axis=1),
    axis=0
) * 100

print("\nConfusion matrix - percentages")
print("-" * 70)
print(broad_confusion_percentages.round(2))


# %% GLOBAL ACCURACY

global_accuracy = np.mean(
    comparison["Jordanas_Equivalence"].astype(str)
    ==
    comparison["scArches_Annotation"].astype(str)
)

print("\n\n")
print("=" * 70)
print("3. GLOBAL ACCURACY")
print("=" * 70)

print(
    f"\nAccuracy: "
    f"{global_accuracy:.6f}"
)

print(
    f"Accuracy (%): "
    f"{global_accuracy * 100:.2f}%"
)

global_confusion_counts = pd.crosstab(
    comparison["Jordanas_Equivalence"],
    comparison["scArches_Annotation"],
    rownames=["Jordanas_Equivalence"],
    colnames=["scArches_Annotation"]
)

print("\nConfusion matrix - cell counts")
print("-" * 70)
print(global_confusion_counts)

global_confusion_percentages = global_confusion_counts.div(
    global_confusion_counts.sum(axis=1),
    axis=0
) * 100

print("\nConfusion matrix - percentages")
print("-" * 70)
print(global_confusion_percentages.round(2))


# %% Save supplementary table outputs

output_dir = (
    project_dir
    / "Resultados_V5"
    / "Supplementary_Table_14"
)

def save_count_table(counts, accuracy, filename):
    output_table = counts.reset_index()
    output_table["Accuracy_percent"] = round(accuracy * 100, 2)
    output_table.to_csv(_output_path(output_dir, filename), index=False)


save_count_table(
    confusion_counts,
    accuracy_original,
    "Supplementary_Table_14_Lorea_equivalence_vs_scArches.csv"
)

save_count_table(
    broad_confusion_counts,
    broad_accuracy,
    "Supplementary_Table_14_broad_accuracy.csv"
)

save_count_table(
    global_confusion_counts,
    global_accuracy,
    "Supplementary_Table_14_global_accuracy.csv"
)


# %% End of script
