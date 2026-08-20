###############################################################################
###############################################################################

# Program: Supplementary_S8C,D,E.py
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


# %% Normalize gene expression

sc.pp.normalize_total(
    adata_final,
    target_sum=1e4
)

sc.pp.log1p(
    adata_final
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


# %% Lorea equivalence

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


# %% Broad annotation

comparison["Jordanas_Broad"] = comparison["Lorea_Annotation"].replace({
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

comparison["scArches_Broad"] = (
    comparison["scArches_Annotation"].replace({
        "CD4 central memory": "CD4",
        "CD8 cytotoxic": "CD8",
        "Proliferative T cells": "Proliferative"
    })
)


# %% Agreement

exact_correct = (
    comparison["Jordanas_Equivalence"]
    ==
    comparison["scArches_Annotation"]
)

broad_correct = (
    comparison["Jordanas_Broad"]
    ==
    comparison["scArches_Broad"]
)

comparison["Agreement"] = "Incorrect"

comparison.loc[
    broad_correct & ~exact_correct,
    "Agreement"
] = "Broad correct"

comparison.loc[
    exact_correct,
    "Agreement"
] = "Exact correct"


# %% Set marker genes

marker_genes = {
    "CD4 central memory": [
        "CD4",
        "TCF7",
        "CCR7",
        "SELL"
    ],

    "CD8 cytotoxic": [
        "CD8A",
        "NKG7",
        "GNLY",
        "IFNG",
        "PRF1",
        "GZMK",
        "GZMB"
    ],

    "Proliferative T cells": [
        "ZWINT",
        "MKI67",
        "TUBA1B",
        "TUBB"
    ]
}


# %% Set path to save figs

output_dir = (
    project_dir
    / "Resultados_Figuras"
    / "Suplementarias"
)

# %% Create stacked violins

group_order = [
    "Atlas",
    "Exact correct",
    "Broad correct",
    "Incorrect"
]

all_labels = adata_final.obs["manual_celltype_annotation_high"].astype(str)
is_jordana = adata_final.obs_names.isin(common_cells)

for celltype, markers in marker_genes.items():
    atlas_cells = adata_final.obs_names[(~is_jordana) & (all_labels == celltype)]
    exact_cells = comparison.index[(comparison["scArches_Annotation"] == celltype) & (comparison["Agreement"] == "Exact correct")]
    broad_cells = comparison.index[(comparison["scArches_Annotation"] == celltype) & (comparison["Agreement"] == "Broad correct")]
    incorrect_cells = comparison.index[(comparison["scArches_Annotation"] == celltype) & (comparison["Agreement"] == "Incorrect")]
    cells_to_plot = adata_final.obs_names.isin(atlas_cells) | adata_final.obs_names.isin(exact_cells) | adata_final.obs_names.isin(broad_cells) | adata_final.obs_names.isin(incorrect_cells)
    adata_plot = adata_final[cells_to_plot].copy()
    adata_plot.obs["Comparison_group"] = "Incorrect"
    adata_plot.obs.loc[adata_plot.obs_names.isin(atlas_cells), "Comparison_group"] = "Atlas"
    adata_plot.obs.loc[adata_plot.obs_names.isin(exact_cells), "Comparison_group"] = "Exact correct"
    adata_plot.obs.loc[adata_plot.obs_names.isin(broad_cells), "Comparison_group"] = "Broad correct"
    adata_plot.obs.loc[adata_plot.obs_names.isin(incorrect_cells), "Comparison_group"] = "Incorrect"
    present_groups = [group for group in group_order if (adata_plot.obs["Comparison_group"] == group).any()]
    adata_plot.obs["Comparison_group"] = pd.Categorical(adata_plot.obs["Comparison_group"], categories=present_groups, ordered=True)
    missing_genes = [gene for gene in markers if gene not in adata_plot.var_names]
    assert len(missing_genes) == 0, f"{celltype}: genes not found: {missing_genes}"
    print("\n" + "=" * 70)
    print(celltype)
    print("=" * 70)
    print(adata_plot.obs["Comparison_group"].value_counts(sort=False))
    plot = sc.pl.stacked_violin(adata_plot, var_names=markers, groupby="Comparison_group", categories_order=present_groups, use_raw=False, yticklabels=False, figsize=(8, 5), title=celltype, cmap="Reds", show=False, return_fig=True)
    safe_name = celltype.replace(" ", "_")
    plot.savefig(_output_path(output_dir, f"S8_Stacked_violin_{safe_name}.pdf"))
    plot.show()


# %% End of script
