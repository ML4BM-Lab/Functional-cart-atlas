###############################################################################
###############################################################################

# Program: Supplementary_S9.py
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
import matplotlib.pyplot as plt

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


# %% Set path to save figs
figures_dir = (
    project_dir
    / "Resultados_Figuras"
    / "Suplementarias"
)


# %% UMAPs: Atlas vs Jordana by scArches cell type

print("\n\n")
print("=" * 70)
print("4. UMAPs: Atlas vs Jordana by scArches cell type")
print("=" * 70)

if "X_umap" not in adata_final.obsm:
    raise ValueError("The UMAP is not present in adata_final.obsm['X_umap'].")

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

color_dict = dict(zip(categorias, palette))

target_labels = [
    cat for cat in categorias
    if cat in set(comparison["scArches_Annotation"])
]

umap = adata_final.obsm["X_umap"]

all_labels = (
    adata_final.obs[
        "manual_celltype_annotation_high"
    ].astype(str)
)

is_jordana = adata_final.obs_names.isin(
    common_cells
)

background_color = "#D9D9D9"
jordana_color = "#3B0D0C"

print("\nUMAPs will be saved in:")
print(figures_dir)

for celltype in target_labels:
    mask_type = all_labels == celltype
    mask_old = mask_type & (~is_jordana)
    mask_jordana = mask_type & is_jordana
    n_old = int(mask_old.sum())
    n_jordana = int(mask_jordana.sum())
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.scatter(umap[:, 0], umap[:, 1], s=4, c=background_color, alpha=0.45, linewidths=0)
    ax.scatter(umap[mask_old, 0], umap[mask_old, 1], s=7, c=color_dict[celltype], alpha=0.90, linewidths=0, label=f"Atlas: {celltype} (n={n_old})")
    ax.scatter(umap[mask_jordana, 0], umap[mask_jordana, 1], s=7, c=jordana_color, alpha=0.95, linewidths=0, label=f"Jordana: {celltype} (n={n_jordana})")
    ax.set_aspect("equal", adjustable="box")
    ax.set_title(celltype, fontsize=14)
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(1.02, 1), borderaxespad=0, markerscale=2.5)
    fig.tight_layout()
    pdf_name = f"S9_UMAP_{celltype.replace(' ', '_')}.pdf"
    tiff_name = f"S9_UMAP_{celltype.replace(' ', '_')}.tiff"
    fig.savefig(_output_path(figures_dir, pdf_name), bbox_inches="tight")
    fig.savefig(_output_path(figures_dir, tiff_name), bbox_inches="tight", dpi=300)
    print(f"Saved: {pdf_name}")
    print(f"Saved: {tiff_name}")
    plt.show()
    plt.close(fig)


# %% UMAPs: Annotation agreement between Lorea and scArches

print("\n\n")
print("=" * 70)
print("5. UMAPs: Annotation agreement between Lorea and scArches")
print("=" * 70)

exact_correct = (
    broad_comparison["Jordanas_Equivalence"].astype(str)
    ==
    broad_comparison["scArches_Annotation"].astype(str)
)

broad_correct = (
    broad_comparison["Jordanas_Broad"].astype(str)
    ==
    broad_comparison["scArches_Broad"].astype(str)
)

broad_comparison["Agreement"] = "Incorrect"

broad_comparison.loc[
    broad_correct,
    "Agreement"
] = "Broad correct"

broad_comparison.loc[
    exact_correct,
    "Agreement"
] = "Exact correct"

background_color = "#D9D9D9"
exact_color = "#7F0000"
broad_color = "#E34A33"
incorrect_color = "#FCAE91"

agreement_celltypes = [
    "CD4 central memory",
    "CD8 cytotoxic",
    "Proliferative T cells"
]

print("\nUMAPs will be saved in:")
print(figures_dir)

for celltype in agreement_celltypes:
    exact_cells = broad_comparison.index[(broad_comparison["scArches_Annotation"] == celltype) & (broad_comparison["Agreement"] == "Exact correct")]
    broad_cells = broad_comparison.index[(broad_comparison["scArches_Annotation"] == celltype) & (broad_comparison["Agreement"] == "Broad correct")]
    incorrect_cells = broad_comparison.index[(broad_comparison["scArches_Annotation"] == celltype) & (broad_comparison["Agreement"] == "Incorrect")]
    mask_exact = adata_final.obs_names.isin(exact_cells)
    mask_broad = adata_final.obs_names.isin(broad_cells)
    mask_incorrect = adata_final.obs_names.isin(incorrect_cells)
    n_exact = int(mask_exact.sum())
    n_broad = int(mask_broad.sum())
    n_incorrect = int(mask_incorrect.sum())
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.scatter(umap[:, 0], umap[:, 1], s=4, c=background_color, alpha=0.40, linewidths=0)
    ax.scatter(umap[mask_incorrect, 0], umap[mask_incorrect, 1], s=7, c=incorrect_color, alpha=0.90, linewidths=0, label=f"Incorrect (n={n_incorrect})")
    ax.scatter(umap[mask_broad, 0], umap[mask_broad, 1], s=7, c=broad_color, alpha=0.95, linewidths=0, label=f"Broad correct (n={n_broad})")
    ax.scatter(umap[mask_exact, 0], umap[mask_exact, 1], s=7, c=exact_color, alpha=1.0, linewidths=0, label=f"Exact correct (n={n_exact})")
    ax.set_aspect("equal", adjustable="box")
    ax.set_title(celltype, fontsize=14)
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(1.02, 1), borderaxespad=0, markerscale=2.5)
    fig.tight_layout()
    pdf_name = f"S9_UMAP_{celltype.replace(' ', '_')}_annotation_agreement.pdf"
    tiff_name = f"S9_UMAP_{celltype.replace(' ', '_')}_annotation_agreement.tiff"
    fig.savefig(_output_path(figures_dir, pdf_name), bbox_inches="tight")
    fig.savefig(_output_path(figures_dir, tiff_name), bbox_inches="tight", dpi=300)
    print(f"Saved: {pdf_name}")
    print(f"Saved: {tiff_name}")
    plt.show()
    plt.close(fig)


# %% End of script
