###############################################################################
###############################################################################

# Program: Supplementary_Table_4.py
# Author: Sergio Cámara Peña
# Date: 09/06/2026
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
    project_dir / "Resultados" / "Joined_datasets" / "Integration_methods_lab" / "Merged_WO_integration" / "Sin_GT" / "Seurat_merged.h5ad",
    project_dir / "Resultados" / "Joined_datasets" / "Integration_methods_lab" / "scVI" / "Sin_GT_With_Python" / "Suplementaria_S2_scVI_state1.h5ad",
    project_dir / "Resultados" / "Joined_datasets" / "Integration_methods_lab" / "Harmony" / "Sin_GT" / "Seurat_harmony.h5ad",
    project_dir / "Resultados" / "Joined_datasets" / "Integration_methods_lab" / "LIGER" / "Sin_GT" / "Seurat_liger.h5ad",
    project_dir / "Resultados" / "Joined_datasets" / "Integration_methods_lab" / "STACAS" / "Sin_GT" / "Seurat_STACAS_integ2.h5ad",
    project_dir / "Resultados" / "Joined_datasets" / "Integration_methods_lab" / "Seurat_RPCA" / "Sin_GT" / "Seurat_RPCA_integ.h5ad",
    project_dir / "Resultados" / "Joined_datasets" / "Integration_methods_lab" / "fastMNN" / "Sin_GT" / "Seurat_fastMNN2.h5ad",
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


#%% Load all the needed libraries
import scanpy as sc
import os
import random
import numpy as np
import pandas as pd
from scib_metrics.benchmark import Benchmarker, BioConservation, BatchCorrection
from rich import print

###############################################################################
###############################################################################

#%% Set seed
os.environ["PYTHONHASHSEED"] = "2504"
random.seed(2504)
np.random.seed(2504)

###############################################################################
###############################################################################

#%% Set PATH to save
output_dir = project_dir / "Resultados" / "Joined_datasets" / "Integration_methods_lab" / "scIB_Merged"

###############################################################################
###############################################################################

#%% Define paths and embeddings

batch_key = "Product_norm"
label_key = "manual_celltype_annotation_high"

methods = {
    "Merged": {
        "path": _input_path(project_dir / "Resultados" / "Joined_datasets" / "Integration_methods_lab" / "Merged_WO_integration" / "Sin_GT", "Seurat_merged.h5ad"),
        "source_embedding": "X_pca_wo_integ",
        "target_embedding": "X_Merged",
    },
    "scVI": {
        "path": _input_path(project_dir / "Resultados" / "Joined_datasets" / "Integration_methods_lab" / "scVI" / "Sin_GT_With_Python", "Suplementaria_S2_scVI_state1.h5ad"),
        "source_embedding": "X_scVI",
        "target_embedding": "X_scVI",
    },
    "Harmony": {
        "path": _input_path(project_dir / "Resultados" / "Joined_datasets" / "Integration_methods_lab" / "Harmony" / "Sin_GT", "Seurat_harmony.h5ad"),
        "source_embedding": "X_harmony",
        "target_embedding": "X_Harmony",
    },
    "LIGER": {
        "path": _input_path(project_dir / "Resultados" / "Joined_datasets" / "Integration_methods_lab" / "LIGER" / "Sin_GT", "Seurat_liger.h5ad"),
        "source_embedding": "X_iNMF",
        "target_embedding": "X_LIGER",
    },
    "STACAS": {
        "path": _input_path(project_dir / "Resultados" / "Joined_datasets" / "Integration_methods_lab" / "STACAS" / "Sin_GT", "Seurat_STACAS_integ2.h5ad"),
        "source_embedding": "X_pca",
        "target_embedding": "X_STACAS",
    },
    "RPCA": {
        "path": _input_path(project_dir / "Resultados" / "Joined_datasets" / "Integration_methods_lab" / "Seurat_RPCA" / "Sin_GT", "Seurat_RPCA_integ.h5ad"),
        "source_embedding": "X_RPCA_pca",
        "target_embedding": "X_RPCA",
    },
    "fastMNN": {
        "path": _input_path(project_dir / "Resultados" / "Joined_datasets" / "Integration_methods_lab" / "fastMNN" / "Sin_GT", "Seurat_fastMNN2.h5ad"),
        "source_embedding": "X_mnn",
        "target_embedding": "X_fastMNN",
    },
}

###############################################################################
###############################################################################

#%% Load reference AnnData

# Use merged object as reference because it contains the common obs metadata
adata_scIB_all = sc.read_h5ad(methods["Merged"]["path"])
adata_scIB_all.obs_names_make_unique()

print("[bold green]Reference object[/bold green]")
print(adata_scIB_all)
print(f"Reference cells: {adata_scIB_all.n_obs}")
print(f"Reference batches/products: {adata_scIB_all.obs[batch_key].nunique()}")
print(f"Reference labels: {adata_scIB_all.obs[label_key].nunique()}")

###############################################################################
###############################################################################

#%% Add all method-specific embeddings to the reference object

embedding_obsm_keys = []

for method, info in methods.items():

    print(f"\n[bold green]Processing {method}[/bold green]")
    adata_method = sc.read_h5ad(info["path"])
    adata_method.obs_names_make_unique()
    source_embedding = info["source_embedding"]
    target_embedding = info["target_embedding"]
    print(f"Available embeddings: {list(adata_method.obsm.keys())}")
    if source_embedding not in adata_method.obsm.keys():
        raise ValueError(
            f"{source_embedding} not found for {method}. "
            f"Available embeddings are: {list(adata_method.obsm.keys())}"
        )
    # Check same cells
    same_cells = set(adata_method.obs_names) == set(adata_scIB_all.obs_names)
    same_number_cells = adata_method.n_obs == adata_scIB_all.n_obs
    print(f"Same number of cells: {same_number_cells}")
    print(f"Same cells: {same_cells}")
    if not same_cells:
        missing_in_method = set(adata_scIB_all.obs_names) - set(adata_method.obs_names)
        extra_in_method = set(adata_method.obs_names) - set(adata_scIB_all.obs_names)
        raise ValueError(
            f"{method} does not contain the same cells as the reference object.\n"
            f"Missing in {method}: {len(missing_in_method)} cells\n"
            f"Extra in {method}: {len(extra_in_method)} cells"
        )
    # Reorder cells to match reference order
    adata_method = adata_method[adata_scIB_all.obs_names].copy()
    # Optional checks for metadata consistency
    same_batch = (
        adata_method.obs[batch_key].astype(str).values
        == adata_scIB_all.obs[batch_key].astype(str).values
    ).all()
    same_label = (
        adata_method.obs[label_key].astype(str).values
        == adata_scIB_all.obs[label_key].astype(str).values
    ).all()
    print(f"Same {batch_key}: {same_batch}")
    print(f"Same {label_key}: {same_label}")
    if not same_batch:
        raise ValueError(f"{batch_key} differs between reference and {method}")
    if not same_label:
        raise ValueError(f"{label_key} differs between reference and {method}")
    # Copy method-specific embedding
    adata_scIB_all.obsm[target_embedding] = adata_method.obsm[source_embedding].copy()
    print(f"Copied {source_embedding} as {target_embedding}")
    print(f"Embedding shape: {adata_scIB_all.obsm[target_embedding].shape}")
    embedding_obsm_keys.append(target_embedding)

###############################################################################
###############################################################################

#%% Save combined AnnData with all embeddings

adata_scIB_all.write_h5ad(_output_path(output_dir, "adata_scIB_all_methods_embeddings.h5ad"))

print("\n[bold green]Embeddings to benchmark[/bold green]")
print(embedding_obsm_keys)

###############################################################################
###############################################################################

#%% Run one single Benchmarker

bm_all = Benchmarker(
    adata_scIB_all,
    batch_key=batch_key,
    label_key=label_key,
    bio_conservation_metrics=BioConservation(),
    batch_correction_metrics=BatchCorrection(),
    embedding_obsm_keys=embedding_obsm_keys,
    n_jobs=32,
)

bm_all.benchmark()

###############################################################################
###############################################################################

#%% Get and save results

df_all_unscaled = bm_all.get_results(min_max_scale=False)
df_all_scaled = bm_all.get_results(min_max_scale=True)

print("\n[bold green]Unscaled results[/bold green]")
print(df_all_unscaled)

print("\n[bold green]Scaled results[/bold green]")
print(df_all_scaled)

df_all_unscaled.to_csv(_output_path(output_dir, "resultados_benchmark_scIB_all_methods_unscaled.csv"), index=True)
df_all_scaled.to_csv(_output_path(output_dir, "resultados_benchmark_scIB_all_methods_scaled.csv"), index=True)

###############################################################################
###############################################################################

#%% Plot results table

table_unscaled = bm_all.plot_results_table(min_max_scale=False)
table_scaled = bm_all.plot_results_table(min_max_scale=True)

table_unscaled.figure.savefig(
    _output_path(output_dir, "scIB_results_table_unscaled.pdf"),
    bbox_inches="tight",
    dpi=300
)

table_scaled.figure.savefig(
    _output_path(output_dir, "scIB_results_table_scaled.pdf"),
    bbox_inches="tight",
    dpi=300
)
###############################################################################
###############################################################################

#%% End of script
