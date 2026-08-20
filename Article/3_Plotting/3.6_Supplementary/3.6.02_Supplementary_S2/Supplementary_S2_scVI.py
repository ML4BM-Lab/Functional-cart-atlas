###############################################################################
###############################################################################

# Program: Supplementary_S2_scVI.py
# Author: Sergio Cámara Peña
# Date: 09/06/2025
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
    project_dir / "Resultados" / "Joined_datasets" / "Integration_methods_lab" / "scVI" / "Sin_GT_With_Python" / "Seurat_merged_RAW_for_Py.h5ad",
}

def _require_path(path):
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

##### Use of scVI in its Python version #####
## More info: https://ccbskillssem.github.io/assets/scvi_notebook.html and https://docs.scvi-tools.org/en/stable/tutorials/notebooks/harmonization.html

#%% Load all the needed libraries
import scanpy as sc
import scvi
from scvi.model.utils import mde
import scanpy.external as sce
import os
import random
import pandas as pd
import seaborn as sns

#%% Set a random seed
random.seed(2504)

# %% Define a 12 color palette
palette = [
"#F8766D",
"#00BA38",
"#FF9900",
"#619CFF",
"#F564E3",
"#B79F00"
]

# %% Switches
Primera_vez = True

#%% Read the data
adata = sc.read(_require_path(project_dir / 'Resultados' / 'Joined_datasets' / 'Integration_methods_lab' / 'scVI' / 'Sin_GT_With_Python' / 'Seurat_merged_RAW_for_Py.h5ad')) # Obtained in Supplementary_S2_Methods.R script
adata.raw = adata  # keep full dimension safe

# %% Obtain the cellular annotation
data_dir = project_dir / 'Resultados' / 'Joined_datasets' / 'Raw_Atlas'
adata_orig = sc.read_h5ad(_input_path(data_dir, "Python_scVI_adata_big_V4_state4.h5ad"))

adata_orig.obs.index = adata_orig.obs.index.str.replace(r"_\d+$", "", regex=True)

common_cells = adata.obs_names.intersection(adata_orig.obs_names)
print(f"Common cells: {len(common_cells)} of {adata.n_obs}")

adata_subset = adata[common_cells].copy()

adata_subset.obs["manual_celltype_annotation_high"] = adata_orig.obs.loc[common_cells, "manual_celltype_annotation_high"]

results_dir = project_dir / 'Resultados' / 'Joined_datasets' / 'Integration_methods_lab'
df_annotation = adata_subset.obs[["manual_celltype_annotation_high"]].copy()
df_annotation["cell_name"] = df_annotation.index
df_annotation = df_annotation.reset_index(drop=True)
df_annotation.to_csv(_output_path(results_dir, "cell_annotation_manual.csv"), index=False)

adata_subset = adata[common_cells].copy()
adata_subset.obs["manual_celltype_annotation_high"] = adata_orig.obs.loc[common_cells, "manual_celltype_annotation_high"]
del adata

adata = adata_subset.copy()

#%% Read the metadata, process and adapt it
df = pd.read_csv(_require_path(project_dir / 'Datasets_Metadata' / 'scRNAseq_metadata_CARTs_v10.csv'), skipfooter=6, sep=";")
Data_rows = adata.obs.Product.unique()
filtered_df = df[df["Sample name"].isin(Data_rows)]
filtered_df = filtered_df[["Sample name", "Age", "Age_Range", "Sex", "CAR_Construct", "CAR_Gen", "ScFv", "Costim_Domain_1", "Costim_Domain_2", "STATUS", "Stimulated"]]

#%% Normalize the data
adata.layers["counts"] = adata.X.copy() # preserve counts
sc.pp.normalize_total(adata, target_sum=1e4) # scale each cell to a common library size
sc.pp.log1p(adata) # log(expression + 1)

#%% Identify and only keep the 2000 most variable genes
sc.pp.highly_variable_genes(
    adata,
    flavor="seurat_v3",
    n_top_genes=2000,
    batch_key="Product_norm",
    subset=True,
    layer="counts"
)

#%% Standard Workflow
sc.pp.scale(adata) # Normalize the columns (genes)
sc.tl.pca(adata)

adata.obsm["X_pca"]

sc.pp.neighbors(adata) # Compute nearest neighbors
sc.tl.umap(adata)

sc.pl.umap(adata, color="Product")

sce.pp.bbknn(adata, batch_key="Product")
sc.tl.umap(adata)

#%% scVI setup
if Primera_vez:
    scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="Product")
    model = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
    scvi.model.SCVI.view_anndata_setup(model)

#%% Train the model
if Primera_vez:
    model.train()

#%% Save/Load trained model
model_dir = project_dir / 'Resultados' / 'Joined_datasets' / 'Integration_methods_lab' / 'scVI' / 'Sin_GT_With_Python'

if Primera_vez:
    model.save(_output_path(model_dir, "Healty_donors_scvi_v1"), overwrite=True, save_anndata=True)
else:
    model = scvi.model.SCVI.load(_input_path(model_dir, 'Healty_donors_scvi_v1'), adata)

#%% Obtain latent representation from scVI to evaluate it
adata.obsm["X_scVI"] = model.get_latent_representation()

# %% Neighbours, leiden (to cluster cells) and UMAP calculation
sc.pp.neighbors(adata, use_rep="X_scVI")
sc.tl.leiden(adata, key_added="leiden_scvi", resolution=1.2)
sc.tl.umap(adata, min_dist=0.4)

# %% Save or read .h5ad
Editar_Figura = False
model_dir = project_dir / 'Resultados' / 'Joined_datasets' / 'Integration_methods_lab' / 'scVI' / 'Sin_GT_With_Python'
if Editar_Figura:
    adata = sc.read_h5ad(_input_path(model_dir, "Suplementaria_S2_scVI_state1.h5ad"))
    print("Correctly loaded")
else:
    adata.write(_output_path(model_dir, "Suplementaria_S2_scVI_state1.h5ad"))
    print("Correctly saved")

#%% UMAP representation and save
figures_dir = project_dir / 'Resultados_Figuras' / 'Suplementarias'
figures_dir.mkdir(parents=True, exist_ok=True)
sc.settings.figdir = figures_dir

sns.palplot(palette)

sc.pl.umap(
    adata,
    color=["orig.ident"],
    frameon=False,
    size=4,
    palette=palette,
    save="Suplementaria_S2_2_scVI.pdf"
)

# %% End of script
