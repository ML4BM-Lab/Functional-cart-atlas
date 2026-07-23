###############################################################################
###############################################################################

# Program: PyDESeq2_IACs_clusters.py
# Author: Sergio Cámara Peña
# Date: 20/05/2026
# Version: FINAL
# More info: Compare against "Dreamlet_IACs_clusters.R"
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


def _output_path(directory, filename):
    directory.mkdir(parents=True, exist_ok=True)
    return directory / filename


# %% Import libraries
import os
import random
import numpy as np
import pandas as pd
import scanpy as sc
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

# %% Set random seed
random.seed(2504)

# %% Read files
path = project_dir / "Input"
file = _input_path(path, "Python_scVI_adata_big_V4_state4.h5ad")

adata = sc.read_h5ad(file)

adata.layers["counts"] = adata.X.copy()

print(adata)
print("Number of cells before filtering:", adata.n_obs)
print("Number of genes:", adata.n_vars)

# Check counts layer
print("Max value in counts layer:", adata.layers["counts"].max())

# %% Filter object
filtered_adata = adata[adata.obs["Antigen"] == "Blood"].copy()
print("After Antigen == Blood:", filtered_adata.n_obs)

filtered_adata = filtered_adata[
    filtered_adata.obs["manual_celltype_annotation_high"] == "Monocyte-like T cells"
].copy()
print("After Monocyte-like T cells:", filtered_adata.n_obs)

filtered_adata = filtered_adata[
    filtered_adata.obs["STATUS"] == "DISEASE"
].copy()
print("After STATUS == DISEASE:", filtered_adata.n_obs)

# %% Set categorical variables
filtered_adata.obs["Max_Response"] = pd.Categorical(
    filtered_adata.obs["Max_Response"],
    categories=["NR", "PR", "CR"],
    ordered=True
)

filtered_adata.obs["ICANS_Grade_Range"] = pd.Categorical(
    filtered_adata.obs["ICANS_Grade_Range"],
    categories=["3-4", "1-2", "0"],
    ordered=True
)

filtered_adata.obs["Time_Point_Ranges"] = pd.Categorical(
    filtered_adata.obs["Time_Point_Ranges"],
    categories=["Infusion_Product", "<2_weeks", "2_weeks-3_months"],
    ordered=True
)

filtered_adata.obs["Age_Range"] = pd.Categorical(
    filtered_adata.obs["Age_Range"],
    categories=["<20", "20-40", "40-60", ">60"],
    ordered=True
)

# %% Create IACs_differenciation column
# Equivalent to:
# filtered_sce$IACs_differenciation <- 'Post_Infusion'
# filtered_sce$IACs_differenciation[filtered_sce$Time_Point_Ranges == 'Infusion_Product'] <- 'Infusion_Product'

filtered_adata.obs["IACs_differenciation"] = "Post_Infusion"

filtered_adata.obs.loc[
    filtered_adata.obs["Time_Point_Ranges"] == "Infusion_Product",
    "IACs_differenciation"
] = "Infusion_Product"

filtered_adata.obs["IACs_differenciation"] = pd.Categorical(
    filtered_adata.obs["IACs_differenciation"],
    categories=["Infusion_Product", "Post_Infusion"],
    ordered=True
)

# %% Create pseudobulk without pseudo replicates
pbs = []
for sample in filtered_adata.obs.Product_norm.unique():
    samp_filtered_adata = filtered_adata[filtered_adata.obs['Product_norm'] == sample]
    
    samp_filtered_adata.X = samp_filtered_adata.layers['counts'] #make sure to use raw data
    
    rep_adata = sc.AnnData(X = samp_filtered_adata.X.sum(axis = 0),
                           var = samp_filtered_adata.var[[]])
    
    rep_adata.obs_names = [sample]
    rep_adata.obs['IACs_differenciation'] = samp_filtered_adata.obs['IACs_differenciation'].iloc[0]
    
    pbs.append(rep_adata)

    print(sample, " DONE")

pb = sc.concat(pbs)
print(pb.obs)

# %% Create PyDESeq2 object
counts = pd.DataFrame(pb.X, columns = pb.var_names) #need to do this to pass var names

dds = DeseqDataSet(
    counts = counts,
    metadata=pb.obs,
    design_factors="IACs_differenciation")

sc.pp.filter_genes(dds, min_cells = 1)

dds

# %% Execute deseq2
dds.deseq2()

# %% Create the comparison
stat_res = DeseqStats(dds, contrast=('IACs-differenciation', 'Post-Infusion', 'Infusion-Product'))
    
stat_res.summary()


# %% Create table with the genes that dreamlet gave to compare results
de  = stat_res.results_df

genes_interes = [
    "BLOC1S1",
    "RIPOR2",
    "LTA4H",
    "AIF1",
    "APBB1IP",
    "NDUFB8",
    "LST1",
    "TXN",
    "C19orf38",
    "ZFP36",
]

de_genes = de.reindex(genes_interes)

print(de_genes)

# Save
output_dir = project_dir / "Resultados" / "Dreamlet_Vs_PyDESeq2"
output_file = _output_path(output_dir, "PyDESeq2_de_genes_IACs_clusters.csv")
de_genes.to_csv(output_file, index=True)

# %% Compare significant genes (padj < 0.05) between both tools
base_dir = project_dir / "Para_Cluster_Profiler"
if not base_dir.is_dir():
    raise FileNotFoundError(f"Required input directory does not exist: {base_dir}")

with open(_input_path(base_dir, "genes_dreamlet.txt"), "r") as f:
    genes_dreamlet = [line.strip() for line in f if line.strip()]

# Significant genes from PyDESeq2
genes_pydeseq2 = de[de.padj < 0.05].index.tolist()

# Convert both lists to sets
set_dreamlet = set(genes_dreamlet)
set_pydeseq2 = set(genes_pydeseq2)

# Genes shared between Dreamlet and PyDESeq2
genes_common = sorted(set_dreamlet & set_pydeseq2)

# Genes only in Dreamlet
genes_only_dreamlet = sorted(set_dreamlet - set_pydeseq2)

# Genes only in PyDESeq2
genes_only_pydeseq2 = sorted(set_pydeseq2 - set_dreamlet)

# Directional overlap percentages
pct_dreamlet_in_pydeseq2 = (
    len(genes_common) / len(set_dreamlet) * 100
    if len(set_dreamlet) > 0 else 0
)

pct_pydeseq2_in_dreamlet = (
    len(genes_common) / len(set_pydeseq2) * 100
    if len(set_pydeseq2) > 0 else 0
)

print(f"Genes Dreamlet padj < 0.05: {len(set_dreamlet)}")
print(f"Genes PyDESeq2 padj < 0.05: {len(set_pydeseq2)}")
print(f"Genes comunes: {len(genes_common)}")
print(f"Solo en Dreamlet: {len(genes_only_dreamlet)}")
print(f"Solo en PyDESeq2: {len(genes_only_pydeseq2)}")
print(f"% de genes Dreamlet presentes en PyDESeq2: {pct_dreamlet_in_pydeseq2:.2f}%")
print(f"% de genes PyDESeq2 presentes en Dreamlet: {pct_pydeseq2_in_dreamlet:.2f}%")

# %% Save comparison summary
comparison_summary = {
    "Dreamlet significant genes": len(set_dreamlet),
    "PyDESeq2 significant genes": len(set_pydeseq2),
    "Common significant genes": len(genes_common),
    "Dreamlet-only genes": len(genes_only_dreamlet),
    "PyDESeq2-only genes": len(genes_only_pydeseq2),
    "Dreamlet genes found in PyDESeq2 (%)": pct_dreamlet_in_pydeseq2,
    "PyDESeq2 genes found in Dreamlet (%)": pct_pydeseq2_in_dreamlet,
}

summary_df = pd.DataFrame(
    comparison_summary.items(),
    columns=["Metric", "Value"]
)

summary_file = _output_path(
    output_dir,
    "PyDESeq2_Dreamlet_IACs_clusters_comparison_summary.csv"
)

summary_df.to_csv(summary_file, index=False)

# %% End of script
