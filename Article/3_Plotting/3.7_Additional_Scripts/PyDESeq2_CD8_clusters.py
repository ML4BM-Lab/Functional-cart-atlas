###############################################################################
###############################################################################

# Program: PyDESeq2_CD8_clusters.py
# Author: Sergio Cámara Peña
# Date: 20/05/2026
# Version: FINAL
# More info: Compare against "Dreamlet_CD8_clusters.R"
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
    if not path.exists():
        raise FileNotFoundError(f"Required input path does not exist: {path}")
    return path


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
if not path.is_dir():
    raise FileNotFoundError(f"Required input directory does not exist: {path}")
file = _input_path(path, "Python_scVI_adata_big_V4_state4.h5ad")

adata = sc.read_h5ad(file)

adata.layers["counts"] = adata.X.copy()

print(adata)
print("Number of cells before filtering:", adata.n_obs)
print("Number of genes:", adata.n_vars)

# Check counts layer
print("Max value in counts layer:", adata.layers["counts"].max())

# Load significant signatures from dreamlet
base_dir = project_dir / "Resultados" / "CD8_clusters"
if not base_dir.is_dir():
    raise FileNotFoundError(f"Required input directory does not exist: {base_dir}")

with open(_input_path(base_dir, "genes_dreamlet_1.txt"), "r") as f:
    genes_dreamlet_1 = [line.strip() for line in f if line.strip()]

with open(_input_path(base_dir, "genes_dreamlet_2.txt"), "r") as f:
    genes_dreamlet_2 = [line.strip() for line in f if line.strip()]

with open(_input_path(base_dir, "genes_dreamlet_3.txt"), "r") as f:
    genes_dreamlet_3 = [line.strip() for line in f if line.strip()]

# %% Filter object
filtered_adata = adata[adata.obs["Antigen"] == "Blood"].copy()
print("After Antigen == Blood:", filtered_adata.n_obs)

filtered_adata = filtered_adata[
    filtered_adata.obs["manual_celltype_annotation_high"] == "CD8 cytotoxic"
].copy()
print("After CD8 cytotoxic:", filtered_adata.n_obs)

filtered_adata = filtered_adata[
    ~filtered_adata.obs["Max_Response"].isna()
].copy()
print("After removing NA Max_Response:", filtered_adata.n_obs)

filtered_adata = filtered_adata[
    filtered_adata.obs["Max_Response"] == "CR"
].copy()
print("After Max_Response == CR:", filtered_adata.n_obs)

filtered_adata = filtered_adata[
    filtered_adata.obs["STATUS"] == "DISEASE"
].copy()
print("After STATUS == DISEASE:", filtered_adata.n_obs)

filtered_adata = filtered_adata[
    ~(
        (filtered_adata.obs["Time_Point_Ranges"] == "Infusion_Product") &
        (filtered_adata.obs["Stimulated"] == "YES")
    )
].copy()
print("After removing stimulated Infusion_Product:", filtered_adata.n_obs)

filtered_adata = filtered_adata[
    ~(
        (filtered_adata.obs["Stimulated"] == "YES") &
        (filtered_adata.obs["Stimulation_Location"] == "In_vitro")
    )
].copy()
print("After removing in vitro stimulated cells:", filtered_adata.n_obs)

# %% Set categorical variables
filtered_adata.obs["Time_Point_Ranges"] = pd.Categorical(
    filtered_adata.obs["Time_Point_Ranges"],
    categories=["Infusion_Product", "<2_weeks", "2_weeks-3_months", ">3_months"],
    ordered=True
)

filtered_adata.obs["Age_Range"] = pd.Categorical(
    filtered_adata.obs["Age_Range"],
    categories=["<20", "20-40", "40-60", ">60"],
    ordered=True
)

print(filtered_adata.obs["ScFv"].value_counts())

# %% Create pseudobulk without pseudo replicates
pbs = []
for sample in filtered_adata.obs.Product_norm.unique():
    samp_filtered_adata = filtered_adata[filtered_adata.obs['Product_norm'] == sample]
    
    samp_filtered_adata.X = samp_filtered_adata.layers['counts'] #make sure to use raw data
    
    rep_adata = sc.AnnData(X = samp_filtered_adata.X.sum(axis = 0),
                           var = samp_filtered_adata.var[[]])
    
    rep_adata.obs_names = [sample]
    rep_adata.obs['Time_Point_Ranges'] = samp_filtered_adata.obs['Time_Point_Ranges'].iloc[0]
    
    pbs.append(rep_adata)

    print(sample, " DONE")

pb = sc.concat(pbs)
print(pb.obs)

# %% Create PyDESeq2 object
counts = pd.DataFrame(pb.X, columns = pb.var_names) #need to do this to pass var names

dds = DeseqDataSet(
    counts = counts,
    metadata=pb.obs,
    design_factors="Time_Point_Ranges")

sc.pp.filter_genes(dds, min_cells = 1)

dds

# %% Execute deseq2
dds.deseq2()

# %% Comparison 1: <2_weeks vs Infusion_Product
stat_res_1 = DeseqStats(dds, contrast=("Time-Point-Ranges", "<2-weeks", "Infusion-Product"))
    
stat_res_1.summary()


## Comparison 1: Create table with the genes that dreamlet gave to compare results
de_1  = stat_res_1.results_df

genes_interes_1 = [
    "SMIM3",
    "MSMO1",
    "LINC01943",
    "FAM174A",
    "MIB1",
    "HMGCS1",
    "DHCR7",
    "PLAAT4",
    "ARAP2",
    "OSM",
]

de_genes_1 = de_1.reindex(genes_interes_1)

print(de_genes_1)

# Save
output_dir = project_dir / "Resultados" / "Dreamlet_Vs_PyDESeq2"
output_file = _output_path(output_dir, "PyDESeq2_de_genes_CD8_Comparison_1.csv")
de_genes_1.to_csv(output_file, index=True)

################################################################################################################################

##### Compare significant genes (padj < 0.05) between both tools #####
# Significant genes from PyDESeq2
genes_pydeseq2_1 = de_1[de_1.padj < 0.05].index.tolist()

# Convert both lists to sets
set_dreamlet_1 = set(genes_dreamlet_1)
set_pydeseq2_1 = set(genes_pydeseq2_1)

# Genes shared between Dreamlet and PyDESeq2
genes_common_1 = sorted(set_dreamlet_1 & set_pydeseq2_1)

# Genes only in Dreamlet
genes_only_dreamlet_1 = sorted(set_dreamlet_1 - set_pydeseq2_1)

# Genes only in PyDESeq2
genes_only_pydeseq2_1 = sorted(set_pydeseq2_1 - set_dreamlet_1)

# Directional overlap percentages
pct_dreamlet_in_pydeseq2_1 = (
    len(genes_common_1) / len(set_dreamlet_1) * 100
    if len(set_dreamlet_1) > 0 else 0
)

pct_pydeseq2_in_dreamlet_1 = (
    len(genes_common_1) / len(set_pydeseq2_1) * 100
    if len(set_pydeseq2_1) > 0 else 0
)

print(f"Genes Dreamlet padj < 0.05: {len(set_dreamlet_1)}")
print(f"Genes PyDESeq2 padj < 0.05: {len(set_pydeseq2_1)}")
print(f"Common genes: {len(genes_common_1)}")
print(f"Only in Dreamlet: {len(genes_only_dreamlet_1)}")
print(f"Only in PyDESeq2: {len(genes_only_pydeseq2_1)}")
print(f"% of genes Dreamlet present in PyDESeq2: {pct_dreamlet_in_pydeseq2_1:.2f}%")
print(f"% of genes PyDESeq2 present in Dreamlet: {pct_pydeseq2_in_dreamlet_1:.2f}%")

# %% Save comparison summary
comparison_summary_1 = {
    "Dreamlet significant genes": len(set_dreamlet_1),
    "PyDESeq2 significant genes": len(set_pydeseq2_1),
    "Common significant genes": len(genes_common_1),
    "Dreamlet-only genes": len(genes_only_dreamlet_1),
    "PyDESeq2-only genes": len(genes_only_pydeseq2_1),
    "Dreamlet genes found in PyDESeq2 (%)": pct_dreamlet_in_pydeseq2_1,
    "PyDESeq2 genes found in Dreamlet (%)": pct_pydeseq2_in_dreamlet_1,
}

summary_df_1 = pd.DataFrame(
    comparison_summary_1.items(),
    columns=["Metric", "Value"]
)

summary_file_1 = _output_path(
    output_dir,
    "PyDESeq2_Dreamlet_CD8_comparison_1_summary.csv"
)

summary_df_1.to_csv(summary_file_1, index=False)

# %% Comparison 2: 2_weeks-3_months vs <2_weeks
stat_res_2 = DeseqStats(dds, contrast=("Time-Point-Ranges", "2-weeks-3-months", "<2-weeks"))
    
stat_res_2.summary()


## Comparison 2: Create table with the genes that dreamlet gave to compare results
de_2  = stat_res_2.results_df

genes_interes_2 = [
    "MRNIP",
    "TOP1MT",
    "SELENOM",
    "LDLRAP1",
    "ARRDC3",
    "PLEKHB1",
    "CLIP4",
    "COTL1",
    "ZNF224",
    "BBS2",
]

de_genes_2 = de_2.reindex(genes_interes_2)

print(de_genes_2)

# Save
output_dir = project_dir / "Resultados" / "Dreamlet_Vs_PyDESeq2"
output_file = _output_path(output_dir, "PyDESeq2_de_genes_CD8_Comparison_2.csv")
de_genes_2.to_csv(output_file, index=True)

################################################################################################################################

##### Compare significant genes (padj < 0.05) between both tools #####
# Significant genes from PyDESeq2
genes_pydeseq2_2 = de_2[de_2.padj < 0.05].index.tolist()

# Convert both lists to sets
set_dreamlet_2 = set(genes_dreamlet_2)
set_pydeseq2_2 = set(genes_pydeseq2_2)

# Genes shared between Dreamlet and PyDESeq2
genes_common_2 = sorted(set_dreamlet_2 & set_pydeseq2_2)

# Genes only in Dreamlet
genes_only_dreamlet_2 = sorted(set_dreamlet_2 - set_pydeseq2_2)

# Genes only in PyDESeq2
genes_only_pydeseq2_2 = sorted(set_pydeseq2_2 - set_dreamlet_2)

# Directional overlap percentages
pct_dreamlet_in_pydeseq2_2 = (
    len(genes_common_2) / len(set_dreamlet_2) * 100
    if len(set_dreamlet_2) > 0 else 0
)

pct_pydeseq2_in_dreamlet_2 = (
    len(genes_common_2) / len(set_pydeseq2_2) * 100
    if len(set_pydeseq2_2) > 0 else 0
)

print(f"Genes Dreamlet padj < 0.05: {len(set_dreamlet_2)}")
print(f"Genes PyDESeq2 padj < 0.05: {len(set_pydeseq2_2)}")
print(f"Common genes: {len(genes_common_2)}")
print(f"Only in Dreamlet: {len(genes_only_dreamlet_2)}")
print(f"Only in PyDESeq2: {len(genes_only_pydeseq2_2)}")
print(f"% of genes Dreamlet present in PyDESeq2: {pct_dreamlet_in_pydeseq2_2:.2f}%")
print(f"% of genes PyDESeq2 present in Dreamlet: {pct_pydeseq2_in_dreamlet_2:.2f}%")

# %% Save comparison summary
comparison_summary_2 = {
    "Dreamlet significant genes": len(set_dreamlet_2),
    "PyDESeq2 significant genes": len(set_pydeseq2_2),
    "Common significant genes": len(genes_common_2),
    "Dreamlet-only genes": len(genes_only_dreamlet_2),
    "PyDESeq2-only genes": len(genes_only_pydeseq2_2),
    "Dreamlet genes found in PyDESeq2 (%)": pct_dreamlet_in_pydeseq2_2,
    "PyDESeq2 genes found in Dreamlet (%)": pct_pydeseq2_in_dreamlet_2,
}

summary_df_2 = pd.DataFrame(
    comparison_summary_2.items(),
    columns=["Metric", "Value"]
)

summary_file_2 = _output_path(
    output_dir,
    "PyDESeq2_Dreamlet_CD8_comparison_2_summary.csv"
)

summary_df_2.to_csv(summary_file_2, index=False)

# %% Comparison 3: >3_months vs 2_weeks-3_months
stat_res_3 = DeseqStats(dds, contrast=("Time-Point-Ranges", ">3-months", "2-weeks-3-months"))
    
stat_res_3.summary()

## Comparison 3: Create table with the genes that dreamlet gave to compare results
de_3  = stat_res_3.results_df

genes_interes_3 = [
    "MEI1",
    "PPP3R1",
    "IFITM2",
    "EPS8L2",
    "FUT8",
    "COQ7",
    "DHRS4",
    "WAKMAR2",
    "IER2",
    "ANTKMT",
]

de_genes_3 = de_3.reindex(genes_interes_3)

print(de_genes_3)

# Save
output_dir = project_dir / "Resultados" / "Dreamlet_Vs_PyDESeq2"
output_file = _output_path(output_dir, "PyDESeq2_de_genes_CD8_Comparison_3.csv")
de_genes_3.to_csv(output_file, index=True)

################################################################################################################################

##### Compare significant genes (padj < 0.05) between both tools #####
# Significant genes from PyDESeq2
genes_pydeseq2_3 = de_3[de_3.padj < 0.05].index.tolist()

# Convert both lists to sets
set_dreamlet_3 = set(genes_dreamlet_3)
set_pydeseq2_3 = set(genes_pydeseq2_3)

# Genes shared between Dreamlet and PyDESeq2
genes_common_3 = sorted(set_dreamlet_3 & set_pydeseq2_3)

# Genes only in Dreamlet
genes_only_dreamlet_3 = sorted(set_dreamlet_3 - set_pydeseq2_3)

# Genes only in PyDESeq2
genes_only_pydeseq2_3 = sorted(set_pydeseq2_3 - set_dreamlet_3)

# Directional overlap percentages
pct_dreamlet_in_pydeseq2_3 = (
    len(genes_common_3) / len(set_dreamlet_3) * 100
    if len(set_dreamlet_3) > 0 else 0
)

pct_pydeseq2_in_dreamlet_3 = (
    len(genes_common_3) / len(set_pydeseq2_3) * 100
    if len(set_pydeseq2_3) > 0 else 0
)

print(f"Genes Dreamlet padj < 0.05: {len(set_dreamlet_3)}")
print(f"Genes PyDESeq2 padj < 0.05: {len(set_pydeseq2_3)}")
print(f"Common genes: {len(genes_common_3)}")
print(f"Only in Dreamlet: {len(genes_only_dreamlet_3)}")
print(f"Only in PyDESeq2: {len(genes_only_pydeseq2_3)}")
print(f"% of genes Dreamlet present in PyDESeq2: {pct_dreamlet_in_pydeseq2_3:.2f}%")
print(f"% of genes PyDESeq2 present in Dreamlet: {pct_pydeseq2_in_dreamlet_3:.2f}%")

# %% Save comparison summary
comparison_summary_3 = {
    "Dreamlet significant genes": len(set_dreamlet_3),
    "PyDESeq2 significant genes": len(set_pydeseq2_3),
    "Common significant genes": len(genes_common_3),
    "Dreamlet-only genes": len(genes_only_dreamlet_3),
    "PyDESeq2-only genes": len(genes_only_pydeseq2_3),
    "Dreamlet genes found in PyDESeq2 (%)": pct_dreamlet_in_pydeseq2_3,
    "PyDESeq2 genes found in Dreamlet (%)": pct_pydeseq2_in_dreamlet_3,
}

summary_df_3 = pd.DataFrame(
    comparison_summary_3.items(),
    columns=["Metric", "Value"]
)

summary_file_3 = _output_path(
    output_dir,
    "PyDESeq2_Dreamlet_CD8_comparison_3_summary.csv"
)

summary_df_3.to_csv(summary_file_3, index=False)

# %% End of script
