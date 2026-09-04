#!/usr/bin/env python3
"""Minimal runtime check for the Margaret Python image."""

import argparse
import importlib.metadata
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import anndata
import bbknn
import celltypist
import gseapy
import invoke
import milopy
import milopy.core as milo
import numpy as np
import pandas as pd
import plotnine
import rich
import scanpy as sc
import scarches
import scgraph
import scib_metrics
import scipy
import scProportionTest
import seaborn
import scvi
import torch
from rpy2.robjects import r
from scib_metrics.benchmark import BatchCorrection, Benchmarker, BioConservation
from skmisc.loess import loess


def package_version(distribution: str) -> str:
    return importlib.metadata.version(distribution)


parser = argparse.ArgumentParser()
parser.add_argument("demo", type=Path, help="Path to Atlas_DEMO.h5ad")
args = parser.parse_args()

if sys.version_info[:3] != (3, 8, 10):
    raise RuntimeError(f"Expected Python 3.8.10, found {sys.version.split()[0]}")
if not args.demo.is_file():
    raise FileNotFoundError(args.demo)

demo = anndata.read_h5ad(args.demo, backed="r")
try:
    if demo.n_obs == 0 or demo.n_vars == 0:
        raise RuntimeError("Atlas_DEMO.h5ad is empty")
    demo_shape = demo.shape
finally:
    demo.file.close()

mini = anndata.AnnData(
    np.array(
        [
            [1.0, 0.0, 3.0, 1.0],
            [0.0, 2.0, 1.0, 1.0],
            [4.0, 1.0, 0.0, 2.0],
            [2.0, 3.0, 1.0, 0.0],
            [1.0, 1.0, 2.0, 3.0],
        ],
        dtype=np.float32,
    )
)
sc.pp.normalize_total(mini, target_sum=1_000)
sc.pp.log1p(mini)
sc.tl.pca(mini, n_comps=2)
if mini.obsm["X_pca"].shape != (5, 2):
    raise RuntimeError("Unexpected Scanpy PCA output")

# Scanpy's Seurat v3 highly-variable-gene method imports scikit-misc lazily.
# Exercise the actual code path rather than checking only that Scanpy imports.
rng = np.random.default_rng(2504)
hvg = anndata.AnnData(rng.poisson(5, size=(80, 100)).astype(np.float32))
hvg.layers["counts"] = hvg.X.copy()
hvg.obs["batch"] = np.repeat(["batch_a", "batch_b"], 40)
sc.pp.highly_variable_genes(
    hvg,
    flavor="seurat_v3",
    n_top_genes=20,
    layer="counts",
    batch_key="batch",
    span=0.6,
)
if int(hvg.var["highly_variable"].sum()) != 20:
    raise RuntimeError("Unexpected Seurat v3 highly-variable-gene result")

product = torch.tensor([[1.0, 2.0]]) @ torch.tensor([[3.0], [4.0]])
if product.item() != 11.0:
    raise RuntimeError("Unexpected PyTorch matrix product")

# The original Margaret environment includes Intel TBB. Force Numba to use
# that conditional backend in a fresh process so an unresolved libtbb cannot
# be hidden by its normal OpenMP fallback.
tbb_check = subprocess.run(
    [
        sys.executable,
        "-c",
        "import numba,numpy as np; from numba import njit,prange; "
        "exec('@njit(parallel=True)\\ndef f(x):\\n s=0.0\\n for i in prange(x.size): s+=x[i]\\n return s',globals()); "
        "assert f(np.arange(10000,dtype=float)) == 49995000.0; "
        "assert numba.threading_layer() == 'tbb'",
    ],
    env={**os.environ, "NUMBA_THREADING_LAYER": "tbb"},
    capture_output=True,
    text=True,
)
if tbb_check.returncode != 0:
    raise RuntimeError(
        "Numba TBB backend failed:\n" + tbb_check.stdout + tbb_check.stderr
    )

edge_r_version = str(r('as.character(utils::packageVersion("edgeR"))')[0])
statmod_version = str(r('as.character(utils::packageVersion("statmod"))')[0])
r('stopifnot(sum(edgeR::cpm(matrix(c(1L, 2L, 3L, 4L), nrow=2L))) > 0)')

# Milo always requests edgeR's robust quasi-likelihood fit. statmod is a
# conditional limma dependency on this path, so imports and basic CPM checks
# cannot detect its absence. Exercise the complete neighbourhood DA route.
rng = np.random.default_rng(2504)
n_samples = 12
cells_per_sample = 20
samples = np.repeat(
    [f"sample_{index:02d}" for index in range(n_samples)],
    cells_per_sample,
)
conditions = np.repeat(["control"] * 6 + ["treated"] * 6, cells_per_sample)
milo_counts = rng.negative_binomial(
    12,
    0.55,
    size=(n_samples * cells_per_sample, 80),
).astype(np.float32)
milo_counts[conditions == "treated", :8] += rng.poisson(
    3,
    size=((conditions == "treated").sum(), 8),
)
milo_data = anndata.AnnData(
    X=milo_counts,
    obs=pd.DataFrame(
        {
            "sample": pd.Categorical(samples),
            "condition": pd.Categorical(
                conditions,
                categories=["control", "treated"],
            ),
        },
        index=[f"cell_{index:03d}" for index in range(len(samples))],
    ),
)
sc.pp.normalize_total(milo_data, target_sum=10_000)
sc.pp.log1p(milo_data)
sc.pp.pca(milo_data, n_comps=20)
sc.pp.neighbors(milo_data, n_neighbors=15, n_pcs=20)
milo.make_nhoods(milo_data, prop=0.2)
milo.count_nhoods(milo_data, sample_col="sample")
milo.DA_nhoods(milo_data, design="~condition")
milo_results = milo_data.uns["nhood_adata"].obs
required_milo_columns = {"logFC", "PValue", "FDR", "SpatialFDR"}
missing_milo_columns = required_milo_columns.difference(milo_results.columns)
if milo_results.empty or missing_milo_columns:
    raise RuntimeError(
        "Invalid Milo DA result; missing columns: "
        f"{sorted(missing_milo_columns)}"
    )

with tempfile.TemporaryDirectory() as directory:
    roundtrip_path = Path(directory) / "marg-python-smoke.h5ad"
    mini.write_h5ad(roundtrip_path)
    roundtrip = anndata.read_h5ad(roundtrip_path)
    if roundtrip.shape != mini.shape:
        raise RuntimeError("Synthetic H5AD round-trip changed the matrix shape")

print(
    "PASS marg-python",
    f"python={sys.version.split()[0]}",
    f"scanpy={package_version('scanpy')}",
    f"anndata={package_version('anndata')}",
    f"bbknn={package_version('bbknn')}",
    f"torch={package_version('torch')}",
    f"scvi-tools={package_version('scvi-tools')}",
    f"scArches={package_version('scArches')}",
    f"scib-metrics={package_version('scib-metrics')}",
    f"scikit-misc={package_version('scikit-misc')}",
    f"tbb={package_version('tbb')}",
    f"edgeR={edge_r_version}",
    f"statmod={statmod_version}",
    f"milo_nhoods={milo_results.shape[0]}",
    f"demo_shape={demo_shape}",
)
