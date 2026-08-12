#!/usr/bin/env python3
"""Minimal runtime check for the Rocinante Python image."""

import argparse
import importlib.metadata
import sys
import tempfile
from pathlib import Path

import anndata
import matplotlib
import numba
import numpy as np
import palantir
import pandas
import pydeseq2
import scanpy as sc
import scipy
import seaborn
import statsmodels
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats


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

# Exercise the same PyDESeq2 fit and contrast path used by the Rocinante
# differential-expression scripts, not just its top-level imports.
rng = np.random.default_rng(2504)
counts = pandas.DataFrame(
    rng.negative_binomial(20, 0.5, size=(8, 30)),
    index=[f"sample_{index}" for index in range(8)],
    columns=[f"gene_{index}" for index in range(30)],
)
counts.loc[counts.index[4:], "gene_0"] += 30
metadata = pandas.DataFrame(
    {"condition": ["control"] * 4 + ["treated"] * 4},
    index=counts.index,
)
dds = DeseqDataSet(
    counts=counts,
    metadata=metadata,
    design_factors="condition",
    refit_cooks=False,
    quiet=True,
)
dds.deseq2()
stats = DeseqStats(
    dds,
    contrast=("condition", "treated", "control"),
    quiet=True,
)
stats.summary()
if stats.results_df.shape[0] != counts.shape[1]:
    raise RuntimeError("PyDESeq2 returned an unexpected number of genes")
if not np.isfinite(stats.results_df["baseMean"]).all():
    raise RuntimeError("PyDESeq2 returned invalid base means")

with tempfile.TemporaryDirectory() as directory:
    roundtrip_path = Path(directory) / "roci-python-smoke.h5ad"
    mini.write_h5ad(roundtrip_path)
    roundtrip = anndata.read_h5ad(roundtrip_path)
    if roundtrip.shape != mini.shape:
        raise RuntimeError("Synthetic H5AD round-trip changed the matrix shape")

print(
    "PASS roci-python",
    f"python={sys.version.split()[0]}",
    f"scanpy={package_version('scanpy')}",
    f"anndata={package_version('anndata')}",
    f"palantir={package_version('palantir')}",
    f"pydeseq2={package_version('pydeseq2')}",
    f"demo_shape={demo_shape}",
)
