#!/usr/bin/env python3
"""Minimal runtime check for the Margaret Python image."""

import argparse
import importlib.metadata
import sys
import tempfile
from pathlib import Path

import anndata
import bbknn
import numpy as np
import scanpy as sc
import scarches
import scvi
import torch
from rpy2.robjects import r


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

product = torch.tensor([[1.0, 2.0]]) @ torch.tensor([[3.0], [4.0]])
if product.item() != 11.0:
    raise RuntimeError("Unexpected PyTorch matrix product")

edge_r_version = str(r('as.character(utils::packageVersion("edgeR"))')[0])
r('stopifnot(sum(edgeR::cpm(matrix(c(1L, 2L, 3L, 4L), nrow=2L))) > 0)')

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
    f"edgeR={edge_r_version}",
    f"demo_shape={demo_shape}",
)
