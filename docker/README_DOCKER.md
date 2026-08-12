# Functional CAR-T Atlas Docker environments

This repository contains four tested Docker environments for the analysis
scripts under `Article/`.

| Image | Main runtime | Auxiliary bridge | Validated size |
|---|---|---|---:|
| `functional-cart-atlas-python-roci:0.2.0` | Python 3.8.10 | — | 1,043,280,309 bytes |
| `functional-cart-atlas-r-roci:0.1.0` | R 4.5.1 / Bioconductor 3.21 | Python 3.12.3 | 3,583,039,602 bytes |
| `functional-cart-atlas-python-marg:0.2.0` | Python 3.8.10 | R 4.1.3 for rpy2/edgeR | 6,577,596,013 bytes |
| `functional-cart-atlas-r-marg:0.1.0` | R 4.1.3 / Bioconductor 3.14 | Python 3.8.10 | 3,185,488,204 bytes |

`roci` and `marg` identify the original Rocinante and Margaret workstation
environments. The image tags are independent of the machine on which they are
run.

The Python `0.2.0` images include additional dependencies required by the
machine-specific analysis scripts and identified during a complete dependency
audit.

## Files required to build

The Dockerfiles are under `docker/roci/` and `docker/marg/`. Their exact tested
runtime closures are under `environments/roci/Docker/` and
`environments/marg/Docker/`:

- `requirements.runtime.txt`: Python packages used by the machine-specific
  scripts and their mandatory runtime dependencies.
- `renv.runtime.lock`: R packages used by the machine-specific R scripts.
- `requirements.r-reticulate.txt`: the small Python/anndata bridge used by each
  R image.
- `environments/marg/Docker/renv.python-bridge.lock`: the R/edgeR bridge used
  through rpy2 by the Margaret Python image.

The original workstation snapshots are retained separately under
`environments/*/Original/` for provenance. They are not consumed by the final
Dockerfiles.

## Pull published images

The tested images are publicly available from GHCR:

```bash
docker --context rootless pull \
  ghcr.io/ml4bm-lab/functional-cart-atlas-python-roci:0.2.0
docker --context rootless pull \
  ghcr.io/ml4bm-lab/functional-cart-atlas-r-roci:0.1.0
docker --context rootless pull \
  ghcr.io/ml4bm-lab/functional-cart-atlas-python-marg:0.2.0
docker --context rootless pull \
  ghcr.io/ml4bm-lab/functional-cart-atlas-r-marg:0.1.0
```

The examples below use the short local tags produced by the build commands.
When using a pulled image, replace its short tag with the corresponding full
`ghcr.io/ml4bm-lab/...` reference shown above.

## Build

Run the builds from the repository root. On a machine using the personal
rootless context, the four direct commands are:

```bash
docker --context rootless build --platform linux/amd64 \
  --file docker/roci/Dockerfile.python \
  --tag functional-cart-atlas-python-roci:0.2.0 \
  .

docker --context rootless build --platform linux/amd64 \
  --file docker/roci/Dockerfile.r \
  --tag functional-cart-atlas-r-roci:0.1.0 \
  .

docker --context rootless build --platform linux/amd64 \
  --file docker/marg/Dockerfile.python \
  --tag functional-cart-atlas-python-marg:0.2.0 \
  .

docker --context rootless build --platform linux/amd64 \
  --file docker/marg/Dockerfile.r \
  --tag functional-cart-atlas-r-marg:0.1.0 \
  .
```

Omit `--context rootless` on a machine that intentionally uses its default
Docker daemon. Change `--platform` only when a different target architecture is
required.

## Run with project files

Host files are made visible with a bind mount. Mounting the project at the same
absolute path is convenient for scripts that still use absolute paths:

```bash
ATLAS=/path/on/the/host/Atlas

docker --context rootless run --rm -it \
  --mount type=bind,source="$ATLAS",target="$ATLAS" \
  --workdir "$ATLAS" \
  functional-cart-atlas-python-marg:0.2.0 \
  bash
```

The mount above is writable. Add `,readonly` to protect it, or mount a separate
output directory with write access.

Initial runtime checks:

```bash
docker --context rootless run --rm \
  functional-cart-atlas-python-roci:0.2.0 python3 --version
docker --context rootless run --rm \
  functional-cart-atlas-python-marg:0.2.0 python3 --version
docker --context rootless run --rm \
  functional-cart-atlas-r-roci:0.1.0 R --version
docker --context rootless run --rm \
  functional-cart-atlas-r-marg:0.1.0 R --version
```

## Local smoke tests

The four smoke tests use `Atlas_DEMO.h5ad` read-only and write temporary output
only inside the disposable containers. Run them from the repository root:

```bash
docker --context rootless run --rm \
  --mount type=bind,source="$PWD",target=/repo,readonly \
  functional-cart-atlas-python-roci:0.2.0 \
  python3 /repo/docker/smoke-tests/roci-python.py /repo/Atlas_DEMO.h5ad

docker --context rootless run --rm \
  --mount type=bind,source="$PWD",target=/repo,readonly \
  functional-cart-atlas-r-roci:0.1.0 \
  Rscript /repo/docker/smoke-tests/roci-r.R /repo/Atlas_DEMO.h5ad

docker --context rootless run --rm \
  --mount type=bind,source="$PWD",target=/repo,readonly \
  functional-cart-atlas-python-marg:0.2.0 \
  python3 /repo/docker/smoke-tests/marg-python.py /repo/Atlas_DEMO.h5ad

docker --context rootless run --rm \
  --mount type=bind,source="$PWD",target=/repo,readonly \
  functional-cart-atlas-r-marg:0.1.0 \
  Rscript /repo/docker/smoke-tests/marg-r.R /repo/Atlas_DEMO.h5ad
```

Each command exits with a non-zero status if a required runtime, package,
R/Python bridge, synthetic operation or demo-file read fails. Omit
`--context rootless` when intentionally using the default Docker daemon.

## Transfer to another machine

For a one-off transfer, save and compress an image:

```bash
docker --context rootless save functional-cart-atlas-python-marg:0.2.0 \
  | zstd -T0 -3 -o functional-cart-atlas-python-marg-0.2.0.tar.zst
```

After transferring the file, load it into the intended Docker daemon:

```bash
zstd -dc functional-cart-atlas-python-marg-0.2.0.tar.zst \
  | docker load
```

Add `--context rootless` to the `docker load` command when the destination also
uses a personal rootless daemon.

## Validation summary

- Rocinante Python: 53 exact distributions; all 14 marked Python scripts parsed
  and synthetic Scanpy, H5AD, plotting, Palantir and PyDESeq2 workflows passed.
- Rocinante R: 312 exact R packages; all 27 marked R scripts parsed and their
  principal Bioconductor, plotting and reticulate/anndata paths passed.
- Margaret Python: 134 exact distributions; all 37 marked Python scripts parsed
  and Scanpy Seurat-v3 HVG, scIB, BBKNN, PyTorch, scVI/scArches, H5AD and
  rpy2/edgeR paths passed.
- Margaret R: 309 exact R packages and 12 Python bridge distributions; all 23
  marked R scripts parsed and Seurat, SeuratDisk, scater, DropletUtils, edgeR,
  enrichment, graphics and reticulate/anndata paths passed.

The final images were built and validated without modifying files under
`Article/`.
