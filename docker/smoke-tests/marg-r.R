#!/usr/bin/env Rscript

# Minimal runtime check for the Margaret R image.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: marg-r.R /path/to/Atlas_DEMO.h5ad", call. = FALSE)
}

demo_path <- normalizePath(args[[1L]], mustWork = TRUE)
expected_r <- "4.1.3"
expected_python <- "3.8.10"
actual_r <- as.character(getRversion())
if (!identical(actual_r, expected_r)) {
  stop(sprintf("Expected R %s, found R %s", expected_r, actual_r), call. = FALSE)
}

packages <- c(
  "Cairo",
  "batchelor",
  "DoubletFinder",
  "DropletUtils",
  "EnrichmentBrowser",
  "GSEABase",
  "Matrix",
  "SMFilter",
  "STACAS",
  "Seurat",
  "SeuratDisk",
  "SeuratWrappers",
  "ShinyCell",
  "SingleCellExperiment",
  "VennDiagram",
  "clusterProfiler",
  "cowplot",
  "data.table",
  "dittoSeq",
  "dplyr",
  "edgeR",
  "enrichplot",
  "future",
  "ggplot2",
  "ggrastr",
  "grid",
  "gtools",
  "harmony",
  "kableExtra",
  "org.Hs.eg.db",
  "patchwork",
  "readr",
  "reticulate",
  "rjson",
  "rliger",
  "scales",
  "scater",
  "scattermore",
  "sceasy",
  "tidyr",
  "tidyverse",
  "zellkonverter"
)
for (package in packages) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(sprintf("Required package is unavailable: %s", package), call. = FALSE)
  }
}

json_value <- rjson::fromJSON('{"functional":true,"count":2}')
if (!isTRUE(json_value$functional) || !identical(as.integer(json_value$count), 2L)) {
  stop("rjson did not parse the synthetic JSON payload", call. = FALSE)
}

set.seed(2504L)
counts <- matrix(
  stats::rpois(50L * 30L, lambda = 5),
  nrow = 50L,
  dimnames = list(paste0("gene", 1:50), paste0("cell", 1:30))
)
seurat_object <- Seurat::CreateSeuratObject(counts = counts)
seurat_object <- Seurat::NormalizeData(seurat_object, verbose = FALSE)
if (!"RNA" %in% names(seurat_object@assays)) {
  stop("Seurat did not create the expected RNA assay", call. = FALSE)
}

raster_plot <- Seurat::VlnPlot(
  seurat_object,
  features = "gene1",
  raster = TRUE,
  pt.size = 0.1
)
if (!inherits(raster_plot, "ggplot")) {
  stop("Seurat VlnPlot did not exercise the ggrastr raster path", call. = FALSE)
}

batch_a <- Seurat::CreateSeuratObject(counts = counts)
batch_b_counts <- matrix(
  stats::rpois(50L * 30L, lambda = 6),
  nrow = 50L,
  dimnames = list(paste0("gene", 1:50), paste0("cell", 31:60))
)
batch_b <- Seurat::CreateSeuratObject(counts = batch_b_counts)
batch_a <- Seurat::NormalizeData(batch_a, verbose = FALSE)
batch_b <- Seurat::NormalizeData(batch_b, verbose = FALSE)
fast_mnn <- SeuratWrappers::RunFastMNN(
  object.list = list(batch_a, batch_b),
  features = rownames(counts),
  k = 5L,
  d = 10L
)
if (!"mnn" %in% names(fast_mnn@reductions)) {
  stop("RunFastMNN did not create the expected mnn reduction", call. = FALSE)
}

sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))
sce <- scater::logNormCounts(sce)
if (!"logcounts" %in% SummarizedExperiment::assayNames(sce)) {
  stop("scater did not create a logcounts assay", call. = FALSE)
}

dge <- edgeR::DGEList(counts = counts)
dge <- edgeR::calcNormFactors(dge)
if (!all(is.finite(dge$samples$norm.factors))) {
  stop("edgeR returned invalid normalization factors", call. = FALSE)
}

anndata <- reticulate::import("anndata", convert = FALSE)
python_version <- reticulate::py_eval(
  "__import__('platform').python_version()",
  convert = TRUE
)
if (!identical(python_version, expected_python)) {
  stop(
    sprintf("Expected reticulate Python %s, found %s", expected_python, python_version),
    call. = FALSE
  )
}
demo <- anndata$read_h5ad(demo_path, backed = "r")
demo_shape <- as.integer(unlist(reticulate::py_to_r(demo$shape), use.names = FALSE))
invisible(demo$file$close())
if (length(demo_shape) != 2L || any(demo_shape == 0L)) {
  stop("Atlas_DEMO.h5ad is empty or has an invalid shape", call. = FALSE)
}

versions <- vapply(
  packages,
  function(package) as.character(utils::packageVersion(package)),
  character(1L)
)
cat(
  "PASS marg-r",
  sprintf("R=%s", actual_r),
  sprintf("reticulate_python=%s", python_version),
  sprintf("demo_shape=(%s)", paste(demo_shape, collapse = ", ")),
  paste(sprintf("%s=%s", names(versions), versions), collapse = " "),
  "\n"
)
