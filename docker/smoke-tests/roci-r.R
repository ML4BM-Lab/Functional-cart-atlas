#!/usr/bin/env Rscript

# Minimal runtime check for the Rocinante R image.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: roci-r.R /path/to/Atlas_DEMO.h5ad", call. = FALSE)
}

demo_path <- normalizePath(args[[1L]], mustWork = TRUE)
expected_r <- "4.5.1"
expected_python <- "3.12.3"
actual_r <- as.character(getRversion())
if (!identical(actual_r, expected_r)) {
  stop(sprintf("Expected R %s, found R %s", expected_r, actual_r), call. = FALSE)
}

packages <- c(
  "AUCell",
  "BiocParallel",
  "Cairo",
  "EnrichmentBrowser",
  "GSEABase",
  "SMFilter",
  "SingleCellExperiment",
  "coin",
  "cowplot",
  "dplyr",
  "dreamlet",
  "edgeR",
  "ggplot2",
  "ggsignif",
  "kableExtra",
  "patchwork",
  "reshape2",
  "reticulate",
  "rjson",
  "scater",
  "scattermore",
  "see",
  "statmod",
  "tidyverse",
  "variancePartition",
  "zellkonverter",
  "zenith"
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

counts <- matrix(
  c(1L, 0L, 3L, 2L, 4L, 1L, 2L, 3L, 1L, 1L, 2L, 4L),
  nrow = 3L,
  dimnames = list(paste0("gene", 1:3), paste0("cell", 1:4))
)
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

set.seed(2504L)
da_counts <- matrix(
  stats::rpois(50L * 12L, lambda = 8),
  nrow = 50L,
  dimnames = list(paste0("gene", 1:50), paste0("sample", 1:12))
)
group <- factor(rep(c("control", "treated"), each = 6L))
design <- stats::model.matrix(~group)
da_dge <- edgeR::DGEList(counts = da_counts)
da_dge <- edgeR::calcNormFactors(da_dge)
da_dge <- edgeR::estimateDisp(da_dge, design)
robust_fit <- edgeR::glmQLFit(da_dge, design, robust = TRUE)
robust_test <- edgeR::glmQLFTest(robust_fit, coef = 2L)
robust_table <- edgeR::topTags(robust_test, n = Inf)$table
if (nrow(robust_table) != nrow(da_counts) || any(!is.finite(robust_table$PValue))) {
  stop("edgeR robust quasi-likelihood path returned invalid results", call. = FALSE)
}

anndata <- reticulate::import("anndata")
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
demo_shape_value <- demo$shape
if (inherits(demo_shape_value, "python.builtin.object")) {
  demo_shape_value <- reticulate::py_to_r(demo_shape_value)
}
demo_shape <- as.integer(unlist(demo_shape_value, use.names = FALSE))
invisible(demo$file$close())
if (length(demo_shape) != 2L || any(demo_shape == 0L)) {
  stop("Atlas_DEMO.h5ad is empty or has an invalid shape", call. = FALSE)
}
demo_memory <- anndata$read_h5ad(demo_path)
demo_sce <- zellkonverter::AnnData2SCE(
  demo_memory,
  "counts",
  uns = FALSE,
  obsm = FALSE,
  obsp = FALSE
)
if (!identical(rev(dim(demo_sce)), demo_shape)) {
  stop("AnnData2SCE changed the demo matrix shape", call. = FALSE)
}

versions <- vapply(
  packages,
  function(package) as.character(utils::packageVersion(package)),
  character(1L)
)
cat(
  "PASS roci-r",
  sprintf("R=%s", actual_r),
  sprintf("reticulate_python=%s", python_version),
  sprintf("demo_shape=(%s)", paste(demo_shape, collapse = ", ")),
  paste(sprintf("%s=%s", names(versions), versions), collapse = " "),
  "\n"
)
